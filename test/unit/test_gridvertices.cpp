// Shared cell corners, rotation, and the precision of written values
#include <cmath>
#include <fstream>
#include <map>
#include <set>
#include <sstream>
#include <string>
#include <utility>
#include <vector>

#include <gtest/gtest.h>

#include "GridToShp.h"
#include "GridToVTK.h"
#include "GridVertices.h"
#include "RoundTrip.h"
#include "def_exporter.h"
#include "def_parser.h"
#include "def_struct.h"
#include "grid_helpers.h"

using namespace cusg;
using gridgen_test::make_mfgrid;
using gridgen_test::make_qtree;

namespace
{

const double PI_DOUBLE = 3.14159265358979323846;
const double XOFF = 543750.0;
const double YOFF = 2792250.0;

// a quadtree with cell sizes that are not exact in binary, at projected
// coordinate offsets, refined to three levels around two points
QuadTree3D * make_refined(double angle)
{
    ModflowGrid * mfgrid = make_mfgrid(1, 6, 6, 812.7, 793.3, 10.0, XOFF, YOFF, angle);
    return make_qtree(
        mfgrid, {{XOFF + 1700.0, YOFF + 2500.0, 3, 0}, {XOFF + 3900.0, YOFF + 1100.0, 2, 0}}, true);
}

// polygon vertices of a shapefile
std::vector<std::pair<double, double>> read_shp_vertices(const std::string & shp)
{
    std::vector<std::pair<double, double>> xy;
    SHPHandle h = SHPOpen(shp.c_str(), "rb");
    if (h == NULL) return xy;
    int n = 0;
    SHPGetInfo(h, &n, NULL, NULL, NULL);
    for (int i = 0; i < n; i++)
    {
        SHPObject * obj = SHPReadObject(h, i);
        for (int v = 0; v < obj->nVertices; v++)
            xy.push_back(std::make_pair(obj->padfX[v], obj->padfY[v]));
        SHPDestroyObject(obj);
    }
    SHPClose(h);
    return xy;
}

// the distinct coordinates written for each corner, keyed by the corner
// rounded to a millimeter
std::map<std::pair<long long, long long>, std::set<std::pair<double, double>>>
group_corners(const std::vector<std::pair<double, double>> & xy)
{
    std::map<std::pair<long long, long long>, std::set<std::pair<double, double>>> groups;
    for (const std::pair<double, double> & p : xy)
        groups[std::make_pair(std::llround(p.first * 1e3), std::llround(p.second * 1e3))].insert(p);
    return groups;
}

void expect_one_coordinate_per_corner(double angle, const std::string & name)
{
    QuadTree3D * qtree = make_refined(angle);
    std::string shp = ::testing::TempDir() + name + ".shp";
    writeQuadtree2Shp(qtree, shp, "polygon");

    std::vector<std::pair<double, double>> xy = read_shp_vertices(shp);
    ASSERT_FALSE(xy.empty());
    for (const auto & g : group_corners(xy))
        EXPECT_EQ(g.second.size(), 1u)
            << "corner near (" << g.first.first / 1e3 << ", " << g.first.second / 1e3 << ")";
    delete qtree;
}

// the corner of a leaf box nearest a point
int nearest_corner(const GridVertices & gv, const Box * b, double x, double y)
{
    const int * ids = gv.corners(b);
    int best = ids[0];
    double dmin = 1e300;
    for (int i = 0; i < 4; i++)
    {
        const Point2d & p = gv.vertices()[ids[i]];
        double d = std::hypot(p[0] - x, p[1] - y);
        if (d < dmin)
        {
            dmin = d;
            best = ids[i];
        }
    }
    return best;
}

// x and y of the points of a vtu file
std::vector<std::pair<double, double>> read_vtu_xy(const std::string & path)
{
    std::ifstream fin(path.c_str());
    std::string text((std::istreambuf_iterator<char>(fin)), std::istreambuf_iterator<char>());
    std::vector<std::pair<double, double>> xy;
    size_t start = text.find("<Points>");
    if (start == std::string::npos) return xy;
    start = text.find(">", text.find("<DataArray", start)) + 1;
    std::istringstream ss(text.substr(start, text.find("</DataArray>", start) - start));
    double x, y, z;
    while (ss >> x >> y >> z) xy.push_back(std::make_pair(x, y));
    return xy;
}

modflow_grid_raw_data * parse_grid(DefinitionParser & parser, const std::string & block)
{
    std::istringstream dfn(block);
    if (!parser.parse(dfn)) return NULL;
    modflow_grid_raw_data * data =
        dynamic_cast<modflow_grid_raw_data *>(parser.getBlockData("mfg"));
    if (data == NULL || !data->initialize(parser.getBlockData())) return NULL;
    return data;
}

// a grid one cell wide, 110 km long, rotated 17.3 degrees
const char * LONG_GRID = "BEGIN MODFLOW_GRID mfg\n"
                         "  ROTATION_ANGLE = 17.3\n"
                         "  X_OFFSET = 543750\n"
                         "  Y_OFFSET = 2792250\n"
                         "  NLAY = 1\n"
                         "  NROW = 1\n"
                         "  NCOL = 1\n"
                         "  DELR = CONSTANT 110000.\n"
                         "  DELC = CONSTANT 1.\n"
                         "  TOP = CONSTANT 0.\n"
                         "  BOTTOM LAYER 1 = CONSTANT -1.\n"
                         "END MODFLOW_GRID\n";

} // namespace

TEST(SharedVertices, ShapefileCornersAreIdenticalUnrotated)
{
    expect_one_coordinate_per_corner(0.0, "gridgen_corners_unrotated");
}

TEST(SharedVertices, ShapefileCornersAreIdenticalRotated)
{
    expect_one_coordinate_per_corner(17.3 * PI_DOUBLE / 180.0, "gridgen_corners_rotated");
}

TEST(SharedVertices, AdjacentCellsShareVertexIndex)
{
    for (double angle : {0.0, 17.3 * PI_DOUBLE / 180.0})
    {
        QuadTree3D * qtree = make_refined(angle);
        GridVertices gv(qtree);

        // every corner of every leaf, by position, maps to one vertex index
        std::map<std::pair<long long, long long>, std::set<int>> groups;
        for (int n = 0; n < qtree->nodes(); n++)
        {
            Box * b = qtree->get_nodeobj(n);
            const int * ids = gv.corners(b);
            for (int i = 0; i < 4; i++)
            {
                const Point2d & p = gv.vertices()[ids[i]];
                groups[std::make_pair(std::llround(p[0] * 1e3), std::llround(p[1] * 1e3))].insert(
                    ids[i]);
            }
        }
        for (const auto & g : groups) EXPECT_EQ(g.second.size(), 1u);
        EXPECT_EQ(groups.size(), gv.vertices().size());

        // a coarse cell and the fine cells along its east edge share the
        // corner they meet at
        for (int n = 0; n < qtree->nodes(); n++)
        {
            Box * b = qtree->get_nodeobj(n);
            const int * ids = gv.corners(b);
            for (int i = 0; i < 4; i++)
            {
                const Point2d & p = gv.vertices()[ids[i]];
                for (int m = 0; m < qtree->nodes(); m++)
                {
                    Box * o = qtree->get_nodeobj(m);
                    int j = nearest_corner(gv, o, p[0], p[1]);
                    const Point2d & q = gv.vertices()[j];
                    if (std::hypot(q[0] - p[0], q[1] - p[1]) < 1e-3) EXPECT_EQ(j, ids[i]);
                }
            }
        }
        delete qtree;
    }
}

TEST(SharedVertices, QuarterTurnAboutLowerLeftCorner)
{
    // one 4 x 2 cell with its lower-left corner at (8, 19)
    ModflowGrid * mfgrid = make_mfgrid(1, 1, 1, 4.0, 2.0, 10.0, 8.0, 19.0, 0.5 * PI_DOUBLE);
    GridVertices gv(mfgrid);
    double x[4], y[4];
    gv.get_corners(mfgrid->nodegroup[0], x, y);

    // a quarter turn counterclockwise maps (x, y) to (-y, x)
    const double ex[4] = {6.0, 6.0, 8.0, 8.0};
    const double ey[4] = {19.0, 23.0, 23.0, 19.0};
    for (int i = 0; i < 4; i++)
    {
        EXPECT_NEAR(x[i], ex[i], 1e-12) << "corner " << i;
        EXPECT_NEAR(y[i], ey[i], 1e-12) << "corner " << i;
    }
    const Point2d & c = gv.center(mfgrid->nodegroup[0]);
    EXPECT_NEAR(c[0], 7.0, 1e-12);
    EXPECT_NEAR(c[1], 21.0, 1e-12);
    delete mfgrid;
}

TEST(SharedVertices, UnrotatedCornersAreEdges)
{
    QuadTree3D * qtree = make_refined(0.0);
    GridVertices gv(qtree);
    ModflowGrid * mfgrid = qtree->getModflowGrid();

    // the base grid corners are the offset plus the summed cell sizes
    double xe = 0.0;
    for (int c = 0; c < mfgrid->ncol; c++) xe += mfgrid->delr[c];
    double ye = 0.0;
    for (int r = 0; r < mfgrid->nrow; r++) ye += mfgrid->delc[r];
    Box * ne = qtree->nodegroup[mfgrid->get_nodeid(0, 0, mfgrid->ncol - 1)];
    while (!ne->isLeaf) ne = ne->pChildren[1];
    const Point2d & p = gv.vertices()[gv.corners(ne)[1]];
    EXPECT_EQ(p[0], XOFF + xe);
    EXPECT_EQ(p[1], YOFF + ye);
    delete qtree;
}

TEST(SharedVertices, RotationUsesExactPi)
{
    DefinitionParser parser;
    modflow_grid_raw_data * data = parse_grid(parser, LONG_GRID);
    ASSERT_NE(data, nullptr);
    EXPECT_DOUBLE_EQ(data->rotation_angle, 17.3 * PI_DOUBLE / 180.0);

    ModflowGrid * mfgrid = dynamic_cast<ModflowGrid *>(data->build());
    ASSERT_NE(mfgrid, nullptr);
    std::string shp = ::testing::TempDir() + "gridgen_long_rotated.shp";
    writeModflowGrid2Shp(mfgrid, shp, "polygon");
    std::vector<std::pair<double, double>> xy = read_shp_vertices(shp);
    ASSERT_GE(xy.size(), 4u);

    // the upper-right corner, 110 km from the origin, rotated with exact pi
    double t = 17.3 * PI_DOUBLE / 180.0;
    double ex = XOFF + 110000.0 * std::cos(t) - 1.0 * std::sin(t);
    double ey = YOFF + 110000.0 * std::sin(t) + 1.0 * std::cos(t);
    double dmin = 1e300;
    for (const std::pair<double, double> & p : xy)
        dmin = std::min(dmin, std::hypot(p.first - ex, p.second - ey));
    EXPECT_LT(dmin, 1e-6);
}

TEST(SharedVertices, RotationAngleRoundTrips)
{
    DefinitionParser parser;
    modflow_grid_raw_data * data = parse_grid(parser, LONG_GRID);
    ASSERT_NE(data, nullptr);

    std::ostringstream out;
    DefinitionExporter exporter;
    exporter.save(out, *data);
    EXPECT_NE(out.str().find("\tROTATION_ANGLE = 17.3\n"), std::string::npos) << out.str();
}

TEST(Precision, ArraysRoundTrip)
{
    double values[] = {0.1, 1.0 / 3.0, 2792250.123456789, 25.396875, 800.0};
    double * ptr = values;
    std::string path = ::testing::TempDir() + "gridgen_array.dat";
    DefinitionExporter exporter;
    exporter.save(path, ptr, 5);

    std::ifstream fin(path.c_str());
    std::string text((std::istreambuf_iterator<char>(fin)), std::istreambuf_iterator<char>());
    std::istringstream in(text);
    for (double v : values)
    {
        double r = 0.0;
        in >> r;
        EXPECT_EQ(r, v);
    }
    // and with no more digits than needed
    EXPECT_EQ(text.substr(0, 4), "0.1 ");
}

TEST(Precision, ConnectionLengthsAreDouble)
{
    // 812.7 / 32 has no float representation
    ModflowGrid * mfgrid = make_mfgrid(1, 1, 2, 812.7 / 16, 50.0);
    CSRData csr = mfgrid->get_usg_csr_data();
    int nconn = 0;
    for (int ipos = 0; ipos < csr.nja; ipos++)
    {
        if (csr.fldr[ipos] == 0) continue; // the diagonal
        EXPECT_EQ(csr.cl1[ipos], 812.7 / 32);
        EXPECT_EQ(csr.cl2[ipos], 812.7 / 32);
        nconn++;
    }
    EXPECT_EQ(nconn, 2);
    csr.destory();
    delete mfgrid;
}

TEST(SharedVertices, VaryingCellSizesRefinedAtGridCorner)
{
    // columns and rows of different sizes, the northeast cell refined to
    // three levels, and the grid rotated
    const double delr[] = {812.7, 400.3, 1001.9, 650.1};
    const double delc[] = {793.3, 1200.7, 333.1};
    std::string dir = ::testing::TempDir();
    {
        std::ofstream fr((dir + "gridgen_delr.dat").c_str());
        for (double v : delr) fr << round_trip(v) << " ";
        std::ofstream fc((dir + "gridgen_delc.dat").c_str());
        for (double v : delc) fc << round_trip(v) << " ";
    }
    std::string block = "BEGIN MODFLOW_GRID mfg\n"
                        "  ROTATION_ANGLE = 17.3\n"
                        "  X_OFFSET = 543750\n"
                        "  Y_OFFSET = 2792250\n"
                        "  NLAY = 1\n"
                        "  NROW = 3\n"
                        "  NCOL = 4\n"
                        "  DELR = OPEN/CLOSE " +
                        dir +
                        "gridgen_delr.dat\n"
                        "  DELC = OPEN/CLOSE " +
                        dir +
                        "gridgen_delc.dat\n"
                        "  TOP = CONSTANT 0.\n"
                        "  BOTTOM LAYER 1 = CONSTANT -1.\n"
                        "END MODFLOW_GRID\n";
    DefinitionParser parser;
    modflow_grid_raw_data * data = parse_grid(parser, block);
    ASSERT_NE(data, nullptr);
    ModflowGrid * mfgrid = dynamic_cast<ModflowGrid *>(data->build());
    ASSERT_NE(mfgrid, nullptr);
    for (int c = 0; c < 4; c++) ASSERT_EQ(mfgrid->delr[c], delr[c]);
    for (int r = 0; r < 3; r++) ASSERT_EQ(mfgrid->delc[r], delc[r]);

    // refinement points in unrotated grid coordinates
    double width = 0.0, height = 0.0;
    for (double v : delr) width += v;
    for (double v : delc) height += v;
    QuadTree3D * qtree = make_qtree(mfgrid, {{XOFF + width - 1.0, YOFF + height - 1.0, 3, 0},
                                             {XOFF + 900.0, YOFF + 400.0, 2, 0}});
    GridVertices gv(qtree);

    // the northeast leaf's upper-right corner is the far grid corner, with
    // the edges summed west to east and south to north; compared to 4 ulps,
    // since a compiler may fuse the multiply-adds of either expression
    double xl = 0.0, yl = 0.0;
    for (int c = 0; c < 4; c++) xl += delr[c];
    for (int r = 2; r >= 0; r--) yl += delc[r];
    double t = 17.3 * PI_DOUBLE / 180.0;
    Box * ne = qtree->nodegroup[mfgrid->get_nodeid(0, 0, 3)];
    ASSERT_FALSE(ne->isLeaf);
    while (!ne->isLeaf) ne = ne->pChildren[1];
    const Point2d & p = gv.vertices()[gv.corners(ne)[1]];
    EXPECT_DOUBLE_EQ(p[0], XOFF + (xl * std::cos(t) - yl * std::sin(t)));
    EXPECT_DOUBLE_EQ(p[1], YOFF + (xl * std::sin(t) + yl * std::cos(t)));

    // and every corner is one vertex
    std::map<std::pair<long long, long long>, std::set<int>> groups;
    for (int n = 0; n < qtree->nodes(); n++)
    {
        const int * ids = gv.corners(qtree->get_nodeobj(n));
        for (int i = 0; i < 4; i++)
        {
            const Point2d & q = gv.vertices()[ids[i]];
            groups[std::make_pair(std::llround(q[0] * 1e3), std::llround(q[1] * 1e3))].insert(
                ids[i]);
        }
    }
    for (const auto & g : groups) EXPECT_EQ(g.second.size(), 1u);
    EXPECT_EQ(groups.size(), gv.vertices().size());
    delete qtree;
}

TEST(SharedVertices, LineShapefileEndsAreGridCorners)
{
    QuadTree3D * qtree = make_refined(17.3 * PI_DOUBLE / 180.0);
    GridVertices gv(qtree);
    std::set<std::pair<double, double>> corners;
    for (const Point2d & p : gv.vertices()) corners.insert(std::make_pair(p[0], p[1]));

    std::string shp = ::testing::TempDir() + "gridgen_lines.shp";
    writeQuadtree2Shp(qtree, shp, "line");
    std::vector<std::pair<double, double>> xy = read_shp_vertices(shp);
    ASSERT_FALSE(xy.empty());
    for (const std::pair<double, double> & p : xy)
        EXPECT_EQ(corners.count(p), 1u) << "line end (" << p.first << ", " << p.second << ")";
    delete qtree;
}

TEST(Precision, PassThroughAreasAreDouble)
{
    // the center cell of layer 2 is inactive, so with pass-through the
    // center cells of layers 1 and 3 connect over the full cell area
    QuadTree3D * qtree = make_qtree(make_mfgrid(3, 3, 3, 812.7 / 16, 793.3 / 16), {});
    ModflowGrid * mfgrid = qtree->getModflowGrid();
    qtree->nodegroup[mfgrid->get_nodeid(1, 1, 1)]->active = false;
    qtree->number_nodes();
    qtree->refresh();
    qtree->vertical_pass_through = true;

    CSRData csr = qtree->get_usg_csr_data();
    int nvert = 0;
    for (int ipos = 0; ipos < csr.nja; ipos++)
    {
        if (std::abs(csr.fldr[ipos]) != 3) continue;
        EXPECT_DOUBLE_EQ(csr.fahl[ipos], (812.7 / 16) * (793.3 / 16));
        nvert++;
    }
    EXPECT_GT(nvert, 0);
    csr.destory();
    delete qtree;
}

TEST(Precision, VtkPointsAreGridCorners)
{
    for (double angle : {0.0, 17.3 * PI_DOUBLE / 180.0})
    {
        QuadTree3D * qtree = make_refined(angle);
        GridVertices gv(qtree);
        std::set<std::pair<double, double>> corners;
        for (const Point2d & p : gv.vertices()) corners.insert(std::make_pair(p[0], p[1]));

        // every x and y written, with and without shared vertices, is a grid
        // corner
        std::string name = ::testing::TempDir() + "gridgen_vtk_precision";
        writeQuadtreeGrid2VTK(qtree, name);
        std::vector<std::pair<double, double>> xy = read_vtu_xy(name + ".vtu");
        EXPECT_EQ(xy.size(), 8u * qtree->nodes());
        for (const std::pair<double, double> & p : xy)
            EXPECT_EQ(corners.count(p), 1u) << "point (" << p.first << ", " << p.second << ")";

        writeQuadtreeGrid2VTK_shareV(qtree, name + "_sv");
        xy = read_vtu_xy(name + "_sv.vtu");
        EXPECT_FALSE(xy.empty());
        for (const std::pair<double, double> & p : xy)
            EXPECT_EQ(corners.count(p), 1u)
                << "shared point (" << p.first << ", " << p.second << ")";
        delete qtree;
    }
}
