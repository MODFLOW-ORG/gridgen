// MODFLOW 6 DISV input of a grid with the same cells in every layer
#include <cmath>
#include <filesystem>
#include <fstream>
#include <set>
#include <sstream>
#include <string>
#include <tuple>
#include <vector>

#include <gtest/gtest.h>

#include "DisvExport.h"
#include "GhostNode.h"
#include "GridgenError.h"
#include "grid_helpers.h"

using namespace cusg;
using gridgen_test::make_mfgrid;
using gridgen_test::make_qtree;

namespace
{

// a 2 x 2 grid of 100 x 100 cells, nlay layers 10 thick, with the northwest
// cell of every layer divided into four
std::unique_ptr<QuadTree3D> make_northwest_refined(int nlay)
{
    std::vector<gridgen_test::RefinePoint> points;
    for (int k = 0; k < nlay; k++) points.push_back({50.0, 150.0, 1, k});
    std::unique_ptr<QuadTree3D> qtree = make_qtree(make_mfgrid(nlay, 2, 2, 100.0, 100.0), points);
    qtree->name = "qtg";
    return qtree;
}

// the leaf cell centered at (x, y) in layer
Box * cell_at(QuadTree3D * qtree, double x, double y, int layer)
{
    for (Box * box : qtree->node_group)
        if (box->is_leaf && box->layer == layer && box->x == x && box->y == y) return box;
    return nullptr;
}

// deactivate the cells centered at the points of layer, then renumber the
// cells and recompute their elevations
void deactivate(QuadTree3D * qtree, const std::vector<std::pair<double, double>> & centers,
                int layer)
{
    for (const auto & center : centers)
    {
        Box * box = cell_at(qtree, center.first, center.second, layer);
        ASSERT_NE(box, nullptr) << center.first << ", " << center.second;
        box->active = false;
    }
    qtree->number_nodes();
    qtree->refresh();
    qtree->get_top();
    qtree->get_bot();
    qtree->set_cell_elevations();
}

// the index of the cell of a layer centered at (x, y), or -1
int cell_index(const DisvData & data, double x, double y)
{
    for (size_t j = 0; j < data.cell2d.size(); j++)
        if (data.cell2d[j].x == x && data.cell2d[j].y == y) return (int)j;
    return -1;
}

// the area of a cell from its vertices, which are clockwise
double cell_area(const DisvData & data, const Cell2d & cell)
{
    double area = 0.0;
    size_t count = cell.vertices.size();
    for (size_t i = 0; i < count; i++)
    {
        const Point2d & a = data.vertices[cell.vertices[i]];
        const Point2d & b = data.vertices[cell.vertices[(i + 1) % count]];
        area += a[1] * b[0] - a[0] * b[1];
    }
    return area / 2.0;
}

std::vector<std::string> read_lines(const std::string & filename)
{
    std::ifstream file(filename);
    std::vector<std::string> lines;
    for (std::string line; std::getline(file, line);) lines.push_back(line);
    return lines;
}

} // namespace

TEST(DisvExport, ModflowGridCellsAreNumberedRowByRow)
{
    // 2 layers of 2 rows and 3 columns: the cells of a layer are numbered from
    // the north row by row, with the elevations of the grid and no ghost nodes
    std::unique_ptr<ModflowGrid> mfgrid = make_mfgrid(2, 2, 3, 100.0, 50.0);
    DisvData data = disv_data(mfgrid.get(), false);
    EXPECT_EQ(data.nlay, 2);
    EXPECT_EQ(data.ncpl, 6);
    EXPECT_EQ(data.vertices.size(), 12u);
    ASSERT_EQ(data.cell2d.size(), 6u);
    EXPECT_EQ(data.cell2d[0].x, 50.0);
    EXPECT_EQ(data.cell2d[0].y, 75.0);
    EXPECT_EQ(data.cell2d[1].x, 150.0);
    EXPECT_EQ(data.cell2d[3].y, 25.0);
    for (const Cell2d & cell : data.cell2d)
    {
        EXPECT_EQ(cell.vertices.size(), 4u);
        EXPECT_DOUBLE_EQ(cell_area(data, cell), 5000.0);
    }
    EXPECT_EQ(data.top, std::vector<double>(6, 0.0));
    ASSERT_EQ(data.botm.size(), 12u);
    EXPECT_EQ(data.botm[0], -10.0);
    EXPECT_EQ(data.botm[6], -20.0);
    EXPECT_EQ(data.idomain, std::vector<int>(12, 1));
    EXPECT_FALSE(data.has_ghost_nodes);
    EXPECT_TRUE(data.ghost_nodes.empty());
}

TEST(DisvExport, LayersShareTheCellsAndHangingVertices)
{
    // the northwest cell divided into four in both layers: 7 cells, whose areas
    // from their vertices fill the grid; the northeast and southwest cells list
    // the corner of the small cells on their edge
    std::unique_ptr<QuadTree3D> qtree = make_northwest_refined(2);
    DisvData data = disv_data(qtree.get(), false);
    EXPECT_EQ(data.ncpl, 7);
    EXPECT_EQ(data.idomain, std::vector<int>(14, 1));
    double total = 0.0;
    for (const Cell2d & cell : data.cell2d) total += cell_area(data, cell);
    EXPECT_DOUBLE_EQ(total, 40000.0);
    int northeast = cell_index(data, 150.0, 150.0);
    int southwest = cell_index(data, 50.0, 50.0);
    int southeast = cell_index(data, 150.0, 50.0);
    ASSERT_GE(northeast, 0);
    ASSERT_GE(southwest, 0);
    ASSERT_GE(southeast, 0);
    EXPECT_EQ(data.cell2d[northeast].vertices.size(), 5u);
    EXPECT_EQ(data.cell2d[southwest].vertices.size(), 5u);
    EXPECT_EQ(data.cell2d[southeast].vertices.size(), 4u);
    EXPECT_EQ(data.tops_not_aligned, 0);
}

TEST(DisvExport, CellsInactiveInOneLayerKeepTheirVertices)
{
    // the small cells inactive in layer 1 and active in layer 2: their corner is
    // still on the edges of the northeast and southwest cells, which MODFLOW 6
    // needs to connect those cells to them in layer 2; inactive cells are 0 in
    // IDOMAIN, or -1 with vertical pass-through
    std::unique_ptr<QuadTree3D> qtree = make_northwest_refined(2);
    deactivate(qtree.get(), {{25.0, 175.0}, {75.0, 175.0}, {25.0, 125.0}, {75.0, 125.0}}, 0);
    for (bool vertical_pass_through : {false, true})
    {
        SCOPED_TRACE(vertical_pass_through ? "pass-through" : "no pass-through");
        DisvData data = disv_data(qtree.get(), vertical_pass_through);
        ASSERT_EQ(data.ncpl, 7);
        EXPECT_EQ(data.cell2d[cell_index(data, 150.0, 150.0)].vertices.size(), 5u);
        EXPECT_EQ(data.cell2d[cell_index(data, 50.0, 50.0)].vertices.size(), 5u);
        int inactive = vertical_pass_through ? -1 : 0;
        for (double x : {25.0, 75.0})
        {
            for (double y : {125.0, 175.0})
            {
                int j = cell_index(data, x, y);
                ASSERT_GE(j, 0);
                EXPECT_EQ(data.idomain[j], inactive);
                EXPECT_EQ(data.idomain[data.ncpl + j], 1);
            }
        }
        int active = 0;
        for (int value : data.idomain) active += (value == 1);
        EXPECT_EQ(active, 10);
    }
}

TEST(DisvExport, InactiveCellAboveAnActiveCellEndsAtItsTop)
{
    // the northeast cell inactive in layer 1: its bottom is the top of the
    // northeast cell of layer 2, which DISV takes as that cell's top, and its
    // top is that of the grid; a top that differs from the bottom of an active
    // cell above is counted
    std::unique_ptr<QuadTree3D> qtree = make_northwest_refined(2);
    deactivate(qtree.get(), {{150.0, 150.0}}, 0);
    Box * below = cell_at(qtree.get(), 150.0, 150.0, 1);
    qtree->top[below->number] = -9.0;
    Box * southeast = cell_at(qtree.get(), 150.0, 50.0, 1);
    qtree->top[southeast->number] = -9.5;

    DisvData data = disv_data(qtree.get(), false);
    int northeast = cell_index(data, 150.0, 150.0);
    ASSERT_GE(northeast, 0);
    EXPECT_EQ(data.top[northeast], 0.0);
    EXPECT_EQ(data.botm[northeast], -9.0);
    EXPECT_EQ(data.botm[data.ncpl + northeast], -20.0);
    EXPECT_EQ(data.tops_not_aligned, 1);
}

TEST(DisvExport, InactiveCellsInConsecutiveLayers)
{
    // the northeast cell inactive in layers 1 and 2 of three: the bottom of the
    // one in layer 2 is the top of the active cell below it, and the bottom of
    // the one in layer 1 and the top are those of the base grid
    std::unique_ptr<QuadTree3D> qtree = make_northwest_refined(3);
    deactivate(qtree.get(), {{150.0, 150.0}}, 0);
    deactivate(qtree.get(), {{150.0, 150.0}}, 1);
    Box * below = cell_at(qtree.get(), 150.0, 150.0, 2);
    qtree->top[below->number] = -19.0;

    DisvData data = disv_data(qtree.get(), true);
    int northeast = cell_index(data, 150.0, 150.0);
    ASSERT_GE(northeast, 0);
    EXPECT_EQ(data.top[northeast], 0.0);
    EXPECT_EQ(data.botm[northeast], -10.0);
    EXPECT_EQ(data.botm[data.ncpl + northeast], -19.0);
    EXPECT_EQ(data.botm[2 * data.ncpl + northeast], -30.0);
    EXPECT_EQ(data.idomain[northeast], -1);
    EXPECT_EQ(data.idomain[data.ncpl + northeast], -1);
    EXPECT_EQ(data.idomain[2 * data.ncpl + northeast], 1);
    EXPECT_EQ(data.tops_not_aligned, 0);
}

TEST(DisvExport, RotatedGridVerticesAreRotatedAndClockwise)
{
    // the grid of LayersShareTheCellsAndHangingVertices rotated 30 degrees about
    // its origin: the cell centers are rotated, and the areas from the
    // vertices are still positive, clockwise, and fill the grid
    // the grid takes its rotation in radians
    double angle = std::acos(-1.0) / 6.0;
    std::vector<gridgen_test::RefinePoint> points = {{50.0, 150.0, 1, 0}, {50.0, 150.0, 1, 1}};
    std::unique_ptr<QuadTree3D> qtree =
        make_qtree(make_mfgrid(2, 2, 2, 100.0, 100.0, 10.0, 0.0, 0.0, angle), points);
    DisvData data = disv_data(qtree.get(), false);
    ASSERT_EQ(data.ncpl, 7);
    double total = 0.0;
    for (const Cell2d & cell : data.cell2d)
    {
        double area = cell_area(data, cell);
        EXPECT_GT(area, 0.0);
        total += area;
    }
    EXPECT_NEAR(total, 40000.0, 1.0e-6);
    // the southeast cell, centered at (150, 50) before rotation
    double x = 150.0 * std::cos(angle) - 50.0 * std::sin(angle);
    double y = 150.0 * std::sin(angle) + 50.0 * std::cos(angle);
    bool found = false;
    for (const Cell2d & cell : data.cell2d)
        found = found || (std::fabs(cell.x - x) < 1.0e-9 && std::fabs(cell.y - y) < 1.0e-9);
    EXPECT_TRUE(found) << x << ", " << y;
}

TEST(DisvExport, GhostNodesAreByLayerAndCellOfTheLayer)
{
    // the southwest cell inactive in layer 1, so that the cells of layer 2 are
    // numbered differently in DISU, and a small cell inactive in layer 2, next
    // to the ghost nodes of that layer; each ghost node is that of the DISU
    // grid, by layer and cell of the layer
    std::unique_ptr<QuadTree3D> qtree = make_northwest_refined(2);
    deactivate(qtree.get(), {{50.0, 50.0}}, 0);
    deactivate(qtree.get(), {{75.0, 175.0}}, 1);
    std::vector<GhostNode> expected_nodes;
    compute_ghost_nodes(qtree.get(), expected_nodes, false);
    DisvData data = disv_data(qtree.get(), false);
    EXPECT_TRUE(data.has_ghost_nodes);
    ASSERT_EQ(data.ghost_nodes.size(), expected_nodes.size());

    using Record = std::tuple<int, double, double, double, double, double, double, double, double,
                              double, double>;
    std::multiset<Record> expected, written;
    for (const GhostNode & ghost_node : expected_nodes)
    {
        Box * n = qtree->get_cell_box(ghost_node.large_node);
        Box * m = qtree->get_cell_box(ghost_node.small_node);
        Box * j1 = qtree->get_cell_box(ghost_node.j1);
        Box * j2 = qtree->get_cell_box(ghost_node.j2);
        expected.insert({n->layer, n->x, n->y, m->x, m->y, j1->x, j1->y, j2->x, j2->y,
                         ghost_node.alpha1, ghost_node.alpha2});
    }
    for (const DisvGhostNode & ghost_node : data.ghost_nodes)
    {
        const Cell2d & n = data.cell2d[ghost_node.large_cell];
        const Cell2d & m = data.cell2d[ghost_node.small_cell];
        const Cell2d & j1 = data.cell2d[ghost_node.contributor1];
        const Cell2d & j2 = data.cell2d[ghost_node.contributor2];
        written.insert({ghost_node.layer, n.x, n.y, m.x, m.y, j1.x, j1.y, j2.x, j2.y,
                        ghost_node.alpha1, ghost_node.alpha2});
    }
    EXPECT_EQ(written, expected);
    std::set<int> layers;
    for (const DisvGhostNode & ghost_node : data.ghost_nodes) layers.insert(ghost_node.layer);
    EXPECT_EQ(layers, (std::set<int>{0, 1}));
}

TEST(DisvExport, LayersRefinedDifferentlyAreAnError)
{
    std::unique_ptr<QuadTree3D> qtree =
        make_qtree(make_mfgrid(2, 2, 2, 100.0, 100.0), {{50.0, 150.0, 1, 0}});
    qtree->name = "qtg";
    EXPECT_THROW(disv_data(qtree.get(), false), GridgenError);
}

TEST(DisvExport, FilesAreWrittenOneBased)
{
    std::unique_ptr<QuadTree3D> qtree = make_northwest_refined(2);
    DisvData data = disv_data(qtree.get(), false);
    std::string prefix =
        (std::filesystem::path(::testing::TempDir()) / "gridgen_unit_disv").generic_string();
    write_disv_data(data, prefix);

    EXPECT_EQ(read_lines(prefix + ".dimensions.dat"),
              (std::vector<std::string>{"NLAY 2", "NCPL 7",
                                        "NVERT " + std::to_string(data.vertices.size())}));
    std::vector<std::string> cell2d = read_lines(prefix + ".cell2d.dat");
    ASSERT_EQ(cell2d.size(), 7u);
    EXPECT_EQ(cell2d[0].rfind("1 ", 0), 0u);
    EXPECT_EQ(read_lines(prefix + ".vertices.dat").size(), data.vertices.size());
    EXPECT_EQ(read_lines(prefix + ".top.dat").size(), 1u);
    EXPECT_EQ(read_lines(prefix + ".botm.dat").size(), 2u);
    EXPECT_EQ(read_lines(prefix + ".idomain.dat"),
              (std::vector<std::string>{"1 1 1 1 1 1 1", "1 1 1 1 1 1 1"}));

    std::vector<std::string> gnc = read_lines(prefix + ".gnc.dat");
    ASSERT_EQ(gnc.size(), data.ghost_nodes.size());
    const DisvGhostNode & first = data.ghost_nodes.front();
    std::istringstream line(gnc.front());
    int values[8];
    for (int & value : values) line >> value;
    EXPECT_EQ(values[0], first.layer + 1);
    EXPECT_EQ(values[1], first.large_cell + 1);
    EXPECT_EQ(values[3], first.small_cell + 1);
    EXPECT_EQ(values[5], first.contributor1 + 1);
    EXPECT_EQ(values[7], first.contributor2 + 1);
    for (const char * suffix : {".dimensions.dat", ".vertices.dat", ".cell2d.dat", ".top.dat",
                                ".botm.dat", ".idomain.dat", ".gnc.dat"})
        std::filesystem::remove(prefix + suffix);
}
