// Shared-vertex VTK export
#include <cmath>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>

#include <gtest/gtest.h>

#include "GridToVTK.h"
#include "grid_helpers.h"

using namespace cusg;
using gridgen_test::make_mfgrid;

namespace
{

// cell top elevation, linear in x and y
double top_at(double x, double y) { return 0.1 * x + 0.2 * y; }

// points of the vtu file as x, y, z triples
std::vector<double> read_vtu_points(const std::string & path)
{
    std::ifstream fin(path.c_str());
    std::string line;
    while (std::getline(fin, line))
        if (line.find("<Points>") != std::string::npos) break;
    std::getline(fin, line); // DataArray tag
    std::string values;
    while (std::getline(fin, line) && line.find("</DataArray>") == std::string::npos)
        values += line + " ";
    std::istringstream ss(values);
    std::vector<double> xyz;
    double v;
    while (ss >> v) xyz.push_back(v);
    return xyz;
}

// z of the center vertices of a 2 x 2 grid with cell (rf, cf) refined and
// cell (ri, ci) inactive, so the vertex is shared by one fine and two coarse
// cells
std::vector<double> center_vertex_z(int rf, int cf, int ri, int ci)
{
    ModflowGrid * mfgrid = make_mfgrid(1, 2, 2, 100.0, 100.0);
    QuadTree3D * qtree = new QuadTree3D(mfgrid);
    Point2d pt(50.0 + 100.0 * cf, 150.0 - 100.0 * rf);
    qtree->refine(pt, 1, 0);
    qtree->nodegroup[mfgrid->get_nodeid(0, ri, ci)]->active = false;
    qtree->set_one_based_node_numbering(false);
    qtree->number_nodes();
    qtree->refresh();

    double * top = qtree->get_top();
    double * bot = qtree->get_bot();
    for (int n = 0; n < qtree->nodes(); n++)
    {
        Box * b = qtree->get_nodeobj(n);
        top[n] = top_at(b->x, b->y);
        bot[n] = top[n] - 10.0;
    }

    std::string name = ::testing::TempDir() + "gridgen_sharedvertex";
    writeQuadtreeGrid2VTK_shareV(qtree, name);
    std::vector<double> xyz = read_vtu_points(name + ".vtu");

    std::vector<double> zs;
    for (size_t i = 0; i + 2 < xyz.size(); i += 3)
    {
        if (std::fabs(xyz[i] - 100.0) < 1e-3 && std::fabs(xyz[i + 1] - 100.0) < 1e-3)
            zs.push_back(xyz[i + 2]);
    }

    delete qtree;
    delete mfgrid;
    return zs;
}

} // namespace

TEST(VTK, SharedVertexPlaneFitAtThreeCells)
{
    // the plane through the three cell centers reproduces the linear surface
    // for every arrangement of the cells (#45)
    double ztop = top_at(100.0, 100.0);
    for (int f = 0; f < 4; f++)
    {
        for (int i = 0; i < 4; i++)
        {
            if (i == f) continue;
            SCOPED_TRACE("refined " + std::to_string(f) + ", inactive " +
                         std::to_string(i));
            std::vector<double> zs = center_vertex_z(f / 2, f % 2, i / 2, i % 2);
            ASSERT_EQ(zs.size(), 2u);
            for (double z : zs)
                EXPECT_NEAR(z, z > ztop - 5.0 ? ztop : ztop - 10.0, 1e-4);
        }
    }
}
