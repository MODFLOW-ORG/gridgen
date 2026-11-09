// Grid specification files, node data tables, and vertex and cell2d files
#include <fstream>
#include <sstream>
#include <string>
#include <vector>

#include <gtest/gtest.h>

#include "GridVertices.h"
#include "VertexExport.h"
#include "grid_helpers.h"

using namespace cusg;
using gridgen_test::make_mfgrid;
using gridgen_test::make_qtree;

namespace
{

// a 2 x 2 grid of 100 x 100 cells, nlay layers 10 thick, with the northwest
// cell of the top layer divided into four
std::unique_ptr<QuadTree3D> make_northwest_refined(int nlay)
{
    return make_qtree(make_mfgrid(nlay, 2, 2, 100.0, 100.0), {{50.0, 150.0, 1, 0}});
}

Box * base_cell(QuadTree3D * qtree, int layer, int row, int column)
{
    return qtree->node_group[qtree->get_modflow_grid()->get_box_index(layer, row, column)];
}

// the coordinates of a list of vertices
std::vector<Point2d> points(const GridVertices & gv, const std::vector<int> & ids)
{
    std::vector<Point2d> xy;
    for (int id : ids) xy.push_back(gv.vertices()[id]);
    return xy;
}

void expect_points(const std::vector<Point2d> & xy, const std::vector<Point2d> & expected)
{
    ASSERT_EQ(xy.size(), expected.size());
    for (size_t i = 0; i < xy.size(); i++)
    {
        EXPECT_EQ(xy[i][0], expected[i][0]) << "vertex " << i;
        EXPECT_EQ(xy[i][1], expected[i][1]) << "vertex " << i;
    }
}

std::vector<std::vector<double>> read_rows(const std::string & path)
{
    std::ifstream fin(path.c_str());
    std::vector<std::vector<double>> rows;
    std::string line;
    while (std::getline(fin, line))
    {
        std::istringstream ss(line);
        std::vector<double> row;
        double v;
        while (ss >> v) row.push_back(v);
        rows.push_back(row);
    }
    return rows;
}

} // namespace

TEST(VertexExport, OutlineAddsCornersOfFinerNeighbors)
{
    std::unique_ptr<QuadTree3D> qtree = make_northwest_refined(1);
    GridVertices gv(qtree.get());

    // northeast: the corner of two fine cells on its west edge
    expect_points(points(gv, gv.outline(base_cell(qtree.get(), 0, 0, 1))),
                  {Point2d(100, 200), Point2d(200, 200), Point2d(200, 100), Point2d(100, 100),
                   Point2d(100, 150)});
    // southwest: the corner on its north edge
    expect_points(
        points(gv, gv.outline(base_cell(qtree.get(), 0, 1, 0))),
        {Point2d(0, 100), Point2d(50, 100), Point2d(100, 100), Point2d(100, 0), Point2d(0, 0)});
    // southeast touches the fine cells only at its corner
    EXPECT_EQ(gv.outline(base_cell(qtree.get(), 0, 1, 1)).size(), 4u);
    // a fine cell has no finer neighbors
    Box * fine = base_cell(qtree.get(), 0, 0, 0)->p_children[1];
    EXPECT_EQ(gv.outline(fine).size(), 4u);
}

TEST(VertexExport, OutlineSkipsInactiveNeighbors)
{
    std::unique_ptr<QuadTree3D> qtree = make_northwest_refined(1);
    GridVertices gv(qtree.get());
    Box * northwest = base_cell(qtree.get(), 0, 0, 0);
    // the fine cells on the west edge of the northeast cell
    northwest->p_children[1]->active = false;
    EXPECT_EQ(gv.outline(base_cell(qtree.get(), 0, 0, 1)).size(), 5u);
    northwest->p_children[2]->active = false;
    EXPECT_EQ(gv.outline(base_cell(qtree.get(), 0, 0, 1)).size(), 4u);
}

TEST(VertexExport, OutlineOfModflowGridIsCorners)
{
    std::unique_ptr<ModflowGrid> mfgrid = make_mfgrid(2, 2, 3, 100.0, 50.0);
    GridVertices gv(mfgrid.get());
    for (auto box : mfgrid->node_group)
    {
        std::vector<int> ids = gv.outline(box);
        ASSERT_EQ(ids.size(), 4u);
        for (int i = 0; i < 4; i++) EXPECT_EQ(ids[i], gv.corners(box)[i]);
    }
}

TEST(VertexExport, GsfListsCornersAtTopAndBottom)
{
    // the fine cells are in the top layer only, so the layers differ
    std::unique_ptr<QuadTree3D> qtree = make_northwest_refined(2);
    GridVertices gv(qtree.get());
    std::string name = ::testing::TempDir() + "gridgen_vertex_export.gsf";
    write_gsf(qtree.get(), name, false);

    std::ifstream fin(name.c_str());
    std::string line;
    std::getline(fin, line);
    EXPECT_EQ(line, "UNSTRUCTURED GWF");
    int nnode, nlay, iz, ic, nvertex;
    fin >> nnode >> nlay >> iz >> ic >> nvertex;
    EXPECT_EQ(nnode, 11);
    EXPECT_EQ(nlay, 2);
    EXPECT_EQ(iz, 1);
    EXPECT_EQ(ic, 1);
    // 14 corners at the top and bottom of layer 1 and 9 at the bottom of layer 2
    EXPECT_EQ(nvertex, 37);
    std::vector<Point3d> vertices(nvertex);
    for (auto & v : vertices) fin >> v[0] >> v[1] >> v[2];

    for (int n = 1; n <= nnode; n++)
    {
        SCOPED_TRACE("node " + std::to_string(n));
        int node, layer, m;
        double x, y, z;
        fin >> node >> x >> y >> z >> layer >> m;
        ASSERT_EQ(node, n);
        ASSERT_EQ(m, 8);
        Box * box = qtree->get_cell_box(n - 1);
        EXPECT_EQ(layer, box->layer + 1);
        EXPECT_EQ(x, gv.center(box)[0]);
        EXPECT_EQ(y, gv.center(box)[1]);
        EXPECT_EQ(z, -5.0 - 10.0 * box->layer);
        for (int i = 0; i < m; i++)
        {
            int vertex;
            fin >> vertex;
            const Point3d & v = vertices[vertex - 1];
            const Point2d & corner = gv.vertices()[gv.corners(box)[i % 4]];
            EXPECT_EQ(v[0], corner[0]);
            EXPECT_EQ(v[1], corner[1]);
            EXPECT_EQ(v[2], -10.0 * (box->layer + i / 4));
        }
    }
}

TEST(VertexExport, GsfHangingVerticesAtTopAndBottom)
{
    std::unique_ptr<QuadTree3D> qtree = make_northwest_refined(2);
    std::string name = ::testing::TempDir() + "gridgen_vertex_export_hanging.gsf";
    write_gsf(qtree.get(), name, true);
    std::vector<std::vector<double>> rows = read_rows(name);
    // 14 vertices at the top and bottom of layer 1 and 9 at the bottom of layer 2
    ASSERT_EQ(rows[2][0], 37.0);

    // the northeast cell of the top layer, node 5, lists (100, 150) on its west edge
    const std::vector<double> & northeast = rows[3 + 37 + 4];
    ASSERT_EQ(northeast[0], 5.0);
    ASSERT_EQ(northeast[5], 10.0);
    for (int ring = 0; ring < 2; ring++)
    {
        const std::vector<double> & v = rows[3 + (int)northeast[6 + 5 * ring + 4] - 1];
        EXPECT_EQ(v[0], 100.0);
        EXPECT_EQ(v[1], 150.0);
        EXPECT_EQ(v[2], -10.0 * ring);
    }
    // the cells of the bottom layer have no finer neighbors
    for (int n = 8; n <= 11; n++) EXPECT_EQ(rows[3 + 37 + n - 1][5], 8.0) << "node " << n;
}

TEST(VertexExport, NodeDataTableNamesTheBaseCell)
{
    std::unique_ptr<QuadTree3D> qtree = make_northwest_refined(2);
    std::string name = ::testing::TempDir() + "gridgen_vertex_export.ndt";
    write_node_data_table(qtree.get(), name);
    std::vector<std::vector<double>> rows = read_rows(name);
    // node, layer, row, and column
    std::vector<std::vector<double>> expected = {
        {1, 1, 1, 1}, {2, 1, 1, 1}, {3, 1, 1, 1}, {4, 1, 1, 1},  {5, 1, 1, 2}, {6, 1, 2, 1},
        {7, 1, 2, 2}, {8, 2, 1, 1}, {9, 2, 1, 2}, {10, 2, 2, 1}, {11, 2, 2, 2}};
    EXPECT_EQ(rows, expected);
}

TEST(VertexExport, Cell2dListsHangingVertices)
{
    std::unique_ptr<QuadTree3D> qtree = make_northwest_refined(2);
    GridVertices gv(qtree.get());
    std::string prefix = ::testing::TempDir() + "gridgen_vertex_export";
    write_vertices_cell2d(qtree.get(), prefix);

    std::vector<std::vector<double>> vertices = read_rows(prefix + ".vertices.dat");
    ASSERT_EQ(vertices.size(), 14u);
    for (size_t i = 0; i < vertices.size(); i++) EXPECT_EQ(vertices[i][0], i + 1.0);

    std::vector<std::vector<double>> cells = read_rows(prefix + ".cell2d.dat");
    ASSERT_EQ(cells.size(), 11u);
    for (int n = 1; n <= 11; n++)
    {
        SCOPED_TRACE("node " + std::to_string(n));
        const std::vector<double> & cell = cells[n - 1];
        Box * box = qtree->get_cell_box(n - 1);
        EXPECT_EQ(cell[0], n);
        EXPECT_EQ(cell[1], gv.center(box)[0]);
        EXPECT_EQ(cell[2], gv.center(box)[1]);
        std::vector<Point2d> expected = points(gv, gv.outline(box));
        ASSERT_EQ(cell[3], (double)expected.size());
        ASSERT_EQ(cell.size(), 4 + expected.size());
        std::vector<Point2d> listed;
        for (size_t i = 0; i < expected.size(); i++)
        {
            const std::vector<double> & v = vertices[(int)cell[4 + i] - 1];
            listed.push_back(Point2d(v[1], v[2]));
        }
        expect_points(listed, expected);
    }
}
