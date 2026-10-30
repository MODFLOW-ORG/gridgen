// Grid intersection with polygons
#include <cmath>
#include <cstddef>
#include <set>
#include <vector>

#include <gtest/gtest.h>

#include "GridIntersection.h"
#include "grid_helpers.h"

using namespace cusg;
using gridgen_test::make_mfgrid;
using gridgen_test::make_qtree;

namespace
{

struct XY
{
    double x;
    double y;
};

// cells intersected by a polygon with the vertices in the order given
template <typename T>
std::vector<Box *> intersect_square(T * grid, const std::vector<XY> & vertices)
{
    c_ply ply(c_ply::POUT);
    ply.beginPoly();
    for (const XY & v : vertices) ply.addVertex(v.x, v.y);
    ply.endPoly();
    c_polygon polygon;
    polygon.push_back(ply);
    polygon.buildBoxAndCenter();

    GridIntersection gi(0, 0);
    gi.intersect(grid, polygon);
    std::vector<Box *> cells = gi.getIntersections();
    polygon.destroy();
    return cells;
}

// a square from lo to hi in x and y
std::vector<XY> square(double lo, double hi)
{
    return {{hi, lo}, {hi, hi}, {lo, hi}, {lo, lo}};
}

// every cell appears once in cells, and there are ncells of them
void expect_each_cell_once(const std::vector<Box *> & cells, std::size_t ncells)
{
    std::set<Box *> unique(cells.begin(), cells.end());
    EXPECT_EQ(unique.size(), cells.size());
    EXPECT_EQ(cells.size(), ncells);
}

bool contains(const Box * b, const XY & p)
{
    return std::abs(p.x - b->x) < 0.5 * b->dx && std::abs(p.y - b->y) < 0.5 * b->dy;
}

} // namespace

TEST(GridIntersection, PolygonCellsStartAtFirstVertex)
{
    // the cells are listed from the polygon's first vertex, whichever vertex
    // that is, rather than in memory-address order (#51)
    ModflowGrid * mfgrid = make_mfgrid(1, 10, 10, 10.0, 10.0);
    QuadTree3D * qtree = make_qtree(mfgrid, {});

    std::vector<XY> square = {{77.0, 23.0}, {77.0, 77.0}, {23.0, 77.0}, {23.0, 23.0}};
    std::set<Box *> first_cells;
    std::size_t ncells = 0;
    for (std::size_t start = 0; start < square.size(); start++)
    {
        SCOPED_TRACE("first vertex " + std::to_string(start));
        std::vector<XY> vertices;
        for (std::size_t i = 0; i < square.size(); i++)
            vertices.push_back(square[(start + i) % square.size()]);

        std::vector<Box *> cells = intersect_square(qtree, vertices);
        ASSERT_FALSE(cells.empty());
        EXPECT_TRUE(contains(cells.front(), vertices.front()));
        if (start == 0) ncells = cells.size();
        EXPECT_EQ(cells.size(), ncells);
        first_cells.insert(cells.front());
    }
    EXPECT_EQ(first_cells.size(), square.size());
}

TEST(GridIntersection, PolygonEnclosingGridSelectsEveryCell)
{
    // a polygon around the whole grid, with no edge crossing it, selects
    // every cell, and one beside the grid selects none (#53)
    ModflowGrid * mfgrid = make_mfgrid(1, 10, 10, 10.0, 10.0);
    expect_each_cell_once(intersect_square(mfgrid, square(-50.0, 150.0)), 100);
    EXPECT_TRUE(intersect_square(mfgrid, square(200.0, 300.0)).empty());

    QuadTree3D * qtree = make_qtree(make_mfgrid(1, 10, 10, 10.0, 10.0),
                                    {{55.0, 55.0, 2, 0}}, true);
    std::size_t nodes = static_cast<std::size_t>(qtree->nodes());
    ASSERT_GT(nodes, 100u);
    expect_each_cell_once(intersect_square(qtree, square(-50.0, 150.0)), nodes);
    EXPECT_TRUE(intersect_square(qtree, square(200.0, 300.0)).empty());
}
