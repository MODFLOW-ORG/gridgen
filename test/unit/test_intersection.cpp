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

// cells intersected by a square polygon with the vertices in the order given
std::vector<Box *> intersect_square(QuadTree3D * qtree, const std::vector<XY> & vertices)
{
    c_ply ply(c_ply::POUT);
    ply.beginPoly();
    for (const XY & v : vertices) ply.addVertex(v.x, v.y);
    ply.endPoly();
    c_polygon polygon;
    polygon.push_back(ply);
    polygon.buildBoxAndCenter();

    GridIntersection gi(0, 0);
    EXPECT_TRUE(gi.intersect(qtree, polygon));
    std::vector<Box *> cells = gi.getIntersections();
    polygon.destroy();
    return cells;
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
