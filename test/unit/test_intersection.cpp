// Grid intersection with polygons
#include <algorithm>
#include <cmath>
#include <cstddef>
#include <map>
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

// a polygon from an outer ring and optional holes, vertices in the order given
c_polygon make_polygon(const std::vector<XY> & outer,
                       const std::vector<std::vector<XY>> & holes = {})
{
    c_polygon polygon;
    for (std::size_t r = 0; r <= holes.size(); r++)
    {
        c_ply ply(r == 0 ? c_ply::POUT : c_ply::PIN);
        ply.beginPoly();
        for (const XY & v : (r == 0 ? outer : holes[r - 1])) ply.addVertex(v.x, v.y);
        ply.endPoly();
        polygon.push_back(ply);
    }
    polygon.buildBoxAndCenter();
    return polygon;
}

// cells intersected by a polygon with the vertices in the order given
template <typename T>
std::vector<Box *> intersect_square(T * grid, const std::vector<XY> & vertices,
                                    const std::vector<std::vector<XY>> & holes = {})
{
    c_polygon polygon = make_polygon(vertices, holes);

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

    QuadTree3D * qtree = make_qtree(make_mfgrid(1, 10, 10, 10.0, 10.0), {{55.0, 55.0, 2, 0}}, true);
    std::size_t nodes = static_cast<std::size_t>(qtree->nodes());
    ASSERT_GT(nodes, 100u);
    expect_each_cell_once(intersect_square(qtree, square(-50.0, 150.0)), nodes);
    EXPECT_TRUE(intersect_square(qtree, square(200.0, 300.0)).empty());
}

TEST(GridIntersection, HoleEnclosingGridSelectsNoCell)
{
    // a hole that also encloses the whole grid leaves it outside the polygon
    ModflowGrid * mfgrid = make_mfgrid(1, 10, 10, 10.0, 10.0);
    std::vector<XY> hole = square(-20.0, 120.0);
    std::reverse(hole.begin(), hole.end());
    EXPECT_TRUE(intersect_square(mfgrid, square(-50.0, 150.0), {hole}).empty());
}

TEST(QuadTree, RefineByPolygonEnclosingGrid)
{
    // a refinement polygon around the whole grid refines every cell, as one
    // crossing the grid refines the cells inside it (#53)
    QuadTree3D * qtree = make_qtree(make_mfgrid(1, 10, 10, 10.0, 10.0), {});
    c_polygon polygon = make_polygon(square(-50.0, 150.0));
    qtree->refine(polygon, 1, 0);
    polygon.destroy();
    qtree->refresh();
    EXPECT_EQ(qtree->nodes(), 400);
}

namespace
{

// the area of a polygon with a hole intersected with each cell of a grid
template <typename T>
std::map<Box *, double> hole_areas(T * grid, const std::vector<XY> & outer,
                                   const std::vector<XY> & hole)
{
    GIS_polygon polygon;
    static_cast<c_polygon &>(polygon) = make_polygon(outer, {hole});
    BoxPlyClip * bpc = convertCPolygonToClip(polygon);

    GridIntersection gi(0, 0);
    gi.set_save_intersection_info(true);
    gi.intersect(grid, polygon, bpc);

    std::map<Box *, double> areas;
    for (const PolygonCellIntersectionInfo & info : gi.getPolygonIntersectionInfo())
    {
        EXPECT_EQ(areas.count(info.box), 0u) << "cell at " << info.box->x << ", " << info.box->y;
        areas[info.box] += info.m_area;
    }
    bpc->destroy();
    delete bpc;
    polygon.destroy();
    return areas;
}

} // namespace

TEST(GridIntersection, PolygonHoleIsNotPartOfTheArea)
{
    // the intersected areas of a 60 by 60 square sum to its area less its
    // hole's, each cell is listed once, and cells in the hole are not, for a
    // hole on the cell edges, one that cuts cells, and one inside a cell (#99)
    struct Case
    {
        double lo;
        double hi;
        double area;
    };
    const Case holes[] = {{40.0, 60.0, 3200.0}, {43.0, 57.0, 3404.0}, {41.0, 48.0, 3551.0}};
    for (const Case & c : holes)
    {
        std::vector<XY> hole = square(c.lo, c.hi);
        std::reverse(hole.begin(), hole.end());
        for (int refined = 0; refined < 2; refined++)
        {
            SCOPED_TRACE(std::string(refined ? "quadtree" : "MODFLOW grid") + ", hole from " +
                         std::to_string(c.lo) + " to " + std::to_string(c.hi));
            ModflowGrid * mfgrid = make_mfgrid(1, 10, 10, 10.0, 10.0);
            std::map<Box *, double> areas;
            if (refined)
                areas = hole_areas(make_qtree(mfgrid, {{45.0, 45.0, 2, 0}}, true),
                                   square(20.0, 80.0), hole);
            else
                areas = hole_areas(mfgrid, square(20.0, 80.0), hole);

            double total = 0.0;
            for (const auto & a : areas)
            {
                const Box * b = a.first;
                total += a.second;
                EXPECT_FALSE(b->x - b->dx / 2 >= c.lo && b->x + b->dx / 2 <= c.hi &&
                             b->y - b->dy / 2 >= c.lo && b->y + b->dy / 2 <= c.hi)
                    << "cell in the hole at " << b->x << ", " << b->y;
            }
            EXPECT_NEAR(total, c.area, 1e-3);
        }
    }
}
