// Grid intersection with polygons
#include <algorithm>
#include <cmath>
#include <cstddef>
#include <map>
#include <set>
#include <string>
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

// the area of a polygon with holes intersected with each cell of a grid, and
// whether each cell is listed once
template <typename T>
std::map<Box *, double> polygon_areas(T * grid, const std::vector<XY> & outer,
                                      const std::vector<std::vector<XY>> & holes, bool once = true)
{
    GIS_polygon polygon;
    static_cast<c_polygon &>(polygon) = make_polygon(outer, holes);
    BoxPlyClip * bpc = convertCPolygonToClip(polygon);

    GridIntersection gi(0, 0);
    gi.set_save_intersection_info(true);
    gi.intersect(grid, polygon, bpc);

    std::map<Box *, double> areas;
    for (const PolygonCellIntersectionInfo & info : gi.getPolygonIntersectionInfo())
    {
        if (once)
            EXPECT_EQ(areas.count(info.box), 0u)
                << "cell at " << info.box->x << ", " << info.box->y;
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
                areas = polygon_areas(make_qtree(mfgrid, {{45.0, 45.0, 2, 0}}, true),
                                      square(20.0, 80.0), {hole});
            else
                areas = polygon_areas(mfgrid, square(20.0, 80.0), {hole});

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

TEST(GridIntersection, PolygonVerticesOnGridLines)
{
    // the intersected areas sum to the polygon's area when its vertices, or
    // those of a hole, lie on grid lines, every time the intersection is
    // repeated; a hole can split a cell in two, so a cell can be listed twice
    // (#109)
    const std::vector<XY> diamond = {{30.0, 40.0}, {40.0, 30.0}, {50.0, 40.0}, {40.0, 50.0}};
    const std::vector<XY> triangle = {{30.0, 30.0}, {36.0, 30.0}, {44.0, 44.0}};
    // as a shapefile stores it, clockwise with the first vertex repeated at
    // the end
    const std::vector<XY> closed = {
        {30.0, 40.0}, {40.0, 50.0}, {50.0, 40.0}, {40.0, 30.0}, {30.0, 40.0}};
    // a vertex off a grid line by rounding, as in a rotated grid, on a grid
    // offset by 1000, where the rounding is one unit in the last place
    const std::vector<XY> rounded = {{1070.0 - 2.2737367544323206e-13, 1070.0},
                                     {1070.0, 1072.0},
                                     {1072.0, 1072.0},
                                     {1072.0, 1070.0}};
    // a hole vertex where a cell edge meets the outer boundary
    const std::vector<XY> notch = {{50.0, 30.0}, {45.0, 22.5}, {45.0, 35.0}};
    const std::vector<XY> star = {{65.0, 60.0}, {55.0, 85.0}, {37.5, 72.5}, {35.0, 50.0},
                                  {32.5, 27.5}, {55.0, 17.5}, {72.5, 32.5}};
    const std::vector<XY> kite = {{60.0, 52.5}, {47.5, 57.5}, {42.5, 45.0}, {50.0, 47.5}};
    // a hole that shares a vertex on a grid line with the outer ring
    const std::vector<XY> shared = {{70.0, 60.0}, {62.5, 87.5}, {27.5, 82.5},
                                    {12.5, 37.5}, {47.5, 27.5}, {67.5, 35.0}};
    const std::vector<XY> sharer = {{70.0, 60.0}, {57.5, 60.0}, {50.0, 60.0},
                                    {42.5, 55.0}, {55.0, 50.0}, {60.0, 52.5}};
    // an edge along a grid line, and one that ends inside a cell edge
    const std::vector<XY> edge = {
        {80.0, 80.0}, {40.0, 60.0}, {30.0, 50.0}, {50.0, 10.0}, {80.0, 20.0}};
    const std::vector<XY> shelf = {{95.0, 60.0}, {30.0, 60.0}, {65.0, 40.0}};
    struct Case
    {
        std::string name;
        std::vector<XY> outer;
        std::vector<std::vector<XY>> holes;
        double area;
        double offset;
        // cells in each row and column, and their size
        int n;
        double size;
    };
    const Case cases[] = {
        {"diamond", diamond, {}, 200.0, 0.0, 10, 10.0},
        {"diamond hole", square(5.0, 95.0), {diamond}, 7900.0, 0.0, 10, 10.0},
        {"triangle hole", square(5.0, 95.0), {triangle}, 8058.0, 0.0, 10, 10.0},
        {"closed diamond", closed, {}, 200.0, 0.0, 10, 10.0},
        {"rounded hole", square(1005.0, 1095.0), {rounded}, 8096.0, 1000.0, 10, 10.0},
        {"notch hole", square(5.0, 95.0), {notch}, 8068.75, 0.0, 10, 10.0},
        {"kite hole", star, {kite}, 1659.375, 0.0, 40, 2.5},
        {"shared vertex", shared, {sharer}, 2390.625, 0.0, 10, 10.0},
        {"edge on a grid line", edge, {}, 2300.0, 0.0, 10, 10.0},
        {"edge ending inside a cell edge", shelf, {}, 650.0, 0.0, 10, 10.0}};
    for (const Case & c : cases)
    {
        for (int refined = 0; refined < 2; refined++)
        {
            for (int repeat = 0; repeat < 5; repeat++)
            {
                SCOPED_TRACE(std::string(refined ? "quadtree" : "MODFLOW grid") + ", " + c.name +
                             ", repeat " + std::to_string(repeat));
                ModflowGrid * mfgrid =
                    make_mfgrid(1, c.n, c.n, c.size, c.size, 10.0, c.offset, c.offset);
                std::map<Box *, double> areas;
                if (refined)
                    areas = polygon_areas(
                        make_qtree(mfgrid, {{45.0 + c.offset, 45.0 + c.offset, 2, 0}}, true),
                        c.outer, c.holes, false);
                else
                    areas = polygon_areas(mfgrid, c.outer, c.holes, false);

                double total = 0.0;
                for (const auto & a : areas) total += a.second;
                EXPECT_NEAR(total, c.area, 1e-3);
            }
        }
    }
}

TEST(GridIntersection, LineAlongGridLine)
{
    // a line along a grid line that ends inside a cell edge selects the cells
    // on both sides of it (#109)
    ModflowGrid * mfgrid = make_mfgrid(1, 10, 10, 10.0, 10.0);
    GridIntersection gi(0, 0);
    gi.intersect(mfgrid, LineSeg2d(95.0, 60.0, 30.0, 60.0));
    expect_each_cell_once(gi.getIntersections(), 14);
}
