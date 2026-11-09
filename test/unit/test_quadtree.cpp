// Quadtree refinement, smoothing, and connectivity
#include <algorithm>
#include <cstdlib>
#include <set>
#include <string>
#include <tuple>
#include <utility>
#include <vector>

#include <gtest/gtest.h>

#include "builder_helpers.h"
#include "csr_checks.h"
#include "grid_helpers.h"
#include "polyline.h"

using namespace cusg;
using gridgen_test::expect_csr_valid;
using gridgen_test::expect_faces_cover_sides;
using gridgen_test::expect_symmetric;
using gridgen_test::make_mfgrid;
using gridgen_test::make_qtree;

namespace
{

// largest refinement level difference across horizontal connections
int max_level_jump(QuadTree3D * qtree, const CSRData & csr, bool vertical)
{
    int jump = 0;
    for (int n = 0; n < csr.nodes; n++)
    {
        for (int ipos = csr.ia[n] + 1; ipos < csr.ia[n + 1]; ipos++)
        {
            bool isvert = std::abs(csr.fldr[ipos]) == 3;
            if (isvert != vertical) continue;
            int dn = qtree->get_cell_box(n)->depth;
            int dm = qtree->get_cell_box(csr.ja[ipos])->depth;
            jump = std::max(jump, std::abs(dn - dm));
        }
    }
    return jump;
}

// refine 3 x 3 cells of 100 once along a line; return the cells refined,
// numbered by row from the top and then by column
std::set<int> line_refined_cells(const std::vector<std::pair<double, double>> & points)
{
    std::unique_ptr<ModflowGrid> mfgrid = make_mfgrid(1, 3, 3, 100.0, 100.0);
    std::unique_ptr<QuadTree3D> qtree = std::make_unique<QuadTree3D>(mfgrid.get());
    Polyline arc;
    arc.begin_poly();
    for (const auto & pt : points) arc.add_vertex(pt.first, pt.second);
    arc.end_poly();
    qtree->refine(arc, 1, 0);
    arc.destroy();
    qtree->number_nodes();
    qtree->refresh();

    std::set<int> cells;
    for (int n = 0; n < qtree->nodes(); n++)
    {
        Box * box = qtree->get_cell_box(n);
        if (box->depth == 0) continue;
        int row = static_cast<int>((300.0 - box->y) / 100.0);
        int col = static_cast<int>(box->x / 100.0);
        cells.insert(3 * row + col);
    }
    return cells;
}

// the level and center of each cell of the 3 x 3 grid of builder_dfn refined
// by shapes with REFINEMENT_LEVEL_BY_ATTRIBUTE = LEVEL
std::vector<std::tuple<int, double, double>>
attribute_refined_cells(const std::string & name, int type, const std::string & feature_type,
                        const std::vector<gridgen_test::Shape> & shapes)
{
    std::string path = gridgen_test::write_shapes(name, type, shapes);
    std::vector<std::tuple<int, double, double>> cells;
    gridgen_test::build_qtree(
        gridgen_test::builder_dfn(path, feature_type, "REFINEMENT_LEVEL_BY_ATTRIBUTE = LEVEL"),
        nullptr, &cells);
    return cells;
}

} // namespace

TEST(QuadTree, RefineCenterCell)
{
    // 3 x 3 cells of 100; refine the center cell once
    std::unique_ptr<ModflowGrid> mfgrid = make_mfgrid(1, 3, 3, 100.0, 100.0);
    std::unique_ptr<QuadTree3D> qtree = make_qtree(mfgrid.get(), {{150.0, 150.0, 1, 0}});
    CSRData csr = qtree->get_usg_csr_data();

    EXPECT_EQ(csr.nodes, 12);
    expect_faces_cover_sides(qtree.get(), csr, 10.0);

    double total = 0.0;
    int nsmall = 0;
    for (int n = 0; n < csr.nodes; n++)
    {
        total += csr.area[n];
        if (csr.area[n] == 2500.0) nsmall++;
    }
    EXPECT_DOUBLE_EQ(total, 90000.0);
    EXPECT_EQ(nsmall, 4);

    // fine-coarse faces have the fine cell width and half lengths of each cell
    for (int n = 0; n < csr.nodes; n++)
    {
        for (int ipos = csr.ia[n] + 1; ipos < csr.ia[n + 1]; ipos++)
        {
            if (csr.area[n] != 2500.0 || csr.area[csr.ja[ipos]] != 10000.0)
            {
                continue;
            }
            EXPECT_DOUBLE_EQ(csr.cl1[ipos], 25.0);
            EXPECT_DOUBLE_EQ(csr.cl2[ipos], 50.0);
            EXPECT_DOUBLE_EQ(csr.fahl[ipos], 500.0);
        }
    }
    expect_symmetric(csr);
    expect_csr_valid(csr);
}

TEST(QuadTree, SmoothingLimitsHorizontalLevelJump)
{
    std::unique_ptr<QuadTree3D> rough =
        make_qtree(make_mfgrid(1, 5, 5, 100.0, 100.0), {{250.0, 250.0, 3, 0}}, false);
    CSRData csr = rough->get_usg_csr_data();
    EXPECT_GT(max_level_jump(rough.get(), csr, false), 1);

    std::unique_ptr<QuadTree3D> smooth =
        make_qtree(make_mfgrid(1, 5, 5, 100.0, 100.0), {{250.0, 250.0, 3, 0}}, true);
    csr = smooth->get_usg_csr_data();
    EXPECT_EQ(max_level_jump(smooth.get(), csr, false), 1);
    expect_faces_cover_sides(smooth.get(), csr, 10.0);
    expect_symmetric(csr);
    expect_csr_valid(csr);
}

TEST(QuadTree, SmoothingLimitsVerticalLevelJump)
{
    // refine only the top layer; smoothing propagates refinement downward
    std::unique_ptr<ModflowGrid> mfgrid = make_mfgrid(3, 4, 4, 100.0, 100.0);
    std::unique_ptr<QuadTree3D> qtree = make_qtree(mfgrid.get(), {{150.0, 150.0, 3, 0}}, true);
    CSRData csr = qtree->get_usg_csr_data();

    EXPECT_EQ(max_level_jump(qtree.get(), csr, true), 1);
    expect_faces_cover_sides(qtree.get(), csr, 10.0);
    expect_symmetric(csr);
    expect_csr_valid(csr);

    // a vertical face is the area of the smaller cell
    for (int n = 0; n < csr.nodes; n++)
    {
        for (int ipos = csr.ia[n] + 1; ipos < csr.ia[n + 1]; ipos++)
        {
            if (std::abs(csr.fldr[ipos]) != 3) continue;
            double amin = std::min(csr.area[n], csr.area[csr.ja[ipos]]);
            EXPECT_DOUBLE_EQ(csr.fahl[ipos], amin);
        }
    }
}

TEST(QuadTree, VerticalPassThroughInactiveCell)
{
    // 3 layers of 3 x 3 cells; the center cell of layer 2 is inactive, so
    // with pass-through the center cells of layers 1 and 3 connect
    std::unique_ptr<QuadTree3D> qtree = make_qtree(make_mfgrid(3, 3, 3, 100.0, 100.0), {});
    Box * center = nullptr;
    for (int n = 0; n < qtree->nodes(); n++)
    {
        Box * box = qtree->get_cell_box(n);
        if (box->layer == 1 && box->x == 150.0 && box->y == 150.0) center = box;
    }
    ASSERT_NE(center, nullptr);
    center->active = false;
    qtree->number_nodes();
    qtree->refresh();

    qtree->vertical_pass_through = true;
    CSRData csr = qtree->get_usg_csr_data();
    EXPECT_EQ(csr.nodes, 26);

    int nlong = 0;
    for (int n = 0; n < csr.nodes; n++)
    {
        for (int ipos = csr.ia[n] + 1; ipos < csr.ia[n + 1]; ipos++)
        {
            Box * bn = qtree->get_cell_box(n);
            Box * bm = qtree->get_cell_box(csr.ja[ipos]);
            if (std::abs(bn->layer - bm->layer) == 2) nlong++;
        }
    }
    // layer 1 to layer 3 and back
    EXPECT_EQ(nlong, 2);
}

TEST(QuadTree, VerticalPassThroughConnectsOnlyOverlappingCells)
{
    // 2 x 2 cells of 100 in 3 layers: the northwest cell of layer 1 refined
    // once, the northwest cell of layer 2 inactive, and layer 3 refined twice;
    // each fine cell of layer 1 connects only to the 4 cells of layer 3 below
    // it, over the area of a layer 3 cell
    std::unique_ptr<ModflowGrid> mfgrid = make_mfgrid(3, 2, 2, 100.0, 100.0);
    std::unique_ptr<QuadTree3D> qtree = make_qtree(mfgrid.get(), {{50.0, 150.0, 1, 0},
                                                                  {50.0, 150.0, 2, 2},
                                                                  {150.0, 150.0, 2, 2},
                                                                  {50.0, 50.0, 2, 2},
                                                                  {150.0, 50.0, 2, 2}});
    qtree->node_group[mfgrid->get_box_index(1, 0, 0)]->active = false;
    qtree->number_nodes();
    qtree->refresh();
    qtree->vertical_pass_through = true;

    CSRData csr = qtree->get_usg_csr_data();
    expect_symmetric(csr);
    expect_csr_valid(csr);

    int across = 0;
    for (int n = 0; n < csr.nodes; n++)
    {
        for (int ipos = csr.ia[n] + 1; ipos < csr.ia[n + 1]; ipos++)
        {
            if (std::abs(csr.fldr[ipos]) != 3) continue;
            Box * bn = qtree->get_cell_box(n);
            Box * bm = qtree->get_cell_box(csr.ja[ipos]);
            Box * smaller = (csr.area[n] < csr.area[csr.ja[ipos]]) ? bn : bm;
            Box * larger = (smaller == bn) ? bm : bn;
            EXPECT_TRUE(larger->in(smaller->x, smaller->y)) << "n " << n << " m " << csr.ja[ipos];
            EXPECT_DOUBLE_EQ(csr.fahl[ipos], std::min(csr.area[n], csr.area[csr.ja[ipos]]));
            if (std::abs(bn->layer - bm->layer) == 2) across++;
        }
    }
    // 4 fine cells of layer 1 to 4 cells of layer 3 each, counted from both
    EXPECT_EQ(across, 32);
}

TEST(QuadTree, VerticalPassThroughAcrossSeveralInactiveCells)
{
    // 2 x 2 cells of 100 in 4 layers. Northwest: layer 1 refined once, layer 2
    // inactive, layer 3 refined once with all but its southeast cell inactive,
    // and layer 4 refined twice. Northeast: layer 2 refined once with its two
    // west cells inactive. Each cell of layer 1 is covered below in full
    std::unique_ptr<ModflowGrid> mfgrid = make_mfgrid(4, 2, 2, 100.0, 100.0);
    std::unique_ptr<QuadTree3D> qtree = make_qtree(
        mfgrid.get(),
        {{50.0, 150.0, 1, 0}, {50.0, 150.0, 1, 2}, {50.0, 150.0, 2, 3}, {150.0, 150.0, 1, 1}});
    qtree->node_group[mfgrid->get_box_index(1, 0, 0)]->active = false;
    for (int n = 0; n < qtree->nodes(); n++)
    {
        Box * box = qtree->get_cell_box(n);
        bool northwest = box->x < 100.0 && box->y > 100.0;
        if (box->layer == 2 && northwest && !(box->x == 75.0 && box->y == 125.0))
            box->active = false;
        if (box->layer == 1 && box->x == 125.0) box->active = false;
    }
    qtree->number_nodes();
    qtree->refresh();
    qtree->vertical_pass_through = true;

    CSRData csr = qtree->get_usg_csr_data();
    expect_symmetric(csr);
    expect_csr_valid(csr);

    int layer_1_to_4 = 0;
    for (int n = 0; n < csr.nodes; n++)
    {
        Box * bn = qtree->get_cell_box(n);
        double covered_below = 0.0;
        for (int ipos = csr.ia[n] + 1; ipos < csr.ia[n + 1]; ipos++)
        {
            if (std::abs(csr.fldr[ipos]) != 3) continue;
            Box * bm = qtree->get_cell_box(csr.ja[ipos]);
            Box * smaller = (csr.area[n] < csr.area[csr.ja[ipos]]) ? bn : bm;
            Box * larger = (smaller == bn) ? bm : bn;
            EXPECT_TRUE(larger->in(smaller->x, smaller->y)) << "n " << n << " m " << csr.ja[ipos];
            EXPECT_LE(csr.fahl[ipos], std::min(csr.area[n], csr.area[csr.ja[ipos]]));
            if (csr.fldr[ipos] == -3) covered_below += csr.fahl[ipos];
            if (bn->layer == 0 && bm->layer == 3) layer_1_to_4++;
            // the two inactive quarter cells of layer 2 below the northeast cell
            if (bn->layer == 0 && bm->layer == 2 && bn->x == 150.0)
                EXPECT_DOUBLE_EQ(csr.fahl[ipos], 5000.0);
        }
        if (bn->layer == 0) EXPECT_DOUBLE_EQ(covered_below, csr.area[n]) << "n " << n;
    }
    // 3 fine cells of layer 1, each above 4 cells of layer 4
    EXPECT_EQ(layer_1_to_4, 12);
}

TEST(QuadTree, BoxWithNoCellNumberThrows)
{
    // 3 x 3 cells of 100; the center cell is split and a corner cell inactive
    std::unique_ptr<ModflowGrid> mfgrid = make_mfgrid(1, 3, 3, 100.0, 100.0);
    std::unique_ptr<QuadTree3D> qtree = make_qtree(mfgrid.get(), {{150.0, 150.0, 1, 0}});
    Box * corner = qtree->get_cell_box(0);
    corner->active = false;
    qtree->number_nodes();
    qtree->refresh();
    qtree->get_cell_array();
    ASSERT_EQ(qtree->nodes(), 11);

    EXPECT_THROW(qtree->get_cell(corner->id), GridgenError);

    int split = -1;
    for (int box_index = 0; box_index < (int)qtree->node_group.size(); box_index++)
    {
        if (!qtree->node_group[box_index]->is_leaf) split = box_index;
    }
    ASSERT_NE(split, -1);
    EXPECT_THROW(qtree->get_cell(split), GridgenError);

    // every active cell maps back to itself
    for (int cell = 0; cell < qtree->nodes(); cell++)
    {
        EXPECT_EQ(qtree->get_cell(qtree->get_cell_box(cell)->id), cell);
    }
}

TEST(QuadTree, LineRefinesCellsItCrosses)
{
    EXPECT_EQ(line_refined_cells({{50.0, 150.0}, {250.0, 150.0}}), (std::set<int>{3, 4, 5}));
}

TEST(QuadTree, LineRefinesCellItIsInside)
{
    EXPECT_EQ(line_refined_cells({{120.0, 120.0}, {180.0, 180.0}}), (std::set<int>{4}));
}

TEST(QuadTree, LineEndingOnEdgeDoesNotRefineNextCell)
{
    EXPECT_EQ(line_refined_cells({{150.0, 250.0}, {150.0, 200.0}}), (std::set<int>{1}));
}

TEST(QuadTree, LineEndingAtCornerDoesNotRefineNextCells)
{
    EXPECT_EQ(line_refined_cells({{150.0, 250.0}, {200.0, 200.0}}), (std::set<int>{1}));
}

TEST(QuadTree, LineThroughCornerRefinesCellsAtCorner)
{
    EXPECT_EQ(line_refined_cells({{150.0, 250.0}, {250.0, 150.0}}), (std::set<int>{1, 2, 4, 5}));
}

TEST(QuadTree, LineAlongEdgeRefinesCellsOnBothSides)
{
    EXPECT_EQ(line_refined_cells({{100.0, 250.0}, {100.0, 150.0}}), (std::set<int>{0, 1, 3, 4}));
}

TEST(QuadTree, LineEqualToEdgeRefinesCellsOnBothSides)
{
    EXPECT_EQ(line_refined_cells({{100.0, 200.0}, {100.0, 100.0}}), (std::set<int>{3, 4}));
}

TEST(QuadTree, PolylineVertexOnEdgeRefinesCellsEachSegmentEnters)
{
    EXPECT_EQ(line_refined_cells({{150.0, 250.0}, {150.0, 200.0}, {250.0, 150.0}}),
              (std::set<int>{1, 4, 5}));
}

TEST(QuadTree, LineWithEndsOutsideGridRefinesCellsItCrosses)
{
    EXPECT_EQ(line_refined_cells({{-50.0, 150.0}, {350.0, 150.0}}), (std::set<int>{3, 4, 5}));
}

TEST(QuadTree, PolylineSegmentsWithEndsOutsideGridRefineCellsTheyCross)
{
    // the second segment is outside the grid
    EXPECT_EQ(line_refined_cells({{-50.0, 250.0}, {350.0, 250.0}, {350.0, 50.0}, {-50.0, 50.0}}),
              (std::set<int>{0, 1, 2, 6, 7, 8}));
}

TEST(QuadTreeBuilder, LineWithEndsOutsideRotatedGridRefinesCellsItCrosses)
{
    // rotated 90 degrees, the grid covers x from -300 to 0, and the line
    // crosses its middle row
    std::string path = gridgen_test::write_test_shapefile("gridgen_unit_rotated_line", SHPT_ARC,
                                                          {{{-150.0, -50.0}, {-150.0, 350.0}}});
    int nodes = 0;
    gridgen_test::build_qtree(gridgen_test::builder_dfn(path, "line", "REFINEMENT_LEVEL = 1", "90"),
                              &nodes);
    EXPECT_EQ(nodes, 18);
}

// each feature of a shape is refined to the level of that shape (issue #133)
TEST(QuadTreeBuilder, MultipartLineIsRefinedToItsShapeLevel)
{
    gridgen_test::Vertices a = {{10.0, 250.0}, {90.0, 250.0}}, b = {{210.0, 250.0}, {290.0, 250.0}},
                           c = {{110.0, 50.0}, {190.0, 50.0}};
    auto cells = attribute_refined_cells("gridgen_unit_level_lines", SHPT_ARC, "line",
                                         {{{a, b}, 2}, {{c}, 1}});
    EXPECT_EQ(cells.size(), 9u + 2 * 15 + 3);
    EXPECT_EQ(cells, attribute_refined_cells("gridgen_unit_level_lines1", SHPT_ARC, "line",
                                             {{{a}, 2}, {{b}, 2}, {{c}, 1}}));
}

TEST(QuadTreeBuilder, MultipointIsRefinedToItsShapeLevel)
{
    gridgen_test::Vertices ab = {{50.0, 250.0}, {250.0, 250.0}}, c = {{150.0, 50.0}};
    auto cells = attribute_refined_cells("gridgen_unit_level_points", SHPT_MULTIPOINT, "point",
                                         {{{ab}, 2}, {{c}, 1}});
    EXPECT_EQ(cells.size(), 9u + 2 * 15 + 3);
    EXPECT_EQ(cells, attribute_refined_cells("gridgen_unit_level_points1", SHPT_MULTIPOINT, "point",
                                             {{{{ab[0]}}, 2}, {{{ab[1]}}, 2}, {{c}, 1}}));
}

TEST(QuadTreeBuilder, MultipatchIsRefinedToItsShapeLevel)
{
    gridgen_test::Vertices a = gridgen_test::square(20.0, 220.0, 40.0),
                           b = gridgen_test::square(220.0, 220.0, 40.0),
                           c = gridgen_test::square(120.0, 20.0, 40.0);
    auto cells = attribute_refined_cells("gridgen_unit_level_polygons", SHPT_MULTIPATCH, "polygon",
                                         {{{a, b}, 2}, {{c}, 1}});
    EXPECT_EQ(cells.size(), 9u + 2 * 15 + 3);
    EXPECT_EQ(cells, attribute_refined_cells("gridgen_unit_level_polygons1", SHPT_MULTIPATCH,
                                             "polygon", {{{a}, 2}, {{b}, 2}, {{c}, 1}}));
}

TEST(QuadTreeBuilder, NullShapeDoesNotShiftLevels)
{
    gridgen_test::Vertices a = {{10.0, 250.0}, {90.0, 250.0}};
    auto cells =
        attribute_refined_cells("gridgen_unit_level_null", SHPT_ARC, "line", {{{}, 2}, {{a}, 1}});
    EXPECT_EQ(cells.size(), 9u + 3);
    EXPECT_EQ(cells,
              attribute_refined_cells("gridgen_unit_level_null1", SHPT_ARC, "line", {{{a}, 1}}));
}

// each outer ring of a polygon shape is a polygon, and each hole belongs to
// the ring it is inside, whatever the order and orientation of the rings
TEST(QuadTreeBuilder, PolygonWithTwoOuterRingsRefinesBoth)
{
    gridgen_test::Vertices a = gridgen_test::square(20.0, 220.0, 40.0),
                           b = gridgen_test::square(220.0, 220.0, 40.0);
    auto cells =
        attribute_refined_cells("gridgen_unit_two_rings", SHPT_POLYGON, "polygon", {{{a, b}, 2}});
    EXPECT_EQ(cells.size(), 9u + 2 * 15);
    EXPECT_EQ(cells, attribute_refined_cells("gridgen_unit_two_rings1", SHPT_POLYGON, "polygon",
                                             {{{a}, 2}, {{b}, 2}}));
}

TEST(QuadTreeBuilder, PolygonHoleIsNotRefined)
{
    // the hole holds the center cell, so the cells around it are refined
    gridgen_test::Vertices outer = gridgen_test::square(10.0, 10.0, 280.0),
                           hole = gridgen_test::square(90.0, 90.0, 120.0);
    auto cells =
        attribute_refined_cells("gridgen_unit_hole", SHPT_POLYGON, "polygon", {{{outer, hole}, 1}});
    EXPECT_EQ(cells.size(), 1u + 8 * 4);
    EXPECT_EQ(cells, attribute_refined_cells("gridgen_unit_hole1", SHPT_POLYGON, "polygon",
                                             {{{hole, outer}, 1}}));
}

TEST(QuadTreeBuilder, PolygonWithTwoOuterRingsAndHoles)
{
    gridgen_test::Vertices a = gridgen_test::square(10.0, 210.0, 80.0),
                           ahole = gridgen_test::square(30.0, 230.0, 40.0),
                           b = gridgen_test::square(210.0, 10.0, 80.0),
                           bhole = gridgen_test::square(230.0, 30.0, 40.0);
    auto separate = attribute_refined_cells("gridgen_unit_holes1", SHPT_POLYGON, "polygon",
                                            {{{a, ahole}, 1}, {{b, bhole}, 1}});
    EXPECT_EQ(separate.size(), 9u + 2 * 3);
    EXPECT_EQ(attribute_refined_cells("gridgen_unit_holes", SHPT_POLYGON, "polygon",
                                      {{{a, ahole, b, bhole}, 1}}),
              separate);
    EXPECT_EQ(attribute_refined_cells("gridgen_unit_holes2", SHPT_POLYGON, "polygon",
                                      {{{a, b, bhole, ahole}, 1}}),
              separate);
}

TEST(QuadTreeBuilder, PolygonHoleTouchingOuterRingIsAHole)
{
    // the hole meets the outer ring along x = 290 and holds the middle right
    // cell, which is not refined whichever vertex the hole starts at
    gridgen_test::Vertices outer = gridgen_test::square(10.0, 10.0, 280.0);
    gridgen_test::Vertices hole = {
        {190.0, 90.0}, {190.0, 210.0}, {290.0, 210.0}, {290.0, 90.0}, {190.0, 90.0}};
    gridgen_test::Vertices hole_from_edge = {
        {290.0, 210.0}, {290.0, 90.0}, {190.0, 90.0}, {190.0, 210.0}, {290.0, 210.0}};
    auto cells = attribute_refined_cells("gridgen_unit_touching_hole", SHPT_POLYGON, "polygon",
                                         {{{outer, hole}, 1}});
    EXPECT_EQ(cells.size(), 1u + 8 * 4);
    EXPECT_NE(std::find(cells.begin(), cells.end(), std::make_tuple(0, 250.0, 150.0)), cells.end());
    EXPECT_EQ(cells, attribute_refined_cells("gridgen_unit_touching_hole1", SHPT_POLYGON, "polygon",
                                             {{{outer, hole_from_edge}, 1}}));
}

TEST(QuadTreeBuilder, PolygonIslandInHoleIsAPolygon)
{
    gridgen_test::Vertices outer = gridgen_test::square(10.0, 10.0, 280.0),
                           hole = gridgen_test::square(60.0, 60.0, 180.0),
                           island = gridgen_test::square(140.0, 140.0, 20.0);
    EXPECT_EQ(attribute_refined_cells("gridgen_unit_island", SHPT_POLYGON, "polygon",
                                      {{{island, outer, hole}, 1}}),
              attribute_refined_cells("gridgen_unit_island1", SHPT_POLYGON, "polygon",
                                      {{{outer, hole}, 1}, {{island}, 1}}));
}
