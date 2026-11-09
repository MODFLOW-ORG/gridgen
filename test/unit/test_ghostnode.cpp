// Horizontal and vertical ghost node correction
#include <cmath>
#include <functional>
#include <set>
#include <utility>
#include <vector>

#include <gtest/gtest.h>

#include "GhostNode.h"
#include "grid_helpers.h"

using namespace cusg;
using gridgen_test::make_mfgrid;
using gridgen_test::make_qtree;

namespace
{

// the connections of a grid as (n, m) pairs
std::set<std::pair<int, int>> connections(QuadTree3D * qtree)
{
    CSRData csr = qtree->get_usg_csr_data();
    std::set<std::pair<int, int>> connected;
    for (int n = 0; n < csr.nodes; n++)
    {
        for (int ipos = csr.ia[n] + 1; ipos < csr.ia[n + 1]; ipos++)
        {
            connected.insert(std::make_pair(n, csr.ja[ipos]));
        }
    }
    return connected;
}

bool is_vertical(QuadTree3D * qtree, const GhostNode & ghost_node)
{
    return qtree->get_cell_box(ghost_node.large_node)->layer !=
           qtree->get_cell_box(ghost_node.small_node)->layer;
}

// the head interpolated at the ghost node from a head field
double ghost_head(QuadTree3D * qtree, const GhostNode & ghost_node,
                  const std::function<double(Box *)> & head)
{
    Box * n = qtree->get_cell_box(ghost_node.large_node);
    Box * j1 = qtree->get_cell_box(ghost_node.j1);
    Box * j2 = qtree->get_cell_box(ghost_node.j2);
    return head(n) + ghost_node.alpha1 * (head(j1) - head(n)) +
           ghost_node.alpha2 * (head(j2) - head(n));
}

// the head of a linear field at the ghost node, which is below or above the
// center of m for a vertical connection and level with the center of n
// across a horizontal one
double ghost_location_head(QuadTree3D * qtree, const GhostNode & ghost_node,
                           const std::function<double(double, double)> & field)
{
    Box * n = qtree->get_cell_box(ghost_node.large_node);
    Box * m = qtree->get_cell_box(ghost_node.small_node);
    if (n->layer != m->layer) return field(m->x, m->y);
    if (std::fabs(m->y - n->y) > std::fabs(m->x - n->x)) return field(m->x, n->y);
    return field(n->x, m->y);
}

double linear_field(double x, double y)
{
    return 3.0 + 0.2 * x - 0.7 * y;
}

double linear_head(Box * box)
{
    return linear_field(box->x, box->y);
}

void expect_linear_head_reproduced(QuadTree3D * qtree, const GhostNode & ghost_node)
{
    EXPECT_NEAR(ghost_head(qtree, ghost_node, linear_head),
                ghost_location_head(qtree, ghost_node, linear_field), 1.0e-9)
        << "n " << ghost_node.large_node << " m " << ghost_node.small_node;
}

std::vector<GhostNode> vertical_ghost_nodes(QuadTree3D * qtree)
{
    std::vector<GhostNode> ghost_nodes, vertical;
    compute_ghost_nodes(qtree, ghost_nodes, true);
    for (const GhostNode & ghost_node : ghost_nodes)
    {
        if (is_vertical(qtree, ghost_node)) vertical.push_back(ghost_node);
    }
    return vertical;
}

} // namespace

TEST(GhostNode, RefinedCenterCell)
{
    // 3 x 3 cells of 100 with the center cell refined once; each of the four
    // side cells connects to two fine cells
    std::unique_ptr<ModflowGrid> mfgrid = make_mfgrid(1, 3, 3, 100.0, 100.0);
    std::unique_ptr<QuadTree3D> qtree = make_qtree(mfgrid.get(), {{150.0, 150.0, 1, 0}});
    CSRData csr = qtree->get_usg_csr_data();

    std::set<std::pair<int, int>> connected;
    for (int n = 0; n < csr.nodes; n++)
    {
        for (int ipos = csr.ia[n] + 1; ipos < csr.ia[n + 1]; ipos++)
        {
            connected.insert(std::make_pair(n, csr.ja[ipos]));
        }
    }

    std::vector<GhostNode> ghost_nodes;
    compute_ghost_nodes(qtree.get(), ghost_nodes);
    ASSERT_EQ(ghost_nodes.size(), 8u);

    for (const GhostNode & ghost_node : ghost_nodes)
    {
        // n is the coarse cell, m the connected fine cell
        EXPECT_EQ(csr.area[ghost_node.large_node], 10000.0);
        EXPECT_EQ(csr.area[ghost_node.small_node], 2500.0);
        EXPECT_TRUE(connected.count(
            std::make_pair((int)ghost_node.large_node, (int)ghost_node.small_node)));
        EXPECT_TRUE(
            connected.count(std::make_pair((int)ghost_node.large_node, (int)ghost_node.j1)));
        EXPECT_TRUE(
            connected.count(std::make_pair((int)ghost_node.large_node, (int)ghost_node.j2)));

        // one contributing cell, written twice with half the factor each
        EXPECT_EQ(ghost_node.j1, ghost_node.j2);
        EXPECT_DOUBLE_EQ(ghost_node.alpha1, 0.125);
        EXPECT_DOUBLE_EQ(ghost_node.alpha2, 0.125);
        expect_linear_head_reproduced(qtree.get(), ghost_node);
    }
}

TEST(GhostNode, TwoFinerNeighborsReproduceLinearHead)
{
    // 3 x 3 cells of 100 with the center and northeast cells refined; the
    // north cell connects to fine cells on its south and east sides, and the
    // neighbors on the other of those sides are two fine cells
    std::unique_ptr<QuadTree3D> qtree = make_qtree(make_mfgrid(1, 3, 3, 100.0, 100.0),
                                                   {{150.0, 150.0, 1, 0}, {250.0, 250.0, 1, 0}});
    std::set<std::pair<int, int>> connected = connections(qtree.get());

    std::vector<GhostNode> ghost_nodes;
    compute_ghost_nodes(qtree.get(), ghost_nodes);
    int two_cells = 0;
    for (const GhostNode & ghost_node : ghost_nodes)
    {
        EXPECT_TRUE(connected.count(
            std::make_pair((int)ghost_node.large_node, (int)ghost_node.small_node)));
        if (ghost_node.j1 != ghost_node.j2) two_cells++;
        expect_linear_head_reproduced(qtree.get(), ghost_node);
    }
    EXPECT_EQ(two_cells, 4);
}

TEST(GhostNode, NoVerticalCorrectionUnlessAsked)
{
    // the center of layer 2 refined; layer 1 is coarse above it
    std::unique_ptr<QuadTree3D> qtree =
        make_qtree(make_mfgrid(2, 3, 3, 100.0, 100.0), {{150.0, 150.0, 1, 1}});
    std::vector<GhostNode> horizontal, all;
    compute_ghost_nodes(qtree.get(), horizontal);
    compute_ghost_nodes(qtree.get(), all, true);

    ASSERT_EQ(horizontal.size(), 8u);
    for (const GhostNode & ghost_node : horizontal)
        EXPECT_FALSE(is_vertical(qtree.get(), ghost_node));

    // the vertical ghost nodes follow the horizontal ones, which are unchanged
    ASSERT_EQ(all.size(), 12u);
    for (size_t i = 0; i < horizontal.size(); i++)
    {
        EXPECT_EQ(all[i].large_node, horizontal[i].large_node);
        EXPECT_EQ(all[i].small_node, horizontal[i].small_node);
    }
}

TEST(GhostNode, NoVerticalCorrectionBetweenCellsOfOneSize)
{
    std::unique_ptr<QuadTree3D> qtree = make_qtree(make_mfgrid(2, 3, 3, 100.0, 100.0), {});
    std::vector<GhostNode> ghost_nodes;
    compute_ghost_nodes(qtree.get(), ghost_nodes, true);
    EXPECT_TRUE(ghost_nodes.empty());
}

TEST(GhostNode, VerticalCorrectionReproducesLinearHead)
{
    // the center cell of layer 1 is above four fine cells of layer 2; each
    // ghost node interpolates from the center cell and its neighbors on the
    // side of the fine cell in x and in y
    std::unique_ptr<QuadTree3D> qtree =
        make_qtree(make_mfgrid(2, 3, 3, 100.0, 100.0), {{150.0, 150.0, 1, 1}});
    std::set<std::pair<int, int>> connected = connections(qtree.get());

    std::vector<GhostNode> vertical = vertical_ghost_nodes(qtree.get());
    ASSERT_EQ(vertical.size(), 4u);
    for (const GhostNode & ghost_node : vertical)
    {
        Box * n = qtree->get_cell_box(ghost_node.large_node);
        EXPECT_EQ(n->layer, 0);
        EXPECT_EQ(n->x, 150.0);
        EXPECT_EQ(n->y, 150.0);
        int n_index = (int)ghost_node.large_node;
        EXPECT_TRUE(connected.count(std::make_pair(n_index, (int)ghost_node.small_node)));
        EXPECT_TRUE(connected.count(std::make_pair(n_index, (int)ghost_node.j1)));
        EXPECT_TRUE(connected.count(std::make_pair(n_index, (int)ghost_node.j2)));
        EXPECT_NE(ghost_node.j1, ghost_node.j2);
        // a quarter of the way to each neighbor
        EXPECT_DOUBLE_EQ(ghost_node.alpha1, 0.25);
        EXPECT_DOUBLE_EQ(ghost_node.alpha2, 0.25);
        expect_linear_head_reproduced(qtree.get(), ghost_node);
    }
}

TEST(GhostNode, VerticalCorrectionWithFinerAndCoarserNeighbors)
{
    // layer 1: the center and east cells refined; layer 2: the center and its
    // northeast quarter refined again, so each layer has coarse cells above or
    // below fine ones, with finer and coarser neighbors beside them
    std::unique_ptr<QuadTree3D> qtree = make_qtree(
        make_mfgrid(2, 4, 4, 100.0, 100.0),
        {{150.0, 250.0, 1, 0}, {250.0, 250.0, 1, 0}, {150.0, 250.0, 1, 1}, {175.0, 275.0, 2, 1}},
        true);
    std::set<std::pair<int, int>> connected = connections(qtree.get());

    std::vector<GhostNode> vertical = vertical_ghost_nodes(qtree.get());
    ASSERT_FALSE(vertical.empty());
    int two_cells = 0;
    int finer_neighbors = 0;
    int coarser_neighbors = 0;
    for (const GhostNode & ghost_node : vertical)
    {
        int n_index = (int)ghost_node.large_node;
        EXPECT_TRUE(connected.count(std::make_pair(n_index, (int)ghost_node.small_node)));
        EXPECT_TRUE(connected.count(std::make_pair(n_index, (int)ghost_node.j1)));
        EXPECT_TRUE(connected.count(std::make_pair(n_index, (int)ghost_node.j2)));
        EXPECT_LT(ghost_node.alpha1 + ghost_node.alpha2, 1.0);

        // a neighbor finer than n is the one of its pair on the side of m
        Box * n = qtree->get_cell_box(ghost_node.large_node);
        Box * m = qtree->get_cell_box(ghost_node.small_node);
        for (long j_node : {ghost_node.j1, ghost_node.j2})
        {
            Box * j = qtree->get_cell_box(j_node);
            if (j->depth < n->depth) coarser_neighbors++;
            if (j->depth <= n->depth) continue;
            finer_neighbors++;
            bool along_x = std::fabs(j->x - n->x) < std::fabs(j->y - n->y);
            double j_offset = along_x ? j->x - n->x : j->y - n->y;
            double m_offset = along_x ? m->x - n->x : m->y - n->y;
            EXPECT_EQ(j_offset > 0, m_offset > 0);
        }

        if (ghost_node.j1 == ghost_node.j2) continue;
        two_cells++;
        expect_linear_head_reproduced(qtree.get(), ghost_node);
    }
    EXPECT_GT(two_cells, 0);
    EXPECT_GT(finer_neighbors, 0);
    EXPECT_GT(coarser_neighbors, 0);
}

TEST(GhostNode, VerticalCorrectionAtTheGridEdge)
{
    // the northwest cell of layer 2 refined: below the northwest corner of
    // the cell above there is no neighbor in x or y, so no ghost node; below
    // its northeast and southwest corners only one direction contributes
    std::unique_ptr<QuadTree3D> qtree =
        make_qtree(make_mfgrid(2, 2, 2, 100.0, 100.0), {{50.0, 150.0, 1, 1}});

    std::vector<GhostNode> vertical = vertical_ghost_nodes(qtree.get());
    ASSERT_EQ(vertical.size(), 3u);
    int one_cell = 0;
    for (const GhostNode & ghost_node : vertical)
    {
        if (ghost_node.j1 != ghost_node.j2)
        {
            expect_linear_head_reproduced(qtree.get(), ghost_node);
            continue;
        }
        one_cell++;
        EXPECT_DOUBLE_EQ(ghost_node.alpha1, ghost_node.alpha2);
        EXPECT_DOUBLE_EQ(ghost_node.alpha1, 0.125);
        // reproduces a head that varies along the direction of the neighbor
        Box * n = qtree->get_cell_box(ghost_node.large_node);
        Box * j = qtree->get_cell_box(ghost_node.j1);
        auto along = [&](double x, double y) { return (j->y == n->y) ? 0.2 * x : -0.7 * y; };
        auto head = [&](Box * box) { return along(box->x, box->y); };
        EXPECT_NEAR(ghost_head(qtree.get(), ghost_node, head),
                    ghost_location_head(qtree.get(), ghost_node, along), 1.0e-9);
    }
    EXPECT_EQ(one_cell, 2);
}

TEST(GhostNode, VerticalCorrectionSkipsInactiveNeighbor)
{
    // as VerticalCorrectionReproducesLinearHead with the east cell of layer 1
    // inactive: below the east half of the center cell only y contributes
    std::unique_ptr<ModflowGrid> mfgrid = make_mfgrid(2, 3, 3, 100.0, 100.0);
    std::unique_ptr<QuadTree3D> qtree = make_qtree(mfgrid.get(), {{150.0, 150.0, 1, 1}});
    Box * east = qtree->node_group[mfgrid->get_box_index(0, 1, 2)];
    east->active = false;
    qtree->number_nodes();
    qtree->refresh();

    std::vector<GhostNode> vertical = vertical_ghost_nodes(qtree.get());
    ASSERT_EQ(vertical.size(), 4u);
    for (const GhostNode & ghost_node : vertical)
    {
        EXPECT_NE(ghost_node.j1, -1);
        EXPECT_NE(ghost_node.j2, -1);
        Box * m = qtree->get_cell_box(ghost_node.small_node);
        if (m->x < 150.0)
        {
            EXPECT_NE(ghost_node.j1, ghost_node.j2);
            continue;
        }
        EXPECT_EQ(ghost_node.j1, ghost_node.j2);
        Box * j = qtree->get_cell_box(ghost_node.j1);
        EXPECT_EQ(j->x, 150.0);
    }
}
