// Horizontal and vertical ghost node correction
#include <algorithm>
#include <cmath>
#include <functional>
#include <iterator>
#include <list>
#include <map>
#include <set>
#include <string>
#include <utility>
#include <vector>

#include <gtest/gtest.h>

#include "GhostNode.h"
#include "GridgenError.h"
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

// the ghost node of the connection between the cells centered at large and
// small, or nullptr
const GhostNode * find_ghost_node(QuadTree3D * qtree, const std::vector<GhostNode> & ghost_nodes,
                                  std::pair<double, double> large, std::pair<double, double> small)
{
    for (const GhostNode & ghost_node : ghost_nodes)
    {
        Box * n = qtree->get_cell_box(ghost_node.large_node);
        Box * m = qtree->get_cell_box(ghost_node.small_node);
        if (n->x == large.first && n->y == large.second && m->x == small.first &&
            m->y == small.second)
            return &ghost_node;
    }
    return nullptr;
}

// a ghost node interpolated from the coarser cell centered at contributor
// and a second cell, with positive factors that sum to less than one
void expect_coarser_contributor(QuadTree3D * qtree, const GhostNode & ghost_node,
                                std::pair<double, double> contributor)
{
    Box * j1 = qtree->get_cell_box(ghost_node.j1);
    EXPECT_EQ(j1->x, contributor.first);
    EXPECT_EQ(j1->y, contributor.second);
    EXPECT_NE(ghost_node.j1, ghost_node.j2);
    EXPECT_GT(ghost_node.alpha1, 0.0);
    EXPECT_GT(ghost_node.alpha2, 0.0);
    EXPECT_LT(ghost_node.alpha1 + ghost_node.alpha2, 1.0);
    expect_linear_head_reproduced(qtree, ghost_node);
}

// the leaf cell centered at (x, y) in layer, or nullptr
Box * cell_at(QuadTree3D * qtree, double x, double y, int layer)
{
    for (Box * box : qtree->node_group)
    {
        if (box->is_leaf && box->layer == layer && box->x == x && box->y == y) return box;
    }
    return nullptr;
}

// deactivate the leaf cells centered at the points of layer and renumber
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
}

bool reproduces_linear_head(QuadTree3D * qtree, const GhostNode & ghost_node)
{
    return std::fabs(ghost_head(qtree, ghost_node, linear_head) -
                     ghost_location_head(qtree, ghost_node, linear_field)) < 1.0e-9;
}

bool next_to_inactive_cell(Box * cell)
{
    for (int side = 0; side < 4; side++)
    {
        std::list<Box *> neighbors;
        cell->get_neighbors(neighbors, side);
        for (Box * neighbor : neighbors)
        {
            if (neighbor->number == -1) return true;
        }
    }
    return false;
}

// every contributing factor is not negative, each pair sums to less than one,
// and the contributing cells are active; returns the number of connections from
// a cell next to an inactive one to a smaller cell, horizontal and, if vertical
// is true, vertical, whose ghost node is missing or does not reproduce a linear
// head field
int expect_interpolated(QuadTree3D * qtree, const std::vector<GhostNode> & ghost_nodes,
                        bool vertical)
{
    std::map<std::pair<long, long>, const GhostNode *> by_connection;
    for (const GhostNode & ghost_node : ghost_nodes)
    {
        by_connection[std::make_pair(ghost_node.large_node, ghost_node.small_node)] = &ghost_node;
        EXPECT_NE(ghost_node.j1, -1);
        EXPECT_NE(ghost_node.j2, -1);
        EXPECT_GE(ghost_node.alpha1, 0.0) << "n " << ghost_node.large_node;
        EXPECT_GE(ghost_node.alpha2, 0.0) << "n " << ghost_node.large_node;
        EXPECT_LT(ghost_node.alpha1 + ghost_node.alpha2, 1.0);
    }

    int inexact_at_domain_edge = 0;
    CSRData csr = qtree->get_usg_csr_data();
    for (int n = 0; n < csr.nodes; n++)
    {
        Box * large = qtree->get_cell_box(n);
        if (!next_to_inactive_cell(large)) continue;
        for (int ipos = csr.ia[n] + 1; ipos < csr.ia[n + 1]; ipos++)
        {
            int m = csr.ja[ipos];
            if (csr.area[n] <= csr.area[m]) continue;
            if (!vertical && large->layer != qtree->get_cell_box(m)->layer) continue;
            auto found = by_connection.find(std::make_pair((long)n, (long)m));
            if (found == by_connection.end() || !reproduces_linear_head(qtree, *found->second))
                inexact_at_domain_edge++;
        }
    }
    return inexact_at_domain_edge;
}

// 8 x 8 cells of cell_size in 3 layers, refined and smoothed in each layer,
// with an irregular active domain in layer 2
std::unique_ptr<QuadTree3D> cross_domain_grid(double cell_size)
{
    double scale = cell_size / 100.0;
    std::unique_ptr<QuadTree3D> qtree = make_qtree(make_mfgrid(3, 8, 8, cell_size, cell_size),
                                                   {{330.0 * scale, 430.0 * scale, 3, 0},
                                                    {510.0 * scale, 270.0 * scale, 2, 0},
                                                    {250.0 * scale, 350.0 * scale, 2, 1},
                                                    {450.0 * scale, 450.0 * scale, 1, 1},
                                                    {410.0 * scale, 390.0 * scale, 3, 2},
                                                    {170.0 * scale, 610.0 * scale, 2, 2}},
                                                   true);
    std::vector<std::pair<double, double>> inactive;
    for (Box * box : qtree->node_group)
    {
        if (!box->is_leaf || box->layer != 1) continue;
        // outside a cross with a notch, so the domain has corners, inlets,
        // and cells with one active neighbor
        double x = box->x / scale;
        double y = box->y / scale;
        bool cross = (x > 200.0 && x < 600.0) || (y > 300.0 && y < 500.0);
        bool notch = x > 300.0 && x < 400.0 && y > 400.0 && y < 500.0;
        if (!cross || notch) inactive.push_back(std::make_pair(box->x, box->y));
    }
    deactivate(qtree.get(), inactive, 1);
    return qtree;
}

// the ghost nodes as the layer and the center, multiplied by scale and
// rounded, of the cells n, m, j1, and j2, with the factors rounded, sorted
std::vector<std::vector<long>>
records_by_location(QuadTree3D * qtree, const std::vector<GhostNode> & ghost_nodes, double scale)
{
    std::vector<std::vector<long>> records;
    for (const GhostNode & ghost_node : ghost_nodes)
    {
        std::vector<long> record;
        for (long node :
             {ghost_node.large_node, ghost_node.small_node, ghost_node.j1, ghost_node.j2})
        {
            Box * box = qtree->get_cell_box(node);
            record.push_back(box->layer);
            record.push_back(std::lround(box->x * scale * 1000.0));
            record.push_back(std::lround(box->y * scale * 1000.0));
        }
        record.push_back(std::lround(ghost_node.alpha1 * 1.0e9));
        record.push_back(std::lround(ghost_node.alpha2 * 1.0e9));
        records.push_back(record);
    }
    std::sort(records.begin(), records.end());
    return records;
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

TEST(GhostNode, CoarserNeighborOffsetTowardSmallCell)
{
    // 3 x 3 cells of 100 with the center cell refined twice in its northeast
    // quarter; the fine cell at (162.5, 187.5) connects west to the northwest
    // quarter, whose neighbor on the north is the coarse cell centered at
    // (150, 250), offset 25 toward the fine cell
    std::unique_ptr<QuadTree3D> qtree =
        make_qtree(make_mfgrid(1, 3, 3, 100.0, 100.0), {{162.5, 187.5, 2, 0}});
    std::set<std::pair<int, int>> connected = connections(qtree.get());
    std::vector<GhostNode> ghost_nodes;
    compute_ghost_nodes(qtree.get(), ghost_nodes);

    const GhostNode * ghost_node =
        find_ghost_node(qtree.get(), ghost_nodes, {125.0, 175.0}, {162.5, 187.5});
    ASSERT_NE(ghost_node, nullptr);
    expect_coarser_contributor(qtree.get(), *ghost_node, {150.0, 250.0});
    EXPECT_TRUE(connected.count(std::make_pair((int)ghost_node->large_node, (int)ghost_node->j2)));
}

TEST(GhostNode, CoarserNeighborOffsetAwayFromSmallCell)
{
    // 3 x 3 cells of 100 with the center cell refined once and the east cell
    // refined twice in its northwest quarter; the fine cell at (212.5, 187.5)
    // connects west to the northeast quarter of the center cell, whose
    // neighbor on the north is the coarse cell centered at (150, 250), offset
    // 25 away from the fine cell
    std::unique_ptr<QuadTree3D> qtree = make_qtree(make_mfgrid(1, 3, 3, 100.0, 100.0),
                                                   {{150.0, 150.0, 1, 0}, {212.5, 187.5, 2, 0}});
    std::set<std::pair<int, int>> connected = connections(qtree.get());
    std::vector<GhostNode> ghost_nodes;
    compute_ghost_nodes(qtree.get(), ghost_nodes);

    const GhostNode * ghost_node =
        find_ghost_node(qtree.get(), ghost_nodes, {175.0, 175.0}, {212.5, 187.5});
    ASSERT_NE(ghost_node, nullptr);
    expect_coarser_contributor(qtree.get(), *ghost_node, {150.0, 250.0});
    EXPECT_TRUE(connected.count(std::make_pair((int)ghost_node->large_node, (int)ghost_node->j2)));
}

TEST(GhostNode, CoarserNeighborWithoutSecondCellUsesTheCellBeyondIt)
{
    // the grid of CoarserNeighborOffsetTowardSmallCell with the west cell
    // inactive; no neighbor of the northwest quarter, with the coarse cell on
    // its north, interpolates the head at the ghost node, but the cell west of
    // the coarse cell, which is not connected to the quarter, does, a fourth
    // of the way from it (#170)
    std::unique_ptr<ModflowGrid> mfgrid = make_mfgrid(1, 3, 3, 100.0, 100.0);
    std::unique_ptr<QuadTree3D> qtree = make_qtree(mfgrid.get(), {{162.5, 187.5, 2, 0}});
    qtree->node_group[mfgrid->get_box_index(0, 1, 0)]->active = false;
    qtree->number_nodes();
    qtree->refresh();
    std::set<std::pair<int, int>> connected = connections(qtree.get());
    std::vector<GhostNode> ghost_nodes;
    int inexact = compute_ghost_nodes(qtree.get(), ghost_nodes);

    const GhostNode * ghost_node =
        find_ghost_node(qtree.get(), ghost_nodes, {125.0, 175.0}, {162.5, 187.5});
    ASSERT_NE(ghost_node, nullptr);
    Box * j1 = qtree->get_cell_box(ghost_node->j1);
    Box * j2 = qtree->get_cell_box(ghost_node->j2);
    EXPECT_EQ(j1->x, 150.0);
    EXPECT_EQ(j1->y, 250.0);
    EXPECT_EQ(j2->x, 50.0);
    EXPECT_EQ(j2->y, 250.0);
    EXPECT_FALSE(connected.count(std::make_pair((int)ghost_node->large_node, (int)ghost_node->j2)));
    EXPECT_DOUBLE_EQ(ghost_node->alpha1, 0.125);
    EXPECT_DOUBLE_EQ(ghost_node->alpha2, 1.0 / 24.0);
    expect_linear_head_reproduced(qtree.get(), *ghost_node);
    EXPECT_EQ(expect_interpolated(qtree.get(), ghost_nodes, false), inexact);
}

TEST(GhostNode, EveryHorizontalGhostNodeReproducesLinearHead)
{
    // 6 x 6 cells of 100 refined by up to three levels at scattered points and
    // smoothed, so cells have coarser neighbors offset along their
    // connections
    std::unique_ptr<QuadTree3D> qtree = make_qtree(make_mfgrid(1, 6, 6, 100.0, 100.0),
                                                   {{130.0, 420.0, 3, 0},
                                                    {310.0, 290.0, 2, 0},
                                                    {470.0, 120.0, 3, 0},
                                                    {250.0, 250.0, 1, 0},
                                                    {90.0, 80.0, 2, 0}},
                                                   true);
    std::set<std::pair<int, int>> connected = connections(qtree.get());
    std::vector<GhostNode> ghost_nodes;
    compute_ghost_nodes(qtree.get(), ghost_nodes);
    ASSERT_FALSE(ghost_nodes.empty());

    int coarser = 0;
    for (const GhostNode & ghost_node : ghost_nodes)
    {
        int n = (int)ghost_node.large_node;
        EXPECT_TRUE(connected.count(std::make_pair(n, (int)ghost_node.j1)));
        EXPECT_TRUE(connected.count(std::make_pair(n, (int)ghost_node.j2)));
        EXPECT_GT(ghost_node.alpha1, 0.0);
        EXPECT_GT(ghost_node.alpha2, 0.0);
        EXPECT_LT(ghost_node.alpha1 + ghost_node.alpha2, 1.0);
        expect_linear_head_reproduced(qtree.get(), ghost_node);
        Box * large = qtree->get_cell_box(ghost_node.large_node);
        Box * j1 = qtree->get_cell_box(ghost_node.j1);
        if (j1->depth < large->depth) coarser++;
    }
    EXPECT_GT(coarser, 0);
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
    // inactive: below the west half of the center cell the ghost nodes are
    // interpolated as before, and below the east half from the north or south
    // cell and the cell beyond it to the east, which is not connected to the
    // center cell (#170)
    std::unique_ptr<ModflowGrid> mfgrid = make_mfgrid(2, 3, 3, 100.0, 100.0);
    std::unique_ptr<QuadTree3D> qtree = make_qtree(mfgrid.get(), {{150.0, 150.0, 1, 1}});
    Box * east = qtree->node_group[mfgrid->get_box_index(0, 1, 2)];
    east->active = false;
    qtree->number_nodes();
    qtree->refresh();
    std::set<std::pair<int, int>> connected = connections(qtree.get());

    std::vector<GhostNode> vertical = vertical_ghost_nodes(qtree.get());
    ASSERT_EQ(vertical.size(), 4u);
    for (const GhostNode & ghost_node : vertical)
    {
        EXPECT_NE(ghost_node.j1, ghost_node.j2);
        expect_linear_head_reproduced(qtree.get(), ghost_node);
        Box * m = qtree->get_cell_box(ghost_node.small_node);
        int n_index = (int)ghost_node.large_node;
        bool both_connected = connected.count(std::make_pair(n_index, (int)ghost_node.j1)) &&
                              connected.count(std::make_pair(n_index, (int)ghost_node.j2));
        EXPECT_EQ(both_connected, m->x < 150.0);
    }
    std::vector<GhostNode> ghost_nodes;
    EXPECT_EQ(compute_ghost_nodes(qtree.get(), ghost_nodes, true), 0);
    EXPECT_EQ(expect_interpolated(qtree.get(), ghost_nodes, true), 0);
}

TEST(GhostNode, CoarserNeighborWithoutSecondCellInterpolatesFromSmallCell)
{
    // the grid of CoarserNeighborOffsetAwayFromSmallCell with the fine cell at
    // (212.5, 187.5)'s neighbor (212.5, 162.5) inactive; the ghost node is
    // interpolated from the coarse cell on the north and the fine cell m of the
    // connection itself (#170)
    std::unique_ptr<QuadTree3D> qtree = make_qtree(make_mfgrid(1, 3, 3, 100.0, 100.0),
                                                   {{150.0, 150.0, 1, 0}, {212.5, 187.5, 2, 0}});
    deactivate(qtree.get(), {{212.5, 162.5}}, 0);
    std::vector<GhostNode> ghost_nodes;
    int inexact = compute_ghost_nodes(qtree.get(), ghost_nodes);

    const GhostNode * ghost_node =
        find_ghost_node(qtree.get(), ghost_nodes, {175.0, 175.0}, {212.5, 187.5});
    ASSERT_NE(ghost_node, nullptr);
    Box * j1 = qtree->get_cell_box(ghost_node->j1);
    Box * j2 = qtree->get_cell_box(ghost_node->j2);
    std::set<std::pair<double, double>> contributors = {{j1->x, j1->y}, {j2->x, j2->y}};
    EXPECT_EQ(contributors, (std::set<std::pair<double, double>>{{150.0, 250.0}, {212.5, 187.5}}));
    expect_linear_head_reproduced(qtree.get(), *ghost_node);
    EXPECT_EQ(expect_interpolated(qtree.get(), ghost_nodes, false), inexact);
}

TEST(GhostNode, HeadBeyondTheActiveCellsIsNotProjected)
{
    // 3 x 3 cells of 100 in 2 layers with the center of layer 1 refined and,
    // in layer 2, the north, east, and south cells inactive: below the west
    // half of the center cell the ghost nodes are interpolated from the west
    // cell and the cells beside it, but below the east half the head could
    // only be projected beyond the active cells, which takes a negative
    // contributing factor, so those two ghost nodes are not written and are
    // counted (#170)
    std::unique_ptr<QuadTree3D> qtree =
        make_qtree(make_mfgrid(2, 3, 3, 100.0, 100.0), {{150.0, 150.0, 1, 0}});
    deactivate(qtree.get(), {{150.0, 250.0}, {250.0, 150.0}, {150.0, 50.0}}, 1);

    std::vector<GhostNode> ghost_nodes;
    int inexact = compute_ghost_nodes(qtree.get(), ghost_nodes, true);
    EXPECT_EQ(inexact, 2);
    EXPECT_EQ(expect_interpolated(qtree.get(), ghost_nodes, true), inexact);
    int below_west_half = 0;
    for (const GhostNode & ghost_node : ghost_nodes)
    {
        if (!is_vertical(qtree.get(), ghost_node)) continue;
        Box * m = qtree->get_cell_box(ghost_node.small_node);
        EXPECT_LT(m->x, 150.0);
        expect_linear_head_reproduced(qtree.get(), ghost_node);
        below_west_half++;
    }
    EXPECT_EQ(below_west_half, 2);
}

TEST(GhostNode, NoNegativeFactorAroundActiveDomains)
{
    // 8 x 8 cells of 100 in 3 layers, refined and smoothed in each layer, with
    // an irregular active domain in layer 2: no contributing factor is
    // negative, and the connections at the edge of the active domain whose
    // ghost node is missing or inexact are the ones counted (#170)
    std::unique_ptr<QuadTree3D> qtree = cross_domain_grid(100.0);
    for (bool pass_through : {false, true})
    {
        SCOPED_TRACE(pass_through ? "pass-through" : "no pass-through");
        qtree->vertical_pass_through = pass_through;
        std::vector<GhostNode> ghost_nodes;
        int inexact = compute_ghost_nodes(qtree.get(), ghost_nodes, true);
        ASSERT_FALSE(ghost_nodes.empty());
        EXPECT_EQ(expect_interpolated(qtree.get(), ghost_nodes, true), inexact);
    }
}

TEST(GhostNode, GhostNodesDoNotDependOnCellSize)
{
    // the grid of NoNegativeFactorAroundActiveDomains with cells of 100 / 3,
    // whose cell centers are not exact in floating point: the same cells
    // contribute with the same factors, and the same connections are counted,
    // so a contributing factor of zero is not rejected as negative
    std::unique_ptr<QuadTree3D> reference = cross_domain_grid(100.0);
    std::unique_ptr<QuadTree3D> scaled = cross_domain_grid(100.0 / 3.0);
    std::vector<GhostNode> reference_nodes, scaled_nodes;
    int reference_inexact = compute_ghost_nodes(reference.get(), reference_nodes, true);
    int scaled_inexact = compute_ghost_nodes(scaled.get(), scaled_nodes, true);
    EXPECT_EQ(scaled_inexact, reference_inexact);
    EXPECT_EQ(expect_interpolated(scaled.get(), scaled_nodes, true), scaled_inexact);
    std::vector<std::vector<long>> scaled_records =
        records_by_location(scaled.get(), scaled_nodes, 3.0);
    std::vector<std::vector<long>> reference_records =
        records_by_location(reference.get(), reference_nodes, 1.0);
    std::vector<std::vector<long>> differ;
    std::set_symmetric_difference(scaled_records.begin(), scaled_records.end(),
                                  reference_records.begin(), reference_records.end(),
                                  std::back_inserter(differ));
    for (const std::vector<long> & record : differ)
    {
        std::string text;
        for (long value : record) text += std::to_string(value) + " ";
        ADD_FAILURE() << "in one grid only: " << text;
    }
}

TEST(GhostNode, NoConnectionCountedWithoutAnActiveDomain)
{
    // every cell active: nothing is counted, and no factor is negative
    std::unique_ptr<QuadTree3D> qtree = make_qtree(
        make_mfgrid(2, 5, 5, 100.0, 100.0), {{250.0, 250.0, 2, 0}, {50.0, 450.0, 1, 1}}, true);
    std::vector<GhostNode> ghost_nodes;
    EXPECT_EQ(compute_ghost_nodes(qtree.get(), ghost_nodes, true), 0);
    EXPECT_EQ(expect_interpolated(qtree.get(), ghost_nodes, true), 0);
}

TEST(GhostNode, VerticalCorrectionAcrossPassThrough)
{
    // 3 layers of 3 x 3 cells with the center cell of layer 3 refined and the
    // center cell of layer 2 inactive; with pass-through the center cell of
    // layer 1 connects to the four fine cells of layer 3 and gets their ghost
    // nodes, and without it there are none
    std::unique_ptr<ModflowGrid> mfgrid = make_mfgrid(3, 3, 3, 100.0, 100.0);
    std::unique_ptr<QuadTree3D> qtree = make_qtree(mfgrid.get(), {{150.0, 150.0, 1, 2}});
    qtree->node_group[mfgrid->get_box_index(1, 1, 1)]->active = false;
    qtree->number_nodes();
    qtree->refresh();

    std::vector<GhostNode> vertical = vertical_ghost_nodes(qtree.get());
    EXPECT_TRUE(vertical.empty());

    qtree->vertical_pass_through = true;
    std::set<std::pair<int, int>> connected = connections(qtree.get());
    vertical = vertical_ghost_nodes(qtree.get());
    ASSERT_EQ(vertical.size(), 4u);
    for (const GhostNode & ghost_node : vertical)
    {
        Box * n = qtree->get_cell_box(ghost_node.large_node);
        Box * m = qtree->get_cell_box(ghost_node.small_node);
        EXPECT_EQ(n->layer, 0);
        EXPECT_EQ(m->layer, 2);
        int n_index = (int)ghost_node.large_node;
        EXPECT_TRUE(connected.count(std::make_pair(n_index, (int)ghost_node.small_node)));
        EXPECT_TRUE(connected.count(std::make_pair(n_index, (int)ghost_node.j1)));
        EXPECT_TRUE(connected.count(std::make_pair(n_index, (int)ghost_node.j2)));
        expect_linear_head_reproduced(qtree.get(), ghost_node);
    }
}

TEST(GhostNode, NoVerticalCorrectionBetweenCellsThatDoNotOverlap)
{
    // 2 x 2 cells of 100 in 3 layers: the northwest cell of layer 1 refined
    // once, the northwest cell of layer 2 inactive, and layer 3 refined twice;
    // pass-through connects each fine cell of layer 1 only to the 4 cells of
    // layer 3 below it, and each of those connections gets a ghost node
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

    int across = 0;
    for (const GhostNode & ghost_node : vertical_ghost_nodes(qtree.get()))
    {
        Box * n = qtree->get_cell_box(ghost_node.large_node);
        Box * m = qtree->get_cell_box(ghost_node.small_node);
        EXPECT_TRUE(n->in(m->x, m->y))
            << "n " << ghost_node.large_node << " m " << ghost_node.small_node;
        EXPECT_LT(ghost_node.alpha1 + ghost_node.alpha2, 1.0);
        if (n->layer == 0) across++;
    }
    // each fine cell of layer 1 is above 4 cells of layer 3, one of which is
    // at the outer corner of the grid with no neighbor in x or y
    EXPECT_EQ(across, 15);
}

TEST(GhostNode, VerticalCorrectionWithMoreThanTwoNeighborsOnASideThrows)
{
    // layer 1: the east cell refined twice, so the center cell has four
    // neighbors on its east side; layer 2: the center cell refined once. The
    // horizontal ghost nodes do not use that side, the vertical ones do
    std::unique_ptr<QuadTree3D> qtree = make_qtree(make_mfgrid(2, 3, 3, 100.0, 100.0),
                                                   {{250.0, 150.0, 2, 0}, {150.0, 150.0, 1, 1}});
    std::vector<GhostNode> ghost_nodes;
    EXPECT_NO_THROW(compute_ghost_nodes(qtree.get(), ghost_nodes));
    ghost_nodes.clear();
    try
    {
        compute_ghost_nodes(qtree.get(), ghost_nodes, true);
        ADD_FAILURE() << "no error";
    }
    catch (const GridgenError & error)
    {
        EXPECT_NE(std::string(error.what()).find("more than two neighbors on one side"),
                  std::string::npos)
            << error.what();
    }
}
