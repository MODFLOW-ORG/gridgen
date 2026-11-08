// Horizontal ghost node correction
#include <set>
#include <utility>
#include <vector>

#include <gtest/gtest.h>

#include "GhostNode.h"
#include "grid_helpers.h"

using namespace cusg;
using gridgen_test::make_mfgrid;
using gridgen_test::make_qtree;

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

        // characterizes the current weights (half the offset of m from n over
        // the n-j distance), which #42 reviews
        EXPECT_DOUBLE_EQ(ghost_node.alpha1, 0.125);
        EXPECT_DOUBLE_EQ(ghost_node.alpha2, 0.125);
    }
}
