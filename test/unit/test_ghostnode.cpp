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
    ModflowGrid * mfgrid = make_mfgrid(1, 3, 3, 100.0, 100.0);
    QuadTree3D * qtree = make_qtree(mfgrid, {{150.0, 150.0, 1, 0}});
    CSRData csr = qtree->get_usg_csr_data();

    std::set<std::pair<int, int>> connected;
    for (int n = 0; n < csr.nodes; n++)
    {
        for (int ipos = csr.ia[n] + 1; ipos < csr.ia[n + 1]; ipos++)
        {
            connected.insert(std::make_pair(n, csr.ja[ipos]));
        }
    }

    std::vector<GhostNode> gns;
    computeGNC(qtree, gns);
    ASSERT_EQ(gns.size(), 8u);

    for (const GhostNode & gn : gns)
    {
        // n is the coarse cell, m the connected fine cell
        EXPECT_EQ(csr.area[gn.nNum], 10000.0);
        EXPECT_EQ(csr.area[gn.mNum], 2500.0);
        EXPECT_TRUE(connected.count(std::make_pair((int)gn.nNum, (int)gn.mNum)));
        EXPECT_TRUE(connected.count(std::make_pair((int)gn.nNum, (int)gn.j1)));
        EXPECT_TRUE(connected.count(std::make_pair((int)gn.nNum, (int)gn.j2)));

        // current weights: half the offset of m from n over the n-j distance
        EXPECT_DOUBLE_EQ(gn.apha1, 0.125);
        EXPECT_DOUBLE_EQ(gn.apha2, 0.125);
    }
    csr.destory();
}
