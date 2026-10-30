// Structured MODFLOW grid connectivity
#include <cstdlib>
#include <list>
#include <string>
#include <utility>
#include <vector>

#include <gtest/gtest.h>

#include "grid_helpers.h"

using namespace cusg;
using gridgen_test::make_mfgrid;

TEST(ModflowGrid, CsrConnectivityOneLayer)
{
    // 2 rows x 3 columns, delr 10, delc 20, 7.5 thick
    ModflowGrid * mfgrid = make_mfgrid(1, 2, 3, 10.0, 20.0, 7.5);
    CSRData csr = mfgrid->get_usg_csr_data();

    // 6 diagonal entries plus 7 connections stored in both directions
    EXPECT_EQ(csr.nodes, 6);
    EXPECT_EQ(csr.nja, 20);

    for (int n = 0; n < csr.nodes; n++)
    {
        EXPECT_DOUBLE_EQ(csr.area[n], 200.0);
    }

    for (int ipos = 0; ipos < csr.nja; ipos++)
    {
        int dir = std::abs(csr.fldr[ipos]);
        if (dir == 0) continue;
        if (dir == 1)
        {
            EXPECT_FLOAT_EQ(csr.cl1[ipos], 5.0f);
            EXPECT_FLOAT_EQ(csr.cl2[ipos], 5.0f);
            EXPECT_FLOAT_EQ(csr.fahl[ipos], 150.0f);
        }
        else if (dir == 2)
        {
            EXPECT_FLOAT_EQ(csr.cl1[ipos], 10.0f);
            EXPECT_FLOAT_EQ(csr.cl2[ipos], 10.0f);
            EXPECT_FLOAT_EQ(csr.fahl[ipos], 75.0f);
        }
        else
        {
            ADD_FAILURE() << "unexpected fldr " << csr.fldr[ipos];
        }
    }
    csr.destory();
}

TEST(ModflowGrid, CsrVerticalConnections)
{
    // 2 layers of 2 x 2 cells; each cell connects to the cell below
    ModflowGrid * mfgrid = make_mfgrid(2, 2, 2, 10.0, 10.0);
    CSRData csr = mfgrid->get_usg_csr_data();

    EXPECT_EQ(csr.nodes, 8);

    int nvert = 0;
    for (int n = 0; n < csr.nodes; n++)
    {
        for (int ipos = csr.ia[n] - 1; ipos < csr.ia[n + 1] - 1; ipos++)
        {
            if (std::abs(csr.fldr[ipos]) != 3) continue;
            nvert++;
            EXPECT_FLOAT_EQ(csr.cl1[ipos], 5.0f);
            EXPECT_FLOAT_EQ(csr.cl2[ipos], 5.0f);
            EXPECT_FLOAT_EQ(csr.fahl[ipos], 100.0f);
            // up is +3 and down is -3, as in the quadtree export
            EXPECT_EQ(csr.fldr[ipos], csr.ja[ipos] - 1 < n ? 3 : -3);
        }
    }
    EXPECT_EQ(nvert, 8);
    csr.destory();
}

TEST(ModflowGrid, CsrFaceAreaUsesBothCellThicknesses)
{
    // 2 layers of 2 x 2 cells, delr 10, delc 20, a different thickness in
    // each cell (#46)
    std::vector<double> thk = {4.5, 6.0, 8.0, 10.5, 3.0, 7.5, 5.0, 2.5};
    double top[4] = {0.0, 0.0, 0.0, 0.0};
    double bot1[4];
    double bot2[4];
    for (int n = 0; n < 4; n++)
    {
        bot1[n] = top[n] - thk[n];
        bot2[n] = bot1[n] - thk[n + 4];
    }
    std::vector<double *> surfaces = {top, bot1, bot2};
    ModflowGrid mfgrid(2, 2, 2, 10.0, 20.0, surfaces, 0.0, 0.0, 0.0);
    CSRData csr = mfgrid.get_usg_csr_data();

    // one-based ia and ja, with ia ending at nja + 1
    ASSERT_EQ(csr.nodes, 8);
    EXPECT_EQ(csr.ia[0], 1);
    EXPECT_EQ(csr.ia[csr.nodes], csr.nja + 1);

    for (int n = 0; n < csr.nodes; n++)
    {
        for (int ipos = csr.ia[n]; ipos < csr.ia[n + 1] - 1; ipos++)
        {
            int m = csr.ja[ipos] - 1;
            SCOPED_TRACE("node " + std::to_string(n + 1) + " to " +
                         std::to_string(m + 1));
            int dir = std::abs(csr.fldr[ipos]);
            if (dir == 3)
            {
                EXPECT_FLOAT_EQ(csr.cl1[ipos], 0.5 * thk[n]);
                EXPECT_FLOAT_EQ(csr.cl2[ipos], 0.5 * thk[m]);
                EXPECT_FLOAT_EQ(csr.fahl[ipos], 200.0f);
                continue;
            }
            // a column connection crosses delr through a face delc wide, a
            // row connection the reverse
            double length = dir == 1 ? 10.0 : 20.0;
            double width = dir == 1 ? 20.0 : 10.0;
            EXPECT_FLOAT_EQ(csr.cl1[ipos], 0.5 * length);
            EXPECT_FLOAT_EQ(csr.cl2[ipos], 0.5 * length);
            EXPECT_FLOAT_EQ(csr.fahl[ipos], width * 0.5 * (thk[n] + thk[m]));
        }
    }
    csr.destory();
}

TEST(ModflowGrid, NodeConnectionsMatchCsr)
{
    // 2 layers of 3 rows x 4 columns; the node file counts these connections
    ModflowGrid * mfgrid = make_mfgrid(2, 3, 4, 10.0, 10.0);
    CSRData csr = mfgrid->get_usg_csr_data();

    for (int n = 0; n < csr.nodes; n++)
    {
        SCOPED_TRACE("node " + std::to_string(n + 1));
        std::list<std::pair<int, int>> con = mfgrid->get_node_connections(n);
        ASSERT_EQ(static_cast<int>(con.size()), csr.iac[n] - 1);
        int ipos = csr.ia[n];
        for (const std::pair<int, int> & c : con)
        {
            EXPECT_EQ(c.first, csr.ja[ipos] - 1);
            EXPECT_EQ(c.second, csr.fldr[ipos]);
            ipos++;
        }
    }
    csr.destory();
}
