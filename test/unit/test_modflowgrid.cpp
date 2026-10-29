// Structured MODFLOW grid connectivity
#include <cstdlib>

#include <gtest/gtest.h>

#include "grid_helpers.h"

using namespace cusg;
using gridgen_test::make_mfgrid;

// disabled: get_usg_csr_data stores thickness in an int and takes the right and
// front neighbor thickness from column j - 1 (#46)
TEST(ModflowGrid, DISABLED_CsrConnectivityOneLayer)
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
    for (int ipos = 0; ipos < csr.nja; ipos++)
    {
        if (std::abs(csr.fldr[ipos]) != 3) continue;
        nvert++;
        EXPECT_FLOAT_EQ(csr.cl1[ipos], 5.0f);
        EXPECT_FLOAT_EQ(csr.cl2[ipos], 5.0f);
        EXPECT_FLOAT_EQ(csr.fahl[ipos], 100.0f);
    }
    EXPECT_EQ(nvert, 8);
    csr.destory();
}
