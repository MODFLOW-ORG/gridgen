// CSR consistency checks shared by the unit tests
#pragma once

#include <cstdlib>
#include <map>
#include <utility>

#include <gtest/gtest.h>

#include "QuadTree3D.h"

namespace gridgen_test
{

using namespace cusg;

// iac matches ia, the diagonal is first, and the off-diagonal ja are sorted
inline void expect_csr_valid(const CSRData & csr)
{
    for (int n = 0; n < csr.nodes; n++)
    {
        EXPECT_EQ(csr.iac[n], csr.ia[n + 1] - csr.ia[n]) << "node " << n;
        EXPECT_EQ(csr.ja[csr.ia[n]], n) << "node " << n;
        for (int ipos = csr.ia[n] + 2; ipos < csr.ia[n + 1]; ipos++)
        {
            EXPECT_LT(csr.ja[ipos - 1], csr.ja[ipos]) << "node " << n;
        }
    }
    EXPECT_EQ(csr.ia[csr.nodes], csr.nja);
}

// every n-m connection has an m-n connection with the same face area and
// direction, and swapped connection lengths
inline void expect_symmetric(const CSRData & csr)
{
    std::map<std::pair<int, int>, int> pos;
    for (int n = 0; n < csr.nodes; n++)
    {
        for (int ipos = csr.ia[n] + 1; ipos < csr.ia[n + 1]; ipos++)
        {
            pos[std::make_pair(n, csr.ja[ipos])] = ipos;
        }
    }
    for (const auto & entry : pos)
    {
        int n = entry.first.first;
        int m = entry.first.second;
        auto it = pos.find(std::make_pair(m, n));
        ASSERT_NE(it, pos.end()) << "no connection " << m << "-" << n;
        int ipos = entry.second;
        int jpos = it->second;
        EXPECT_EQ(csr.fldr[ipos], -csr.fldr[jpos]);
        EXPECT_FLOAT_EQ(csr.fahl[ipos], csr.fahl[jpos]);
        EXPECT_FLOAT_EQ(csr.cl1[ipos], csr.cl2[jpos]);
        EXPECT_FLOAT_EQ(csr.cl2[ipos], csr.cl1[jpos]);
    }
}

// on a fully active grid of thickness thk, the faces on each side of a cell
// that has neighbors on that side add up to the side length times thk (or
// the cell area for vertical faces)
inline void expect_faces_cover_sides(QuadTree3D * qtree, const CSRData & csr, double thk)
{
    for (int n = 0; n < csr.nodes; n++)
    {
        Box * box = qtree->get_nodeobj(n);
        std::map<int, double> total;
        for (int ipos = csr.ia[n] + 1; ipos < csr.ia[n + 1]; ipos++)
        {
            total[csr.fldr[ipos]] += csr.fahl[ipos];
        }
        for (const auto & side : total)
        {
            int dir = std::abs(side.first);
            double expected = (dir == 1)   ? box->dy * thk
                              : (dir == 2) ? box->dx * thk
                                           : box->dx * box->dy;
            EXPECT_NEAR(side.second, expected, 1e-4 * expected)
                << "node " << n << " fldr " << side.first;
        }
    }
}

// two CSR structures hold the same values
inline void expect_csr_equal(const CSRData & a, const CSRData & b)
{
    ASSERT_EQ(a.nodes, b.nodes);
    ASSERT_EQ(a.nja, b.nja);
    for (int n = 0; n <= a.nodes; n++) EXPECT_EQ(a.ia[n], b.ia[n]);
    for (int n = 0; n < a.nodes; n++) EXPECT_EQ(a.area[n], b.area[n]);
    for (int i = 0; i < a.nja; i++)
    {
        EXPECT_EQ(a.ja[i], b.ja[i]);
        EXPECT_EQ(a.fldr[i], b.fldr[i]);
        EXPECT_EQ(a.cl1[i], b.cl1[i]);
        EXPECT_EQ(a.cl2[i], b.cl2[i]);
        EXPECT_EQ(a.fahl[i], b.fahl[i]);
    }
}

} // namespace gridgen_test
