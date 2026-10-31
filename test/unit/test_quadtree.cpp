// Quadtree refinement, smoothing, and connectivity
#include <cstdlib>

#include <gtest/gtest.h>

#include "csr_checks.h"
#include "grid_helpers.h"

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
            int dn = qtree->get_nodeobj(n)->depth;
            int dm = qtree->get_nodeobj(csr.ja[ipos])->depth;
            jump = std::max(jump, std::abs(dn - dm));
        }
    }
    return jump;
}

} // namespace

TEST(QuadTree, RefineCenterCell)
{
    // 3 x 3 cells of 100; refine the center cell once
    ModflowGrid * mfgrid = make_mfgrid(1, 3, 3, 100.0, 100.0);
    QuadTree3D * qtree = make_qtree(mfgrid, {{150.0, 150.0, 1, 0}});
    CSRData csr = qtree->get_usg_csr_data();

    EXPECT_EQ(csr.nodes, 12);
    expect_faces_cover_sides(qtree, csr, 10.0);

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
    csr.destory();
}

TEST(QuadTree, SmoothingLimitsHorizontalLevelJump)
{
    QuadTree3D * rough =
        make_qtree(make_mfgrid(1, 5, 5, 100.0, 100.0), {{250.0, 250.0, 3, 0}}, false);
    CSRData csr = rough->get_usg_csr_data();
    EXPECT_GT(max_level_jump(rough, csr, false), 1);
    csr.destory();

    QuadTree3D * smooth =
        make_qtree(make_mfgrid(1, 5, 5, 100.0, 100.0), {{250.0, 250.0, 3, 0}}, true);
    csr = smooth->get_usg_csr_data();
    EXPECT_EQ(max_level_jump(smooth, csr, false), 1);
    expect_faces_cover_sides(smooth, csr, 10.0);
    expect_symmetric(csr);
    expect_csr_valid(csr);
    csr.destory();
}

TEST(QuadTree, SmoothingLimitsVerticalLevelJump)
{
    // refine only the top layer; smoothing propagates refinement downward
    ModflowGrid * mfgrid = make_mfgrid(3, 4, 4, 100.0, 100.0);
    QuadTree3D * qtree = make_qtree(mfgrid, {{150.0, 150.0, 3, 0}}, true);
    CSRData csr = qtree->get_usg_csr_data();

    EXPECT_EQ(max_level_jump(qtree, csr, true), 1);
    expect_faces_cover_sides(qtree, csr, 10.0);
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
    csr.destory();
}

TEST(QuadTree, VerticalPassThroughInactiveCell)
{
    // 3 layers of 3 x 3 cells; the center cell of layer 2 is inactive, so
    // with pass-through the center cells of layers 1 and 3 connect
    QuadTree3D * qtree = make_qtree(make_mfgrid(3, 3, 3, 100.0, 100.0), {});
    Box * center = nullptr;
    for (int n = 0; n < qtree->nodes(); n++)
    {
        Box * box = qtree->get_nodeobj(n);
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
            Box * bn = qtree->get_nodeobj(n);
            Box * bm = qtree->get_nodeobj(csr.ja[ipos]);
            if (std::abs(bn->layer - bm->layer) == 2) nlong++;
        }
    }
    // layer 1 to layer 3 and back
    EXPECT_EQ(nlong, 2);
    csr.destory();
}
