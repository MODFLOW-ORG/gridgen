// Box split and rotation
#include <cmath>

#include <gtest/gtest.h>

#include "Box.h"

using namespace cusg;

TEST(Box, SplitCreatesFourQuarterChildren)
{
    Box box(0.0, 0.0, 8.0, 4.0);
    ASSERT_TRUE(box.split(1e-10));
    EXPECT_FALSE(box.isLeaf);

    // children are numbered NW, NE, SE, SW
    const double xc[4] = {-2.0, 2.0, 2.0, -2.0};
    const double yc[4] = {1.0, 1.0, -1.0, -1.0};
    for (int i = 0; i < 4; i++)
    {
        Box * child = box.pChildren[i];
        ASSERT_NE(child, nullptr);
        EXPECT_TRUE(child->isLeaf);
        EXPECT_EQ(child->depth, 1);
        EXPECT_DOUBLE_EQ(child->x, xc[i]);
        EXPECT_DOUBLE_EQ(child->y, yc[i]);
        EXPECT_DOUBLE_EQ(child->dx, 4.0);
        EXPECT_DOUBLE_EQ(child->dy, 2.0);
    }
}

TEST(Box, SplitRefusedForLeafSmallerThanEpsilon)
{
    Box box(0.0, 0.0, 1.0, 1.0);
    EXPECT_FALSE(box.split(2.0));
    EXPECT_TRUE(box.isLeaf);
}

TEST(Box, RotateCornersAndCenter)
{
    Box box(10.0, 20.0, 4.0, 2.0);
    double ldx, ldy, lux, luy, rdx, rdy, rux, ruy;
    double angle = 0.5 * std::acos(-1.0);
    box.rotate(0.0, 0.0, angle, ldx, ldy, lux, luy, rdx, rdy, rux, ruy);

    // a quarter turn counterclockwise maps (x, y) to (-y, x)
    EXPECT_NEAR(ldx, -19.0, 1e-12);
    EXPECT_NEAR(ldy, 8.0, 1e-12);
    EXPECT_NEAR(rux, -21.0, 1e-12);
    EXPECT_NEAR(ruy, 12.0, 1e-12);

    double cx, cy;
    box.rotateCtr(0.0, 0.0, angle, cx, cy);
    EXPECT_NEAR(cx, -20.0, 1e-12);
    EXPECT_NEAR(cy, 10.0, 1e-12);
}
