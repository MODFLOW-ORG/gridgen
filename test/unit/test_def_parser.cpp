// Definition file parsing
#include <sstream>
#include <string>

#include <gtest/gtest.h>

#include "def_parser.h"
#include "def_struct.h"

using namespace cusg;

namespace
{

const char * MFG_BLOCK =
    "BEGIN MODFLOW_GRID mfg\n"
    "  ROTATION_ANGLE = 30\n"
    "  X_OFFSET = 1000.5\n"
    "  Y_OFFSET = 2000.25\n"
    "  LENGTH_UNIT = meters\n"
    "  NLAY = 2\n"
    "  NROW = 3\n"
    "  NCOL = 4\n"
    "  DELR = CONSTANT 50.\n"
    "  DELC = CONSTANT 25.\n"
    "  TOP = CONSTANT 10.\n"
    "  BOTTOM LAYER 1 = CONSTANT 0.\n"
    "  BOTTOM LAYER 2 = CONSTANT -10.\n"
    "END MODFLOW_GRID\n";

} // namespace

TEST(DefinitionParser, ModflowGridBlock)
{
    std::istringstream dfn(MFG_BLOCK);
    DefinitionParser parser;
    ASSERT_TRUE(parser.parse(dfn));

    modflow_grid_raw_data * data =
        dynamic_cast<modflow_grid_raw_data *>(parser.getBlockData("mfg"));
    ASSERT_NE(data, nullptr);
    ASSERT_TRUE(data->initialize(parser.getBlockData()));

    EXPECT_EQ(data->nlay, 2);
    EXPECT_EQ(data->nrow, 3);
    EXPECT_EQ(data->ncol, 4);
    EXPECT_DOUBLE_EQ(data->x_offset, 1000.5);
    EXPECT_DOUBLE_EQ(data->y_offset, 2000.25);
    for (int j = 0; j < data->ncol; j++) EXPECT_DOUBLE_EQ(data->delr[j], 50.0);
    for (int i = 0; i < data->nrow; i++) EXPECT_DOUBLE_EQ(data->delc[i], 25.0);
    for (int i = 0; i < 12; i++)
    {
        EXPECT_DOUBLE_EQ(data->bot[0][i], 10.0);
        EXPECT_DOUBLE_EQ(data->bot[1][i], 0.0);
        EXPECT_DOUBLE_EQ(data->bot[2][i], -10.0);
    }

    // degrees are converted with pi = 3.1415926
    EXPECT_DOUBLE_EQ(data->rotation_angle, 30.0 * 2 * 3.1415926 / 360);
}

TEST(DefinitionParser, UnknownBlockIsNull)
{
    std::istringstream dfn(MFG_BLOCK);
    DefinitionParser parser;
    ASSERT_TRUE(parser.parse(dfn));
    EXPECT_EQ(parser.getBlockData("missing"), nullptr);
}

TEST(DefinitionParser, LayerRecordsReadIndependently)
{
    // an active domain on layer 2 only and refinement on layer 3 only, with
    // no records for layer 1 (#49)
    std::istringstream dfn(
        "BEGIN QUADTREE_BUILDER qtb\n"
        "  MODFLOW_GRID = mfg\n"
        "  ACTIVE_DOMAIN LAYER 2 = domain\n"
        "  REFINEMENT_FEATURES LAYER 3 = wells\n"
        "  GRID_DEFINITION_FILE = qtg.dfn\n"
        "END QUADTREE_BUILDER\n");
    DefinitionParser parser;
    ASSERT_TRUE(parser.parse(dfn));

    quadtree_builder_raw_data * data =
        dynamic_cast<quadtree_builder_raw_data *>(parser.getBlockData("qtb"));
    ASSERT_NE(data, nullptr);

    EXPECT_EQ(data->active_domain_str.count(0), 0u);
    ASSERT_EQ(data->active_domain_str.count(1), 1u);
    EXPECT_EQ(data->active_domain_str[1].front(), "domain");
    EXPECT_EQ(data->active_domain_str.count(2), 0u);

    ASSERT_EQ(data->refine_by_features_str.count(2), 1u);
    EXPECT_EQ(data->refine_by_features_str[2].front(), "wells");
}
