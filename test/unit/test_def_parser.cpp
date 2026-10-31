// Definition file parsing
#include <cmath>
#include <list>
#include <sstream>
#include <string>

#include <gtest/gtest.h>

#include "def_parser.h"
#include "def_struct.h"

using namespace cusg;

namespace
{

const char * MFG_BLOCK = "BEGIN MODFLOW_GRID mfg\n"
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

    EXPECT_DOUBLE_EQ(data->rotation_angle, 30.0 * std::acos(-1.0) / 180);
    EXPECT_EQ(data->rotation_degrees, 30.0);
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
    std::istringstream dfn("BEGIN QUADTREE_BUILDER qtb\n"
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

TEST(DefinitionParser, SurfaceRecordsAfterGapAreRead)
{
    // bottom and top records for layer 3 only, and an out-of-range layer number
    std::istringstream dfn("BEGIN QUADTREE_BUILDER qtb\n"
                           "  MODFLOW_GRID = mfg\n"
                           "  TOP LAYER 3 = REPLICATE mfg\n"
                           "  BOTTOM LAYER 3 = INTERPOLATE mfg\n"
                           "  ACTIVE_DOMAIN LAYER 99999999 = typo\n"
                           "  GRID_DEFINITION_FILE = qtg.dfn\n"
                           "END QUADTREE_BUILDER\n"
                           "BEGIN QUADTREE qtg\n"
                           "  MODFLOW_GRID = mfg\n"
                           "  STRUCTURE_FILE = OPEN/CLOSE qtg.tsf\n"
                           "  TOP LAYER 2 = CONSTANT 0.\n"
                           "  BOTTOM LAYER 2 = CONSTANT -10.\n"
                           "END QUADTREE\n");
    DefinitionParser parser;
    ASSERT_TRUE(parser.parse(dfn));

    quadtree_builder_raw_data * qtb =
        dynamic_cast<quadtree_builder_raw_data *>(parser.getBlockData("qtb"));
    ASSERT_NE(qtb, nullptr);
    ASSERT_EQ(qtb->top_str.count(2), 1u);
    ASSERT_EQ(qtb->bot_str.count(2), 1u);
    EXPECT_EQ(qtb->bot_str[2].front(), "INTERPOLATE");
    EXPECT_EQ(qtb->bot_str.count(0), 0u);
    EXPECT_EQ(qtb->active_domain_str.size(), 1u);

    quadtree_grid_raw_data * qtg =
        dynamic_cast<quadtree_grid_raw_data *>(parser.getBlockData("qtg"));
    ASSERT_NE(qtg, nullptr);
    EXPECT_EQ(qtg->top_str.count(0), 0u);
    ASSERT_EQ(qtg->top_str.count(1), 1u);
    ASSERT_EQ(qtg->bot_str.count(1), 1u);
}

TEST(DefinitionParser, RefinementAliasesAndTopOnlyAreRead)
{
    // both spellings of the refinement records, refine_by_* taking precedence,
    // and a top record with no bottom record for the same layer
    std::istringstream dfn("BEGIN QUADTREE_BUILDER qtb\n"
                           "  MODFLOW_GRID = mfg\n"
                           "  REFINE_BY_FEATURES LAYER 1 = wells\n"
                           "  REFINEMENT_FEATURES LAYER 1 = rivers\n"
                           "  REFINE_BY_ARRAY LAYER 2 = CONSTANT 1\n"
                           "  REFINEMENT_ARRAY LAYER 3 = CONSTANT 2\n"
                           "  TOP LAYER 2 = REPLICATE mfg\n"
                           "  GRID_DEFINITION_FILE = qtg.dfn\n"
                           "END QUADTREE_BUILDER\n");
    DefinitionParser parser;
    ASSERT_TRUE(parser.parse(dfn));

    quadtree_builder_raw_data * data =
        dynamic_cast<quadtree_builder_raw_data *>(parser.getBlockData("qtb"));
    ASSERT_NE(data, nullptr);

    ASSERT_EQ(data->refine_by_features_str.count(0), 1u);
    EXPECT_EQ(data->refine_by_features_str[0].front(), "wells");

    EXPECT_EQ(data->refine_by_array_str.count(0), 0u);
    ASSERT_EQ(data->refine_by_array_str.count(1), 1u);
    EXPECT_EQ(data->refine_by_array_str[1].back(), "1");
    ASSERT_EQ(data->refine_by_array_str.count(2), 1u);
    EXPECT_EQ(data->refine_by_array_str[2].back(), "2");

    ASSERT_EQ(data->top_str.count(1), 1u);
    EXPECT_EQ(data->bot_str.count(1), 0u);
}

namespace
{

quadtree_builder_raw_data * parse_builder(DefinitionParser & parser, const std::string & records)
{
    std::istringstream dfn("BEGIN QUADTREE_BUILDER qtb\n"
                           "  MODFLOW_GRID = mfg\n"
                           "  GRID_DEFINITION_FILE = qtg.dfn\n" +
                           records + "END QUADTREE_BUILDER\n");
    EXPECT_TRUE(parser.parse(dfn));
    return dynamic_cast<quadtree_builder_raw_data *>(parser.getBlockData("qtb"));
}

refinement_raw_data * parse_refinement(DefinitionParser & parser, const std::string & records)
{
    std::istringstream dfn("BEGIN REFINEMENT_FEATURES rf\n"
                           "  SHAPEFILE = wells\n"
                           "  FEATURE_TYPE = point\n" +
                           records + "END REFINEMENT_FEATURES\n");
    EXPECT_TRUE(parser.parse(dfn));
    return dynamic_cast<refinement_raw_data *>(parser.getBlockData("rf"));
}

} // namespace

TEST(DefinitionParser, SmoothingLevelSpellings)
{
    // the report's spelling, the code's spelling, and both (#52)
    const char * spellings[][2] = {{"HORIZONTAL_SMOOTHING_LEVEL", "VERTICAL_SMOOTHING_LEVEL"},
                                   {"SMOOTHING_LEVEL_HORIZONTAL", "SMOOTHING_LEVEL_VERTICAL"}};
    for (const auto & s : spellings)
    {
        SCOPED_TRACE(s[0]);
        DefinitionParser parser;
        quadtree_builder_raw_data * data =
            parse_builder(parser, std::string("  ") + s[0] + " = 3\n  " + s[1] + " = 2\n");
        ASSERT_NE(data, nullptr);
        EXPECT_EQ(data->smoothing_level_horizontal, 3);
        EXPECT_EQ(data->smoothing_level_vertical, 2);
    }

    DefinitionParser parser;
    quadtree_builder_raw_data * data = parse_builder(parser, "  SMOOTHING_LEVEL_HORIZONTAL = 2\n"
                                                             "  HORIZONTAL_SMOOTHING_LEVEL = 4\n");
    ASSERT_NE(data, nullptr);
    EXPECT_EQ(data->smoothing_level_horizontal, 4);
}

TEST(DefinitionParser, RefineLevelByAttributeSpellings)
{
    // either spelling, and the attribute overrides REFINEMENT_LEVEL (#52)
    for (const char * name : {"REFINE_LEVEL_BY_ATTRIBUTE", "REFINEMENT_LEVEL_BY_ATTRIBUTE"})
    {
        SCOPED_TRACE(name);
        DefinitionParser parser;
        refinement_raw_data * data = parse_refinement(
            parser, std::string("  REFINEMENT_LEVEL = 2\n  ") + name + " = LEVEL\n");
        ASSERT_NE(data, nullptr);
        EXPECT_EQ(data->refinement_level_by_attribute, "LEVEL");
    }

    DefinitionParser parser;
    refinement_raw_data * data = parse_refinement(parser, "  REFINEMENT_LEVEL = 2\n");
    ASSERT_NE(data, nullptr);
    EXPECT_EQ(data->refinement_level, 2);
    EXPECT_TRUE(data->refinement_level_by_attribute.empty());
}

TEST(DefinitionParser, UnrecognizedRecordIsReported)
{
    // a misspelled record is reported, and accepted records are not (#52)
    DefinitionParser parser;
    parse_builder(parser, "  SMOOTHING = full\n"
                          "  TOP LAYER 1 = REPLICATE mfg\n"
                          "  REFINEMENT_FEATURES LAYER 2 = wells\n"
                          "  HORIZONTAL_SMOOTHING_LEVLE = 2\n");
    const std::list<std::string> & warnings = parser.getWarnings();
    ASSERT_EQ(warnings.size(), 1u);
    EXPECT_NE(warnings.front().find("HORIZONTAL_SMOOTHING_LEVLE"), std::string::npos)
        << warnings.front();
    EXPECT_NE(warnings.front().find("QUADTREE_BUILDER block qtb"), std::string::npos)
        << warnings.front();
}
