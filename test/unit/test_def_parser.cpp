// Definition file parsing
#include <cmath>
#include <functional>
#include <list>
#include <memory>
#include <sstream>
#include <string>
#include <utility>
#include <vector>

#include <gtest/gtest.h>

#include "GridgenError.h"
#include "Warnings.h"
#include "def_parser.h"
#include "def_struct.h"

using namespace cusg;

namespace
{

const char * mfg_block = "BEGIN MODFLOW_GRID mfg\n"
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

// mfg_block with one record replaced
std::string mfg_block_with(const std::string & from, const std::string & to)
{
    std::string block = mfg_block;
    block.replace(block.find(from), from.size(), to);
    return block;
}

// the message of the GridgenError that action throws, or "" if it throws none
std::string error_message(const std::function<void()> & action)
{
    try
    {
        action();
    }
    catch (const GridgenError & e)
    {
        return e.what();
    }
    return "";
}

// the message of the error from initializing block name of the definition file
std::string initialize_error(const std::string & definition, const std::string & name)
{
    std::istringstream dfn(definition);
    DefinitionParser parser;
    parser.parse(dfn);
    DefinitionBlock * data = parser.get_block_data(name);
    if (data == nullptr) return "no block " + name;
    return error_message([&] { data->initialize(parser.get_block_data()); });
}

} // namespace

TEST(DefinitionParser, ModflowGridBlock)
{
    std::istringstream dfn(mfg_block);
    DefinitionParser parser;
    ASSERT_TRUE(parser.parse(dfn));

    ModflowGridBlock * data = dynamic_cast<ModflowGridBlock *>(parser.get_block_data("mfg"));
    ASSERT_NE(data, nullptr);
    ASSERT_TRUE(data->initialize(parser.get_block_data()));

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

TEST(DefinitionParser, MissingRequiredRecordIsAnError)
{
    std::istringstream dfn(mfg_block_with("  NROW = 3\n", ""));
    DefinitionParser parser;
    EXPECT_EQ(error_message([&] { parser.parse(dfn); }),
              "MODFLOW_GRID block mfg has no NROW record");

    std::istringstream bottom(mfg_block_with("  BOTTOM LAYER 2 = CONSTANT -10.\n", ""));
    EXPECT_EQ(error_message([&] { DefinitionParser().parse(bottom); }),
              "MODFLOW_GRID block mfg has no BOTTOM LAYER 2 record");
}

TEST(DefinitionParser, ReferenceToAnUndefinedBlockIsAnError)
{
    EXPECT_EQ(initialize_error("BEGIN GRID_TO_SHAPEFILE shp\n"
                               "  GRID = nogrid\n"
                               "  SHAPEFILE = out\n"
                               "  FEATURE_TYPE = polygon\n"
                               "END GRID_TO_SHAPEFILE\n",
                               "shp"),
              "record GRID in GRID_TO_SHAPEFILE block shp names block nogrid, which is not "
              "defined");

    EXPECT_EQ(initialize_error("BEGIN QUADTREE_BUILDER qtb\n"
                               "  MODFLOW_GRID = nogrid\n"
                               "  GRID_DEFINITION_FILE = qtg.dfn\n"
                               "END QUADTREE_BUILDER\n",
                               "qtb"),
              "record MODFLOW_GRID in QUADTREE_BUILDER block qtb names block nogrid, which is "
              "not defined");

    EXPECT_EQ(initialize_error(std::string(mfg_block) + "BEGIN QUADTREE_BUILDER qtb\n"
                                                        "  MODFLOW_GRID = mfg\n"
                                                        "  REFINEMENT_FEATURES LAYER 2 = typo\n"
                                                        "  GRID_DEFINITION_FILE = qtg.dfn\n"
                                                        "END QUADTREE_BUILDER\n",
                               "qtb"),
              "record REFINEMENT_FEATURES LAYER 2 in QUADTREE_BUILDER block qtb names block "
              "typo, which is not defined");
}

TEST(DefinitionParser, ReferenceToABlockOfTheWrongTypeIsAnError)
{
    const std::string domain = "BEGIN ACTIVE_DOMAIN domain\n"
                               "  SHAPEFILE = domain\n"
                               "  FEATURE_TYPE = polygon\n"
                               "END ACTIVE_DOMAIN\n";

    EXPECT_EQ(initialize_error(domain + "BEGIN GRID_TO_SHAPEFILE shp\n"
                                        "  GRID = domain\n"
                                        "  SHAPEFILE = out\n"
                                        "  FEATURE_TYPE = polygon\n"
                                        "END GRID_TO_SHAPEFILE\n",
                               "shp"),
              "record GRID in GRID_TO_SHAPEFILE block shp names ACTIVE_DOMAIN block domain, "
              "which is not a MODFLOW_GRID, QUADTREE, or QUADTREE_BUILDER block");

    EXPECT_EQ(initialize_error(std::string(mfg_block) + domain +
                                   "BEGIN QUADTREE_BUILDER qtb\n"
                                   "  MODFLOW_GRID = mfg\n"
                                   "  REFINEMENT_FEATURES LAYER 1 = domain\n"
                                   "  GRID_DEFINITION_FILE = qtg.dfn\n"
                                   "END QUADTREE_BUILDER\n",
                               "qtb"),
              "record REFINEMENT_FEATURES LAYER 1 in QUADTREE_BUILDER block qtb names "
              "ACTIVE_DOMAIN block domain, which is not a REFINEMENT_FEATURES block");
}

TEST(DefinitionParser, ArrayErrors)
{
    // DELC is read after DELR, which an error must not leave behind
    for (const char * delc : {"OPEN/CLOSE gridgen_unit_missing.dat", "BOGUS 25."})
    {
        SCOPED_TRACE(delc);
        std::istringstream dfn(mfg_block_with("CONSTANT 25.", delc));
        DefinitionParser parser;
        ASSERT_TRUE(parser.parse(dfn));
        ModflowGridBlock * data = dynamic_cast<ModflowGridBlock *>(parser.get_block_data("mfg"));
        ASSERT_NE(data, nullptr);
        EXPECT_THROW(data->initialize(parser.get_block_data()), GridgenError);
        EXPECT_EQ(data->delr, nullptr);
        EXPECT_TRUE(data->bot.empty());
    }
}

TEST(DefinitionParser, UnknownBlockIsNull)
{
    std::istringstream dfn(mfg_block);
    DefinitionParser parser;
    ASSERT_TRUE(parser.parse(dfn));
    EXPECT_EQ(parser.get_block_data("missing"), nullptr);
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

    QuadtreeBuilderBlock * data =
        dynamic_cast<QuadtreeBuilderBlock *>(parser.get_block_data("qtb"));
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

    QuadtreeBuilderBlock * qtb = dynamic_cast<QuadtreeBuilderBlock *>(parser.get_block_data("qtb"));
    ASSERT_NE(qtb, nullptr);
    ASSERT_EQ(qtb->top_str.count(2), 1u);
    ASSERT_EQ(qtb->bot_str.count(2), 1u);
    EXPECT_EQ(qtb->bot_str[2].front(), "INTERPOLATE");
    EXPECT_EQ(qtb->bot_str.count(0), 0u);
    EXPECT_EQ(qtb->active_domain_str.size(), 1u);

    QuadtreeBlock * qtg = dynamic_cast<QuadtreeBlock *>(parser.get_block_data("qtg"));
    ASSERT_NE(qtg, nullptr);
    EXPECT_EQ(qtg->top_str.count(0), 0u);
    ASSERT_EQ(qtg->top_str.count(1), 1u);
    ASSERT_EQ(qtg->bot_str.count(1), 1u);
}

TEST(DefinitionParser, RefinementAliasesAndTopOnlyAreRead)
{
    // both spellings of the refinement features record, refine_by_features taking precedence,
    // and a top record with no bottom record for the same layer
    std::istringstream dfn("BEGIN QUADTREE_BUILDER qtb\n"
                           "  MODFLOW_GRID = mfg\n"
                           "  REFINE_BY_FEATURES LAYER 1 = wells\n"
                           "  REFINEMENT_FEATURES LAYER 1 = rivers\n"
                           "  TOP LAYER 2 = REPLICATE mfg\n"
                           "  GRID_DEFINITION_FILE = qtg.dfn\n"
                           "END QUADTREE_BUILDER\n");
    DefinitionParser parser;
    ASSERT_TRUE(parser.parse(dfn));

    QuadtreeBuilderBlock * data =
        dynamic_cast<QuadtreeBuilderBlock *>(parser.get_block_data("qtb"));
    ASSERT_NE(data, nullptr);

    ASSERT_EQ(data->refine_by_features_str.count(0), 1u);
    EXPECT_EQ(data->refine_by_features_str[0].front(), "wells");

    ASSERT_EQ(data->top_str.count(1), 1u);
    EXPECT_EQ(data->bot_str.count(1), 0u);
}

namespace
{

QuadtreeBuilderBlock * parse_builder(DefinitionParser & parser, const std::string & records)
{
    std::istringstream dfn("BEGIN QUADTREE_BUILDER qtb\n"
                           "  MODFLOW_GRID = mfg\n"
                           "  GRID_DEFINITION_FILE = qtg.dfn\n" +
                           records + "END QUADTREE_BUILDER\n");
    EXPECT_TRUE(parser.parse(dfn));
    return dynamic_cast<QuadtreeBuilderBlock *>(parser.get_block_data("qtb"));
}

RefinementFeaturesBlock * parse_refinement(DefinitionParser & parser, const std::string & records)
{
    std::istringstream dfn("BEGIN REFINEMENT_FEATURES rf\n"
                           "  SHAPEFILE = wells\n"
                           "  FEATURE_TYPE = point\n" +
                           records + "END REFINEMENT_FEATURES\n");
    EXPECT_TRUE(parser.parse(dfn));
    return dynamic_cast<RefinementFeaturesBlock *>(parser.get_block_data("rf"));
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
        QuadtreeBuilderBlock * data =
            parse_builder(parser, std::string("  ") + s[0] + " = 3\n  " + s[1] + " = 2\n");
        ASSERT_NE(data, nullptr);
        EXPECT_EQ(data->smoothing_level_horizontal, 3);
        EXPECT_EQ(data->smoothing_level_vertical, 2);
    }

    DefinitionParser parser;
    QuadtreeBuilderBlock * data = parse_builder(parser, "  SMOOTHING_LEVEL_HORIZONTAL = 2\n"
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
        RefinementFeaturesBlock * data = parse_refinement(
            parser, std::string("  REFINEMENT_LEVEL = 2\n  ") + name + " = LEVEL\n");
        ASSERT_NE(data, nullptr);
        EXPECT_EQ(data->refinement_level_by_attribute, "LEVEL");
    }

    DefinitionParser parser;
    RefinementFeaturesBlock * data = parse_refinement(parser, "  REFINEMENT_LEVEL = 2\n");
    ASSERT_NE(data, nullptr);
    EXPECT_EQ(data->refinement_level, 2);
    EXPECT_TRUE(data->refinement_level_by_attribute.empty());
}

TEST(DefinitionParser, RefinementArrayIsReportedAsIgnored)
{
    // REFINEMENT_ARRAY was read but never used, without a message (#93)
    const char * names[] = {"REFINEMENT_ARRAY", "REFINE_BY_ARRAY"};
    for (const char * name : names)
    {
        DefinitionParser parser;
        parse_builder(parser, std::string("  ") + name + " LAYER 1 = OPEN/CLOSE refine.dat\n");
        const std::list<std::string> & warnings = parser.get_warnings();
        ASSERT_EQ(warnings.size(), 1u) << name;
        EXPECT_NE(warnings.front().find(std::string(name) + " LAYER 1"), std::string::npos)
            << warnings.front();
        EXPECT_NE(warnings.front().find("is ignored"), std::string::npos) << warnings.front();
    }
}

TEST(DefinitionParser, RefinementArrayFileIsNotRead)
{
    // an ignored refinement array whose file does not exist does not stop initialization (#93)
    std::istringstream dfn(std::string(mfg_block) +
                           "BEGIN QUADTREE_BUILDER qtb\n"
                           "  MODFLOW_GRID = mfg\n"
                           "  GRID_DEFINITION_FILE = qtg.dfn\n"
                           "  REFINEMENT_ARRAY LAYER 1 = OPEN/CLOSE no_such_refine.dat\n"
                           "END QUADTREE_BUILDER\n");
    DefinitionParser parser;
    ASSERT_TRUE(parser.parse(dfn));
    QuadtreeBuilderBlock * data =
        dynamic_cast<QuadtreeBuilderBlock *>(parser.get_block_data("qtb"));
    ASSERT_NE(data, nullptr);
    EXPECT_TRUE(data->initialize(parser.get_block_data()));
}

TEST(DefinitionParser, UnrecognizedRecordIsReported)
{
    // a misspelled record is reported, and accepted records are not (#52)
    DefinitionParser parser;
    parse_builder(parser, "  SMOOTHING = full\n"
                          "  TOP LAYER 1 = REPLICATE mfg\n"
                          "  REFINEMENT_FEATURES LAYER 2 = wells\n"
                          "  HORIZONTAL_SMOOTHING_LEVLE = 2\n");
    const std::list<std::string> & warnings = parser.get_warnings();
    ASSERT_EQ(warnings.size(), 1u);
    EXPECT_NE(warnings.front().find("HORIZONTAL_SMOOTHING_LEVLE"), std::string::npos)
        << warnings.front();
    EXPECT_NE(warnings.front().find("QUADTREE_BUILDER block qtb"), std::string::npos)
        << warnings.front();
}

TEST(DefinitionParser, IntersectionVtkfileIsReportedAsIgnored)
{
    // OUTPUT_VTKFILE was read but no VTK file was written, without a message (#102)
    std::istringstream dfn("BEGIN GRID_INTERSECTION gi\n"
                           "  GRID = mfg\n"
                           "  LAYER = 1\n"
                           "  SHAPEFILE = wells\n"
                           "  FEATURE_TYPE = point\n"
                           "  OUTPUT_FILE = wells.ifo\n"
                           "  OUTPUT_VTKFILE = wells\n"
                           "END GRID_INTERSECTION\n");
    DefinitionParser parser;
    ASSERT_TRUE(parser.parse(dfn));
    const std::list<std::string> & warnings = parser.get_warnings();
    ASSERT_EQ(warnings.size(), 1u);
    EXPECT_NE(warnings.front().find("OUTPUT_VTKFILE"), std::string::npos) << warnings.front();
    EXPECT_NE(warnings.front().find("is ignored"), std::string::npos) << warnings.front();
}

namespace
{

// GRID_TO_SHAPEFILE and GRID_TO_VTKFILE blocks of a grid, with the records
// that only a MODFLOW_GRID uses
std::string export_blocks(const std::string & grid)
{
    return "BEGIN QUADTREE_BUILDER qtb\n"
           "  MODFLOW_GRID = mfg\n"
           "  GRID_DEFINITION_FILE = qtg.dfn\n"
           "END QUADTREE_BUILDER\n"
           "BEGIN GRID_TO_SHAPEFILE shp\n"
           "  GRID = " +
           grid +
           "\n"
           "  SHAPEFILE = grid\n"
           "  FEATURE_TYPE = polygon\n"
           "  NODE_NUMBERING = zero-based\n"
           "  WITHOUT_INACTIVE = false\n"
           "END GRID_TO_SHAPEFILE\n"
           "BEGIN GRID_TO_VTKFILE vtk\n"
           "  GRID = " +
           grid +
           "\n"
           "  VTKFILE = grid\n"
           "  NODE_NUMBERING = zero-based\n"
           "  WITHOUT_INACTIVE = false\n"
           "END GRID_TO_VTKFILE\n";
}

} // namespace

TEST(DefinitionParser, QuadtreeExportRecordsAreReportedAsIgnored)
{
    // a quadtree grid is written with its own numbering and without inactive cells (#102)
    std::istringstream dfn(std::string(mfg_block) + export_blocks("qtb"));
    DefinitionParser parser;
    ASSERT_TRUE(parser.parse(dfn));
    const char * blocks[] = {"shp", "vtk"};
    for (const char * name : blocks)
    {
        DefinitionBlock * data = parser.get_block_data(name);
        ASSERT_NE(data, nullptr) << name;
        ASSERT_TRUE(data->initialize(parser.get_block_data())) << name;
        const std::list<std::string> & warnings = data->get_warnings();
        ASSERT_EQ(warnings.size(), 2u) << name;
        EXPECT_NE(warnings.front().find("NODE_NUMBERING"), std::string::npos) << warnings.front();
        EXPECT_NE(warnings.back().find("WITHOUT_INACTIVE"), std::string::npos) << warnings.back();
    }
}

TEST(DefinitionParser, ModflowGridExportRecordsAreNotReported)
{
    // the same records are used for a MODFLOW_GRID (#102)
    std::istringstream dfn(std::string(mfg_block) + export_blocks("mfg"));
    DefinitionParser parser;
    ASSERT_TRUE(parser.parse(dfn));
    const char * blocks[] = {"shp", "vtk"};
    for (const char * name : blocks)
    {
        DefinitionBlock * data = parser.get_block_data(name);
        ASSERT_NE(data, nullptr) << name;
        ASSERT_TRUE(data->initialize(parser.get_block_data())) << name;
        EXPECT_TRUE(data->get_warnings().empty()) << name;
    }
}

TEST(DefinitionParser, QuadtreeExportWithoutTheRecordsIsNotReported)
{
    // only records that are given are reported (#102)
    std::istringstream dfn(std::string(mfg_block) + "BEGIN QUADTREE_BUILDER qtb\n"
                                                    "  MODFLOW_GRID = mfg\n"
                                                    "  GRID_DEFINITION_FILE = qtg.dfn\n"
                                                    "END QUADTREE_BUILDER\n"
                                                    "BEGIN GRID_TO_SHAPEFILE shp\n"
                                                    "  GRID = qtb\n"
                                                    "  SHAPEFILE = grid\n"
                                                    "  FEATURE_TYPE = point\n"
                                                    "END GRID_TO_SHAPEFILE\n");
    DefinitionParser parser;
    ASSERT_TRUE(parser.parse(dfn));
    DefinitionBlock * data = parser.get_block_data("shp");
    ASSERT_NE(data, nullptr);
    ASSERT_TRUE(data->initialize(parser.get_block_data()));
    EXPECT_TRUE(data->get_warnings().empty());
}

// a block of type T, named after its type, with no records
template <typename T> std::unique_ptr<DefinitionBlock> empty_block(const std::string & type)
{
    DefinitionBlock block;
    block.type = type;
    block.name = type;
    return std::make_unique<T>(block);
}

TEST(DefinitionParser, RunningABlockThatIsNotAnActionWarns)
{
    // running an object block writes nothing (#151); running an action does
    std::vector<std::pair<std::unique_ptr<DefinitionBlock>, bool>> blocks;
    blocks.emplace_back(empty_block<ModflowGridBlock>("modflow_grid"), false);
    blocks.emplace_back(empty_block<QuadtreeBlock>("quadtree"), false);
    blocks.emplace_back(empty_block<RefinementFeaturesBlock>("refinement_features"), false);
    blocks.emplace_back(empty_block<ActiveDomainBlock>("active_domain"), false);
    blocks.emplace_back(empty_block<QuadtreeBuilderBlock>("quadtree_builder"), true);
    blocks.emplace_back(empty_block<GridToShapefileBlock>("grid_to_shapefile"), true);
    blocks.emplace_back(empty_block<GridToUsgdataBlock>("grid_to_usgdata"), true);
    blocks.emplace_back(empty_block<GridToVtkfileBlock>("grid_to_vtkfile"), true);
    blocks.emplace_back(empty_block<GridIntersectionBlock>("grid_intersection"), true);
    blocks.emplace_back(empty_block<GridToGsfBlock>("grid_to_gsf"), true);
    for (auto & [block, is_action] : blocks)
    {
        clear_warnings();
        block->warn_if_not_action();
        if (is_action)
        {
            EXPECT_TRUE(get_warnings().empty()) << block->label();
            continue;
        }
        ASSERT_EQ(get_warnings().size(), 1u) << block->label();
        EXPECT_EQ(get_warnings()[0].rfind(block->label() + " is not an action", 0), 0u)
            << get_warnings()[0];
    }
    clear_warnings();
}

TEST(DefinitionParser, GridToGsfBlock)
{
    std::istringstream dfn(std::string(mfg_block) + "BEGIN GRID_TO_GSF gsf\n"
                                                    "  GRID = mfg\n"
                                                    "  GSF_FILE = output\\mfg.gsf\n"
                                                    "  HANGING_VERTICES = TRUE\n"
                                                    "  NODE_DATA_TABLE_FILE = mfg.ndt\n"
                                                    "END GRID_TO_GSF\n"
                                                    "BEGIN GRID_TO_GSF defaults\n"
                                                    "  GRID = mfg\n"
                                                    "  GSF_FILE = mfg.gsf\n"
                                                    "END GRID_TO_GSF\n");
    DefinitionParser parser;
    ASSERT_TRUE(parser.parse(dfn));

    auto * gsf = dynamic_cast<GridToGsfBlock *>(parser.get_block_data("gsf"));
    ASSERT_NE(gsf, nullptr);
    ASSERT_TRUE(gsf->initialize(parser.get_block_data()));
    EXPECT_NE(gsf->grid, nullptr);
    EXPECT_EQ(gsf->gsf_file, "output/mfg.gsf");
    EXPECT_TRUE(gsf->hanging_vertices);
    EXPECT_EQ(gsf->node_data_table_file, "mfg.ndt");
    // a MODFLOW grid has no hanging vertices
    ASSERT_EQ(gsf->get_warnings().size(), 1u);
    EXPECT_NE(gsf->get_warnings().front().find("HANGING_VERTICES"), std::string::npos);

    auto * defaults = dynamic_cast<GridToGsfBlock *>(parser.get_block_data("defaults"));
    ASSERT_NE(defaults, nullptr);
    ASSERT_TRUE(defaults->initialize(parser.get_block_data()));
    EXPECT_FALSE(defaults->hanging_vertices);
    EXPECT_TRUE(defaults->node_data_table_file.empty());
    EXPECT_TRUE(defaults->get_warnings().empty());
}

TEST(DefinitionParser, GridToGsfHangingVerticesOfQuadtreeIsNotReported)
{
    std::istringstream dfn(std::string(mfg_block) + "BEGIN QUADTREE_BUILDER qtb\n"
                                                    "  MODFLOW_GRID = mfg\n"
                                                    "  GRID_DEFINITION_FILE = qtg.dfn\n"
                                                    "END QUADTREE_BUILDER\n"
                                                    "BEGIN GRID_TO_GSF gsf\n"
                                                    "  GRID = qtb\n"
                                                    "  GSF_FILE = qtg.gsf\n"
                                                    "  HANGING_VERTICES = True\n"
                                                    "END GRID_TO_GSF\n");
    DefinitionParser parser;
    ASSERT_TRUE(parser.parse(dfn));
    DefinitionBlock * data = parser.get_block_data("gsf");
    ASSERT_NE(data, nullptr);
    ASSERT_TRUE(data->initialize(parser.get_block_data()));
    EXPECT_TRUE(data->get_warnings().empty());
}

TEST(DefinitionParser, GridToGsfWithoutFileIsAnError)
{
    std::istringstream dfn(std::string(mfg_block) + "BEGIN GRID_TO_GSF gsf\n"
                                                    "  GRID = mfg\n"
                                                    "END GRID_TO_GSF\n");
    EXPECT_EQ(error_message([&] { DefinitionParser().parse(dfn); }),
              "GRID_TO_GSF block gsf has no GSF_FILE record");
}
