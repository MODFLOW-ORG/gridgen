// File names in the grid definition files that gridgen writes
#include <filesystem>
#include <fstream>
#include <memory>
#include <sstream>
#include <string>
#include <vector>

#include <gtest/gtest.h>

#include "def_exporter.h"
#include "def_parser.h"
#include "Warnings.h"
#include "def_struct.h"
#include "grid_helpers.h"

using namespace cusg;
using gridgen_test::make_mfgrid;
using gridgen_test::make_qtree;

namespace
{

// a MODFLOW_GRID block read from a definition file in folder "in"
const char * mfg_block = "BEGIN MODFLOW_GRID mfg\n"
                         "  ROTATION_ANGLE = 0\n"
                         "  X_OFFSET = 0\n"
                         "  Y_OFFSET = 0\n"
                         "  NLAY = 1\n"
                         "  NROW = 2\n"
                         "  NCOL = 2\n"
                         "  DELR = CONSTANT 50.\n"
                         "  DELC = OPEN/CLOSE arrays\\delc.dat\n"
                         "  TOP = OPEN/CLOSE top.dat\n"
                         "  BOTTOM LAYER 1 = open/close ../bottom.dat\n"
                         "END MODFLOW_GRID\n";

// the MODFLOW_GRID block of mfg_block as written for a definition file in output_folder
std::string written_block(const std::string & output_folder)
{
    std::istringstream dfn(mfg_block);
    DefinitionParser parser;
    EXPECT_TRUE(parser.parse(dfn, "in"));
    ModflowGridBlock * data = dynamic_cast<ModflowGridBlock *>(parser.get_block_data("mfg"));
    EXPECT_NE(data, nullptr);
    if (data == nullptr) return "";
    DefinitionExporter exporter;
    exporter.set_output_folder(output_folder);
    std::ostringstream out;
    exporter.save(out, *data);
    return out.str();
}

} // namespace

TEST(DefinitionExporter, ArraysInTheSameFolderAreWrittenAsGiven)
{
    std::string block = written_block("in");
    EXPECT_NE(block.find("DELR =  CONSTANT 50.\n"), std::string::npos) << block;
    EXPECT_NE(block.find("DELC =  OPEN/CLOSE arrays\\delc.dat\n"), std::string::npos) << block;
    EXPECT_NE(block.find("TOP =  OPEN/CLOSE top.dat\n"), std::string::npos) << block;
    EXPECT_NE(block.find("LAYER 1 =  open/close ../bottom.dat\n"), std::string::npos) << block;
}

TEST(DefinitionExporter, ArraysInAnotherFolderAreRelativeToTheFile)
{
    // written in the folder gridgen is run from
    std::string block = written_block("");
    EXPECT_NE(block.find("DELR =  CONSTANT 50.\n"), std::string::npos) << block;
    EXPECT_NE(block.find("DELC =  OPEN/CLOSE in/arrays/delc.dat\n"), std::string::npos) << block;
    EXPECT_NE(block.find("TOP =  OPEN/CLOSE in/top.dat\n"), std::string::npos) << block;
    EXPECT_NE(block.find("LAYER 1 =  open/close in/../bottom.dat\n"), std::string::npos) << block;

    // written in a folder two levels down
    block = written_block("out/grids");
    EXPECT_NE(block.find("DELC =  OPEN/CLOSE ../../in/arrays/delc.dat\n"), std::string::npos)
        << block;
    EXPECT_NE(block.find("TOP =  OPEN/CLOSE ../../in/top.dat\n"), std::string::npos) << block;
    EXPECT_NE(block.find("LAYER 1 =  open/close ../../bottom.dat\n"), std::string::npos) << block;
}

TEST(DefinitionExporter, QuadtreeFilesAreRelativeToTheFile)
{
    // a grid named for a definition file in a folder below the temporary folder
    std::filesystem::path folder =
        std::filesystem::path(::testing::TempDir()) / "gridgen_unit_grid_definition";
    std::filesystem::create_directories(folder);
    std::unique_ptr<ModflowGrid> mfgrid = make_mfgrid(1, 2, 2, 100.0, 100.0);
    mfgrid->name = "mfg";
    std::unique_ptr<QuadTree3D> qtree = make_qtree(mfgrid.get(), {{50.0, 150.0, 1, 0}});
    qtree->name = (folder / "qtg").generic_string();

    DefinitionExporter exporter;
    exporter.set_output_folder(folder.generic_string());
    std::ostringstream out;
    exporter.save(out, qtree.get());
    std::string block = out.str();
    EXPECT_NE(block.find("BEGIN QUADTREE qtg\n"), std::string::npos) << block;
    EXPECT_NE(block.find("STRUCTURE_File  = OPEN/CLOSE qtg.tsf\n"), std::string::npos) << block;
    EXPECT_NE(block.find("TOP LAYER 1 = OPEN/CLOSE qtg.top1.dat\n"), std::string::npos) << block;
    EXPECT_NE(block.find("BOTTOM LAYER 1 = OPEN/CLOSE qtg.bot1.dat\n"), std::string::npos) << block;
    // the files are written where the name of the grid says
    EXPECT_TRUE(std::filesystem::exists(folder / "qtg.tsf"));
    std::filesystem::remove_all(folder);
}

TEST(DefinitionExporter, GridFilesAreWrittenBesideTheFileInAFolderWithAPeriod)
{
    // a definition file without an extension in a folder with a period in its name (#184)
    std::filesystem::path folder =
        std::filesystem::path(::testing::TempDir()) / "gridgen_unit_out.v2";
    std::filesystem::create_directories(folder);
    std::string definition_file = (folder / "qtg").generic_string();
    std::istringstream dfn("BEGIN MODFLOW_GRID mfg\n"
                           "  ROTATION_ANGLE = 0\n"
                           "  X_OFFSET = 0\n"
                           "  Y_OFFSET = 0\n"
                           "  NLAY = 1\n"
                           "  NROW = 2\n"
                           "  NCOL = 2\n"
                           "  DELR = CONSTANT 50.\n"
                           "  DELC = CONSTANT 50.\n"
                           "  TOP = CONSTANT 10.\n"
                           "  BOTTOM LAYER 1 = CONSTANT 0.\n"
                           "END MODFLOW_GRID\n"
                           "BEGIN QUADTREE_BUILDER qtb\n"
                           "  MODFLOW_GRID = mfg\n"
                           "  GRID_DEFINITION_FILE = " +
                           definition_file + "\nEND QUADTREE_BUILDER\n");
    DefinitionParser parser;
    ASSERT_TRUE(parser.parse(dfn));
    QuadtreeBuilderBlock * data =
        dynamic_cast<QuadtreeBuilderBlock *>(parser.get_block_data("qtb"));
    ASSERT_NE(data, nullptr);
    for (const char * name : {"mfg", "qtb"})
        parser.get_block_data(name)->initialize(parser.get_block_data());
    ASSERT_NE(data->build(), nullptr);

    DefinitionExporter exporter;
    ASSERT_TRUE(exporter.save_grid_definition(definition_file, *data));
    std::stringstream block;
    {
        //closed before the folder is removed, which Windows does not allow for an open file
        std::ifstream written(definition_file);
        block << written.rdbuf();
    }
    EXPECT_NE(block.str().find("BEGIN QUADTREE qtg\n"), std::string::npos) << block.str();
    EXPECT_NE(block.str().find("STRUCTURE_File  = OPEN/CLOSE qtg.tsf\n"), std::string::npos)
        << block.str();
    EXPECT_TRUE(std::filesystem::exists(folder / "qtg.tsf"));
    EXPECT_TRUE(std::filesystem::exists(folder / "qtg.top1.dat"));
    EXPECT_TRUE(std::filesystem::exists(folder / "qtg.bot1.dat"));
    std::filesystem::remove_all(folder);
}

namespace
{

// the warnings of exporting the USG data of qtree, with vertical ghost nodes, to
// a prefix in the temporary folder
std::vector<std::string> usg_data_warnings(QuadTree3D * qtree)
{
    std::istringstream dfn(
        "BEGIN GRID_TO_USGDATA usg\n"
        "  GRID = qtg\n"
        "  USG_DATA_PREFIX = " +
        (std::filesystem::path(::testing::TempDir()) / "gridgen_unit_gnc_edge").generic_string() +
        "\n  GNC_VERTICAL = True\nEND GRID_TO_USGDATA\n");
    DefinitionParser parser;
    EXPECT_TRUE(parser.parse(dfn));
    auto * usg = dynamic_cast<GridToUsgdataBlock *>(parser.get_block_data("usg"));
    EXPECT_NE(usg, nullptr);
    if (usg == nullptr) return {};
    usg->grid = qtree;
    clear_warnings();
    DefinitionExporter exporter;
    exporter.save(*usg);
    std::vector<std::string> warnings = get_warnings();
    clear_warnings();
    return warnings;
}

} // namespace

TEST(DefinitionExporter, GhostNodesBeyondTheActiveCellsAreReported)
{
    // 3 x 3 cells of 100 in 2 layers with the center of layer 1 refined and,
    // in layer 2, the north, east, and south cells inactive: below the east
    // half of the center cell the head at the two ghost nodes could only be
    // projected beyond the active cells (#170)
    std::unique_ptr<ModflowGrid> mfgrid = make_mfgrid(2, 3, 3, 100.0, 100.0);
    mfgrid->name = "mfg";
    std::unique_ptr<QuadTree3D> qtree = make_qtree(mfgrid.get(), {{150.0, 150.0, 1, 0}});
    qtree->name = "qtg";
    for (int row_col : {1, 5, 7})
        qtree->node_group[mfgrid->get_box_index(1, row_col / 3, row_col % 3)]->active = false;
    qtree->number_nodes();
    qtree->refresh();

    std::vector<std::string> warnings = usg_data_warnings(qtree.get());
    ASSERT_EQ(warnings.size(), 1u);
    EXPECT_NE(warnings.front().find("GRID_TO_USGDATA block usg, 2 connections"), std::string::npos)
        << warnings.front();
    EXPECT_NE(warnings.front().find("next to an inactive cell"), std::string::npos);
    EXPECT_NE(warnings.front().find("XT3D"), std::string::npos);
}

TEST(DefinitionExporter, NoGhostNodeWarningWithoutAnActiveDomain)
{
    std::unique_ptr<ModflowGrid> mfgrid = make_mfgrid(2, 3, 3, 100.0, 100.0);
    mfgrid->name = "mfg";
    std::unique_ptr<QuadTree3D> qtree = make_qtree(mfgrid.get(), {{150.0, 150.0, 1, 0}});
    qtree->name = "qtg";
    EXPECT_TRUE(usg_data_warnings(qtree.get()).empty());
}

namespace
{

// the warnings of exporting the DISV data of grid, with vertical pass-through
// cells if vertical_pass_through is true, to a prefix in the temporary folder
std::vector<std::string> disv_data_warnings(Grid * grid, bool vertical_pass_through)
{
    std::string prefix =
        (std::filesystem::path(::testing::TempDir()) / "gridgen_unit_disv_warnings")
            .generic_string();
    std::istringstream dfn("BEGIN GRID_TO_DISVDATA disv\n"
                           "  GRID = qtg\n"
                           "  DISV_DATA_PREFIX = " +
                           prefix + "\n  VERTICAL_PASS_THROUGH = " +
                           (vertical_pass_through ? "True" : "False") + "\nEND GRID_TO_DISVDATA\n");
    DefinitionParser parser;
    EXPECT_TRUE(parser.parse(dfn));
    auto * disv = dynamic_cast<GridToDisvdataBlock *>(parser.get_block_data("disv"));
    EXPECT_NE(disv, nullptr);
    if (disv == nullptr) return {};
    disv->grid = grid;
    clear_warnings();
    DefinitionExporter exporter;
    exporter.save(*disv);
    std::vector<std::string> warnings = get_warnings();
    clear_warnings();
    return warnings;
}

} // namespace

TEST(DefinitionExporter, DisvDataWithoutInactiveCellsHasNoWarning)
{
    std::unique_ptr<ModflowGrid> mfgrid = make_mfgrid(2, 3, 3, 100.0, 100.0);
    std::unique_ptr<QuadTree3D> qtree =
        make_qtree(mfgrid.get(), {{150.0, 150.0, 1, 0}, {150.0, 150.0, 1, 1}});
    qtree->name = "qtg";
    EXPECT_TRUE(disv_data_warnings(qtree.get(), false).empty());
}

TEST(DefinitionExporter, DisvPassThroughWithoutInactiveCellsWarns)
{
    std::unique_ptr<ModflowGrid> mfgrid = make_mfgrid(2, 3, 3, 100.0, 100.0);
    std::vector<std::string> warnings = disv_data_warnings(mfgrid.get(), true);
    ASSERT_EQ(warnings.size(), 1u);
    EXPECT_NE(warnings.front().find("VERTICAL_PASS_THROUGH in GRID_TO_DISVDATA block disv"),
              std::string::npos)
        << warnings.front();
    EXPECT_NE(warnings.front().find("no inactive cells"), std::string::npos);
}

TEST(DefinitionExporter, DisvTopsThatDifferFromTheBottomAboveAreReported)
{
    // the center cell of layer 2 starts 1 below the bottom of the cell above it
    std::unique_ptr<ModflowGrid> mfgrid = make_mfgrid(2, 3, 3, 100.0, 100.0);
    std::unique_ptr<QuadTree3D> qtree = make_qtree(mfgrid.get(), {});
    qtree->name = "qtg";
    Box * center = qtree->node_group[mfgrid->get_box_index(1, 1, 1)];
    qtree->top[center->number] -= 1.0;
    std::vector<std::string> warnings = disv_data_warnings(qtree.get(), false);
    ASSERT_EQ(warnings.size(), 1u);
    EXPECT_NE(warnings.front().find("GRID_TO_DISVDATA block disv, 1 active cell has a top"),
              std::string::npos)
        << warnings.front();
}

TEST(DefinitionExporter, DisvLayersRefinedDifferentlyNameTheBlock)
{
    std::unique_ptr<ModflowGrid> mfgrid = make_mfgrid(2, 3, 3, 100.0, 100.0);
    std::unique_ptr<QuadTree3D> qtree = make_qtree(mfgrid.get(), {{150.0, 150.0, 1, 0}});
    qtree->name = "qtg";
    try
    {
        disv_data_warnings(qtree.get(), false);
        ADD_FAILURE() << "no error";
    }
    catch (const GridgenError & error)
    {
        std::string message = error.what();
        EXPECT_EQ(message.rfind("In GRID_TO_DISVDATA block disv, layer 2 of grid qtg", 0), 0u)
            << message;
    }
}
