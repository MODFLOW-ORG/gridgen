// File names in the grid definition files that gridgen writes
#include <filesystem>
#include <memory>
#include <sstream>
#include <string>

#include <gtest/gtest.h>

#include "def_exporter.h"
#include "def_parser.h"
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
