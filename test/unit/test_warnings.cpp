// Warnings and errors for input that is ignored or has no effect
#include <cstdio>
#include <sstream>
#include <string>
#include <vector>

#include <gtest/gtest.h>

#include "GridDomain.h"
#include "GridIntersection.h"
#include "GridgenError.h"
#include "builder_helpers.h"
#include "Warnings.h"
#include "def_parser.h"
#include "def_struct.h"
#include "grid_helpers.h"
#include "shapefile/shpReader.h"

using namespace cusg;
using gridgen_test::build_qtree;
using gridgen_test::builder_dfn;
using gridgen_test::make_mfgrid;
using gridgen_test::make_qtree;
using gridgen_test::square;
using gridgen_test::Vertices;
using gridgen_test::write_test_shapefile;

namespace
{

// whether any recorded warning contains text
bool warned(const std::string & text)
{
    for (const std::string & w : get_warnings())
    {
        if (w.find(text) != std::string::npos) return true;
    }
    return false;
}

} // namespace

TEST(Warnings, FeaturesOutsideAreCounted)
{
    clear_warnings();
    warn_features_outside(0, 3, "a", "ACTIVE_DOMAIN d", 0, "no cells are made active by them");
    EXPECT_TRUE(get_warnings().empty());

    warn_features_outside(1, 1, "a", "ACTIVE_DOMAIN d", 0, "no cells are made active by them");
    warn_features_outside(2, 2, "a", "ACTIVE_DOMAIN d", 1, "no cells are made active by them");
    warn_features_outside(1, 3, "a", "ACTIVE_DOMAIN d", 0, "no cells are made active by them");
    ASSERT_EQ(get_warnings().size(), 3u);
    EXPECT_EQ(get_warnings()[0], "The feature in shapefile a (ACTIVE_DOMAIN d) does not intersect "
                                 "layer 1 of the grid; no cells are made active by it.");
    EXPECT_EQ(get_warnings()[1], "None of the 2 features in shapefile a (ACTIVE_DOMAIN d) "
                                 "intersect layer 2 of the grid; no cells are made active by "
                                 "them.");
    EXPECT_EQ(get_warnings()[2], "1 of the 3 features in shapefile a (ACTIVE_DOMAIN d) do not "
                                 "intersect layer 1 of the grid; no cells are made active by it.");

    clear_warnings();
    EXPECT_TRUE(get_warnings().empty());
}

TEST(QuadTree, RefineReportsFeatureOutsideGrid)
{
    std::unique_ptr<QuadTree3D> qtree = make_qtree(make_mfgrid(1, 3, 3, 100.0, 100.0), {});

    Point2d inside(150.0, 150.0), outside(450.0, 150.0);
    EXPECT_TRUE(qtree->refine(inside, 1));
    EXPECT_FALSE(qtree->refine(outside, 1));

    Polyline line_in, line_out;
    line_in.beginPoly();
    line_in.addVertex(10.0, 10.0);
    line_in.addVertex(290.0, 10.0);
    line_in.endPoly();
    line_out.beginPoly();
    line_out.addVertex(400.0, 10.0);
    line_out.addVertex(500.0, 10.0);
    line_out.endPoly();
    EXPECT_TRUE(qtree->refine(line_in, 1));
    EXPECT_FALSE(qtree->refine(line_out, 1));
    line_in.destroy();
    line_out.destroy();

    // a cell already at the level is still intersected
    EXPECT_TRUE(qtree->refine(inside, 1));

    // a line that crosses the grid with both ends outside intersects it
    Polyline line_across;
    line_across.beginPoly();
    line_across.addVertex(-50.0, 250.0);
    line_across.addVertex(350.0, 250.0);
    line_across.endPoly();
    EXPECT_TRUE(qtree->refine(line_across, 1));
    line_across.destroy();
}

TEST(QuadTreeBuilder, RefinementFeatureOutsideGridIsReported)
{
    std::string path = write_test_shapefile("gridgen_unit_refine_outside", SHPT_POLYGON,
                                            {square(50.0, 50.0, 20.0), square(500.0, 500.0, 20.0)});
    clear_warnings();
    build_qtree(builder_dfn(path, "Polygon"));
    EXPECT_TRUE(warned("1 of the 2 features in shapefile " + path +
                       " (REFINEMENT_FEATURES rf) do not intersect layer 1 of the grid; the grid "
                       "is not refined for it."));
}

TEST(QuadTreeBuilder, RefinementShapefileWithoutFeaturesOfTypeIsReported)
{
    std::string path =
        write_test_shapefile("gridgen_unit_refine_type", SHPT_POLYGON, {square(50.0, 50.0, 20.0)});
    clear_warnings();
    build_qtree(builder_dfn(path, "line"));
    EXPECT_TRUE(warned("has no line features"));
}

TEST(QuadTreeBuilder, RefinementErrors)
{
    std::string path = write_test_shapefile("gridgen_unit_refine_errors", SHPT_POLYGON,
                                            {square(50.0, 50.0, 20.0)});
    std::string missing = testing::TempDir() + "gridgen_unit_no_such_shapefile";
    EXPECT_THROW(build_qtree(builder_dfn(missing, "polygon")), GridgenError);
    EXPECT_THROW(
        build_qtree(builder_dfn(path, "polygon", "REFINEMENT_LEVEL_BY_ATTRIBUTE = NOSUCH")),
        GridgenError);
    EXPECT_THROW(build_qtree(builder_dfn(path, "area")), GridgenError);
    // the field exists, so this refines
    build_qtree(builder_dfn(path, "polygon", "REFINEMENT_LEVEL_BY_ATTRIBUTE = LEVEL"));
}

// lines from beside the 3x3 grid of 100 m cells across it, within one cell, and beside the grid
std::vector<Vertices> lines_one_outside()
{
    return {{{-50.0, 10.0}, {290.0, 10.0}},
            {{10.0, 10.0}, {20.0, 10.0}},
            {{400.0, 10.0}, {500.0, 10.0}}};
}

TEST(GridDomain, FeatureOutsideGridIsReported)
{
    std::string path = write_test_shapefile("gridgen_unit_domain_outside", SHPT_POLYGON,
                                            {square(50.0, 50.0, 20.0), square(500.0, 500.0, 20.0)});
    std::unique_ptr<QuadTree3D> qtree = make_qtree(make_mfgrid(1, 3, 3, 100.0, 100.0), {});

    DefinitionBlock block;
    block.name = "domain";
    ActiveDomainBlock domain(block);
    domain.shapefile = path;
    domain.feature_type = "polygon";

    clear_warnings();
    GridDomain().build(qtree.get(), 0, &domain);
    EXPECT_TRUE(warned("1 of the 2 features in shapefile " + path +
                       " (ACTIVE_DOMAIN domain) do not intersect layer 1 of the grid; no cells "
                       "are made active by it."));

    // a line within one cell is inside; one beside the grid is not
    domain.shapefile =
        write_test_shapefile("gridgen_unit_domain_lines", SHPT_ARC, lines_one_outside());
    domain.feature_type = "line";
    clear_warnings();
    GridDomain().build(qtree.get(), 0, &domain);
    EXPECT_TRUE(warned("1 of the 3 features in shapefile " + domain.shapefile +
                       " (ACTIVE_DOMAIN domain) do not intersect layer 1 of the grid; no cells "
                       "are made active by it."));

    domain.shapefile = testing::TempDir() + "gridgen_unit_no_such_shapefile";
    EXPECT_THROW(GridDomain().build(qtree.get(), 0, &domain), GridgenError);
}

TEST(GridIntersection, FeatureOutsideGridIsReported)
{
    std::string path = write_test_shapefile("gridgen_unit_intersect_outside", SHPT_POINT,
                                            {{{150.0, 150.0}}, {{450.0, 150.0}}});
    std::unique_ptr<QuadTree3D> qtree = make_qtree(make_mfgrid(1, 3, 3, 100.0, 100.0), {});

    DefinitionBlock block;
    block.name = "points";
    GridIntersectionBlock data(block);
    data.grid = qtree.get();
    data.shapefile = path;
    data.feature_type = "point";

    clear_warnings();
    GridIntersection gi(0, 0);
    gi.build(data);
    EXPECT_TRUE(warned("1 of the 2 features in shapefile " + path +
                       " (GRID_INTERSECTION points) do not intersect layer 1 of the grid; no "
                       "intersections are written for it."));

    data.shapefile =
        write_test_shapefile("gridgen_unit_intersect_lines", SHPT_ARC, lines_one_outside());
    data.feature_type = "line";
    clear_warnings();
    GridIntersection gi_lines(0, 0);
    gi_lines.build(data);
    EXPECT_TRUE(warned("1 of the 3 features in shapefile " + data.shapefile +
                       " (GRID_INTERSECTION points) do not intersect layer 1 of the grid; no "
                       "intersections are written for it."));

    data.shapefile = testing::TempDir() + "gridgen_unit_no_such_shapefile";
    GridIntersection gi2(0, 0);
    EXPECT_THROW(gi2.build(data), GridgenError);
}

TEST(ShapeReader, ShapefileWithoutFieldsIsRead)
{
    std::string path =
        write_test_shapefile("gridgen_unit_no_fields", SHPT_POLYGON,
                             {square(0.0, 0.0, 10.0), square(20.0, 0.0, 10.0)}, false);
    ShapeReader reader;
    ASSERT_TRUE(reader.read(path));
    EXPECT_EQ(reader.getPolyList().size(), 2u);
}

TEST(DefinitionParser, FeatureTypeIsCheckedAndCaseInsensitive)
{
    const std::string block = "BEGIN ACTIVE_DOMAIN domain\n"
                              "  SHAPEFILE = domain\n"
                              "  FEATURE_TYPE = ";
    {
        std::istringstream dfn(block + "Polygon\n  INCLUDE_BOUNDARY = True\nEND ACTIVE_DOMAIN\n");
        DefinitionParser parser;
        ASSERT_TRUE(parser.parse(dfn));
        auto * data = dynamic_cast<ActiveDomainBlock *>(parser.getBlockData("domain"));
        ASSERT_NE(data, nullptr);
        EXPECT_EQ(data->feature_type, "polygon");
    }
    {
        std::istringstream dfn(block + "area\nEND ACTIVE_DOMAIN\n");
        DefinitionParser parser;
        EXPECT_THROW(parser.parse(dfn), GridgenError);
    }
}

// CLOSE_POLYGONS of a GRID_TO_SHAPEFILE block given the records
bool close_polygons(const std::string & records)
{
    std::istringstream dfn("BEGIN GRID_TO_SHAPEFILE shp\n"
                           "  GRID = qtg\n"
                           "  SHAPEFILE = qtg_shp\n"
                           "  FEATURE_TYPE = polygon\n" +
                           records + "END GRID_TO_SHAPEFILE\n");
    DefinitionParser parser;
    EXPECT_TRUE(parser.parse(dfn));
    auto * data = dynamic_cast<GridToShapefileBlock *>(parser.getBlockData("shp"));
    EXPECT_NE(data, nullptr);
    return data != nullptr && data->close_polygons;
}

TEST(DefinitionParser, ConcideIsADeprecatedNameOfClosePolygons)
{
    const std::string deprecated =
        "CONCIDE (GRID_TO_SHAPEFILE shp) is deprecated; use CLOSE_POLYGONS.";

    clear_warnings();
    EXPECT_TRUE(close_polygons(""));
    EXPECT_TRUE(close_polygons("  CLOSE_POLYGONS = True\n"));
    EXPECT_FALSE(close_polygons("  CLOSE_POLYGONS = False\n"));
    EXPECT_TRUE(get_warnings().empty());

    clear_warnings();
    EXPECT_FALSE(close_polygons("  CONCIDE = False\n"));
    EXPECT_TRUE(warned(deprecated));

    // CLOSE_POLYGONS is used when both are given
    clear_warnings();
    EXPECT_TRUE(close_polygons("  CONCIDE = False\n  CLOSE_POLYGONS = True\n"));
    EXPECT_TRUE(close_polygons("  CLOSE_POLYGONS = True\n  CONCIDE = False\n"));
    EXPECT_TRUE(warned(deprecated));
}

TEST(DefinitionParser, IgnoredRecordsAndBlocksAreReported)
{
    std::string block = "BEGIN ACTIVE_DOMAIN domain\n"
                        "  SHAPEFILE = domain\n"
                        "  FEATURE_TYPE = polygon\n";
    std::istringstream dfn(block + "  NO_SUCH_RECORD = 1\nEND ACTIVE_DOMAIN\n" + block +
                           "END ACTIVE_DOMAIN\nBEGIN NO_SUCH_TYPE other\nEND NO_SUCH_TYPE\n");
    DefinitionParser parser;
    clear_warnings();
    ASSERT_TRUE(parser.parse(dfn));
    ASSERT_EQ(get_warnings().size(), 3u);
    EXPECT_NE(get_warnings()[0].find("NO_SUCH_RECORD"), std::string::npos);
    EXPECT_TRUE(
        warned("Block domain is defined more than once; only the last definition is used."));
    EXPECT_TRUE(warned("Block other has the unrecognized type no_such_type and is ignored."));
}

TEST(DefinitionParser, MalformedInputIsAnError)
{
    {
        // a record with no value
        std::istringstream dfn("BEGIN ACTIVE_DOMAIN domain\n"
                               "  SHAPEFILE =\n"
                               "  FEATURE_TYPE = polygon\n"
                               "END ACTIVE_DOMAIN\n");
        DefinitionParser parser;
        EXPECT_THROW(parser.parse(dfn), GridgenError);
    }
    {
        std::istringstream dfn("BEGIN ACTIVE_DOMAIN domain\n"
                               "  SHAPEFILE = domain\n"
                               "  FEATURE_TYPE = polygon\n"
                               "END REFINEMENT_FEATURES\n");
        DefinitionParser parser;
        EXPECT_THROW(parser.parse(dfn), GridgenError);
    }
    {
        DefinitionParser parser;
        EXPECT_THROW(parser.parse(testing::TempDir() + "gridgen_unit_no_such.dfn"), GridgenError);
    }
}

TEST(DefinitionParser, ShortArrayFileIsAnError)
{
    std::string path = testing::TempDir() + "gridgen_unit_short.dat";
    FILE * f = std::fopen(path.c_str(), "w");
    ASSERT_NE(f, nullptr);
    std::fputs("50. 50. 50.\n", f);
    std::fclose(f);

    std::istringstream dfn("BEGIN MODFLOW_GRID mfg\n"
                           "  ROTATION_ANGLE = 0\n"
                           "  X_OFFSET = 0\n"
                           "  Y_OFFSET = 0\n"
                           "  NLAY = 1\n"
                           "  NROW = 3\n"
                           "  NCOL = 4\n"
                           "  DELR = OPEN/CLOSE " +
                           path +
                           "\n"
                           "  DELC = CONSTANT 25.\n"
                           "  TOP = CONSTANT 10.\n"
                           "  BOTTOM LAYER 1 = CONSTANT 0.\n"
                           "END MODFLOW_GRID\n");
    DefinitionParser parser;
    ASSERT_TRUE(parser.parse(dfn));
    DefinitionBlock * data = parser.getBlockData("mfg");
    ASSERT_NE(data, nullptr);
    EXPECT_THROW(data->initialize(parser.getBlockData()), GridgenError);
}

// the cell with the largest negative thickness comes before a thinner one
TEST(QuadTreeBuilder, NegativeThicknessNamesTheThickestCell)
{
    std::string path = testing::TempDir() + "gridgen_unit_negative_bottom.dat";
    FILE * f = std::fopen(path.c_str(), "w");
    ASSERT_NE(f, nullptr);
    std::fputs("15. 0. 0.\n0. 12. 0.\n0. 0. 0.\n", f);
    std::fclose(f);

    std::string dfn = "BEGIN MODFLOW_GRID mfg\n"
                      "  ROTATION_ANGLE = 0\n"
                      "  X_OFFSET = 0\n"
                      "  Y_OFFSET = 0\n"
                      "  NLAY = 1\n"
                      "  NROW = 3\n"
                      "  NCOL = 3\n"
                      "  DELR = CONSTANT 100.\n"
                      "  DELC = CONSTANT 100.\n"
                      "  TOP = CONSTANT 10.\n"
                      "  BOTTOM LAYER 1 = OPEN/CLOSE " +
                      path +
                      "\n"
                      "END MODFLOW_GRID\n"
                      "BEGIN QUADTREE_BUILDER qtb\n"
                      "  MODFLOW_GRID = mfg\n"
                      "  GRID_DEFINITION_FILE = qtg.dfn\n"
                      "END QUADTREE_BUILDER\n";
    clear_warnings();
    build_qtree(dfn);
    ASSERT_EQ(get_warnings().size(), 1u);
    EXPECT_EQ(get_warnings()[0], "There are 2 cells with negative thickness. Node 1 in layer 1 "
                                 "has the largest negative thickness, -5.");
}
