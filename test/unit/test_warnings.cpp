// Warnings and errors for input that is ignored or has no effect
#include <cstdio>
#include <sstream>
#include <string>
#include <vector>

#include <gtest/gtest.h>

#include "GridDomain.h"
#include "GridIntersection.h"
#include "GridgenError.h"
#include "Warnings.h"
#include "def_parser.h"
#include "def_struct.h"
#include "grid_helpers.h"
#include "shapefile/shpReader.h"

using namespace cusg;
using gridgen_test::make_mfgrid;
using gridgen_test::make_qtree;

namespace
{

using Ring = std::vector<std::pair<double, double>>;

// a closed square ring with its lower left corner at (x0, y0)
Ring square(double x0, double y0, double size)
{
    return {{x0, y0}, {x0, y0 + size}, {x0 + size, y0 + size}, {x0 + size, y0}, {x0, y0}};
}

// write a shapefile of one shape for each ring, each with an integer field
// LEVEL of 1 unless fields is false
std::string write_shapefile(const std::string & name, int type, const std::vector<Ring> & rings,
                            bool fields = true)
{
    std::string path = testing::TempDir() + name;
    SHPHandle shp = SHPCreate(path.c_str(), type);
    DBFHandle dbf = DBFCreate(path.c_str());
    if (fields) DBFAddField(dbf, "LEVEL", FTInteger, 4, 0);
    for (size_t i = 0; i < rings.size(); i++)
    {
        std::vector<double> x, y;
        for (const auto & p : rings[i])
        {
            x.push_back(p.first);
            y.push_back(p.second);
        }
        SHPObject * obj =
            SHPCreateSimpleObject(type, static_cast<int>(x.size()), x.data(), y.data(), nullptr);
        SHPWriteObject(shp, -1, obj);
        SHPDestroyObject(obj);
        if (fields) DBFWriteIntegerAttribute(dbf, static_cast<int>(i), 0, 1);
    }
    SHPClose(shp);
    DBFClose(dbf);
    return path;
}

// whether any recorded warning contains text
bool warned(const std::string & text)
{
    for (const std::string & w : get_warnings())
    {
        if (w.find(text) != std::string::npos) return true;
    }
    return false;
}

// a quadtree builder over a 3 x 3 grid of 100 m cells from (0, 0), refined by
// a shapefile
std::string builder_dfn(const std::string & shapefile, const std::string & feature_type,
                        const std::string & level = "REFINEMENT_LEVEL = 1")
{
    return "BEGIN MODFLOW_GRID mfg\n"
           "  ROTATION_ANGLE = 0\n"
           "  X_OFFSET = 0\n"
           "  Y_OFFSET = 0\n"
           "  NLAY = 1\n"
           "  NROW = 3\n"
           "  NCOL = 3\n"
           "  DELR = CONSTANT 100.\n"
           "  DELC = CONSTANT 100.\n"
           "  TOP = CONSTANT 10.\n"
           "  BOTTOM LAYER 1 = CONSTANT 0.\n"
           "END MODFLOW_GRID\n"
           "BEGIN QUADTREE_BUILDER qtb\n"
           "  MODFLOW_GRID = mfg\n"
           "  REFINEMENT_FEATURES LAYER 1 = rf\n"
           "  GRID_DEFINITION_FILE = qtg.dfn\n"
           "END QUADTREE_BUILDER\n"
           "BEGIN REFINEMENT_FEATURES rf\n"
           "  SHAPEFILE = " +
           shapefile + "\n  FEATURE_TYPE = " + feature_type + "\n  " + level +
           "\nEND REFINEMENT_FEATURES\n";
}

// build the quadtree of builder_dfn
void build_qtree(const std::string & dfn_text)
{
    std::istringstream dfn(dfn_text);
    DefinitionParser parser;
    ASSERT_TRUE(parser.parse(dfn));
    auto * data = dynamic_cast<quadtree_builder_raw_data *>(parser.getBlockData("qtb"));
    ASSERT_NE(data, nullptr);
    ASSERT_TRUE(data->initialize(parser.getBlockData()));
    data->build();
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

    c_plyline line_in, line_out;
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
}

TEST(QuadTreeBuilder, RefinementFeatureOutsideGridIsReported)
{
    std::string path = write_shapefile("gridgen_unit_refine_outside", SHPT_POLYGON,
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
        write_shapefile("gridgen_unit_refine_type", SHPT_POLYGON, {square(50.0, 50.0, 20.0)});
    clear_warnings();
    build_qtree(builder_dfn(path, "line"));
    EXPECT_TRUE(warned("has no line features"));
}

TEST(QuadTreeBuilder, RefinementErrors)
{
    std::string path =
        write_shapefile("gridgen_unit_refine_errors", SHPT_POLYGON, {square(50.0, 50.0, 20.0)});
    std::string missing = testing::TempDir() + "gridgen_unit_no_such_shapefile";
    EXPECT_THROW(build_qtree(builder_dfn(missing, "polygon")), GridgenError);
    EXPECT_THROW(
        build_qtree(builder_dfn(path, "polygon", "REFINEMENT_LEVEL_BY_ATTRIBUTE = NOSUCH")),
        GridgenError);
    // the field exists, so this refines
    build_qtree(builder_dfn(path, "polygon", "REFINEMENT_LEVEL_BY_ATTRIBUTE = LEVEL"));
}

TEST(GridDomain, FeatureOutsideGridIsReported)
{
    std::string path = write_shapefile("gridgen_unit_domain_outside", SHPT_POLYGON,
                                       {square(50.0, 50.0, 20.0), square(500.0, 500.0, 20.0)});
    std::unique_ptr<QuadTree3D> qtree = make_qtree(make_mfgrid(1, 3, 3, 100.0, 100.0), {});

    block_raw_data block;
    block.name = "domain";
    active_domain_raw_data domain(block);
    domain.shapefile = path;
    domain.feature_type = "polygon";

    clear_warnings();
    GridDomain().build(qtree.get(), 0, &domain);
    EXPECT_TRUE(warned("1 of the 2 features in shapefile " + path +
                       " (ACTIVE_DOMAIN domain) do not intersect layer 1 of the grid; no cells "
                       "are made active by it."));

    domain.shapefile = testing::TempDir() + "gridgen_unit_no_such_shapefile";
    EXPECT_THROW(GridDomain().build(qtree.get(), 0, &domain), GridgenError);
}

TEST(GridIntersection, FeatureOutsideGridIsReported)
{
    std::string path = write_shapefile("gridgen_unit_intersect_outside", SHPT_POINT,
                                       {{{150.0, 150.0}}, {{450.0, 150.0}}});
    std::unique_ptr<QuadTree3D> qtree = make_qtree(make_mfgrid(1, 3, 3, 100.0, 100.0), {});

    block_raw_data block;
    block.name = "points";
    grid_intersection_raw_data data(block);
    data.grid = qtree.get();
    data.shapefile = path;
    data.feature_type = "point";

    clear_warnings();
    GridIntersection gi(0, 0);
    gi.build(data);
    EXPECT_TRUE(warned("1 of the 2 features in shapefile " + path +
                       " (GRID_INTERSECTION points) do not intersect layer 1 of the grid; no "
                       "intersections are written for it."));

    data.shapefile = testing::TempDir() + "gridgen_unit_no_such_shapefile";
    GridIntersection gi2(0, 0);
    EXPECT_THROW(gi2.build(data), GridgenError);
}

TEST(ShapeReader, ShapefileWithoutFieldsIsRead)
{
    std::string path = write_shapefile("gridgen_unit_no_fields", SHPT_POLYGON,
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
        auto * data = dynamic_cast<active_domain_raw_data *>(parser.getBlockData("domain"));
        ASSERT_NE(data, nullptr);
        EXPECT_EQ(data->feature_type, "polygon");
    }
    {
        std::istringstream dfn(block + "area\nEND ACTIVE_DOMAIN\n");
        DefinitionParser parser;
        EXPECT_THROW(parser.parse(dfn), GridgenError);
    }
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
    block_raw_data * data = parser.getBlockData("mfg");
    ASSERT_NE(data, nullptr);
    EXPECT_THROW(data->initialize(parser.getBlockData()), GridgenError);
}
