// Field types of the .dbf files of the shapefiles gridgen writes
#include <filesystem>
#include <memory>
#include <sstream>
#include <string>
#include <vector>

#include <gtest/gtest.h>

#include "GridIntersection.h"
#include "GridToShp.h"
#include "dbfFields.h"
#include "def_parser.h"
#include "def_struct.h"
#include "grid_helpers.h"

using namespace cusg;
using gridgen_test::make_mfgrid;
using gridgen_test::make_qtree;

namespace
{

// every numeric field with no decimals in the .dbf file at path is an integer field of the
// width gridgen writes, and the integer field names are those expected
void expect_integer_fields(const std::string & path, const std::vector<std::string> & expected)
{
    DBFHandle dbf = DBFOpen(path.c_str(), "rb");
    ASSERT_NE(dbf, nullptr) << path;
    std::vector<std::string> integer_fields;
    for (int i = 0; i < DBFGetFieldCount(dbf); i++)
    {
        char name[12];
        int width = 0;
        int decimals = 0;
        DBFFieldType type = DBFGetFieldInfo(dbf, i, name, &width, &decimals);
        char native_type = DBFGetNativeFieldType(dbf, i);
        if ((native_type != 'N' && native_type != 'F') || decimals != 0) continue;
        SCOPED_TRACE(name);
        EXPECT_EQ(type, FTInteger);
        EXPECT_EQ(width, DBF_INTEGER_WIDTH);
        integer_fields.push_back(name);
    }
    EXPECT_EQ(integer_fields, expected);
    DBFClose(dbf);
}

// the integer value of the field named field of record in the .dbf file at path
int integer_value(const std::string & path, int record, const char * field)
{
    DBFHandle dbf = DBFOpen(path.c_str(), "rb");
    if (dbf == nullptr) return -1;
    int value = DBFReadIntegerAttribute(dbf, record, DBFGetFieldIndex(dbf, field));
    DBFClose(dbf);
    return value;
}

} // namespace

TEST(ShapefileFields, GridIntegerFieldsAreIntegers)
{
    // node, layer, row, and column fields are read as integers, not doubles (#188)
    std::string folder = ::testing::TempDir();
    std::unique_ptr<ModflowGrid> mfgrid = make_mfgrid(2, 3, 4, 10.0, 10.0);
    const std::vector<std::string> cell_fields = {"nodenumber", "layer", "row", "col"};
    //the node number is one-based and the layer, row, and column zero-based
    for (const char * feature_type : {"polygon", "point"})
    {
        SCOPED_TRACE(feature_type);
        std::string path = folder + "gridgen_unit_fields_mfg_" + feature_type;
        write_shapefile(mfgrid.get(), path, feature_type);
        expect_integer_fields(path + ".dbf", cell_fields);
        EXPECT_EQ(integer_value(path + ".dbf", 23, "nodenumber"), 24);
        EXPECT_EQ(integer_value(path + ".dbf", 23, "layer"), 1);
        EXPECT_EQ(integer_value(path + ".dbf", 23, "row"), 2);
        EXPECT_EQ(integer_value(path + ".dbf", 23, "col"), 3);
    }
    // a line shapefile has only the layer field
    std::string line_path = folder + "gridgen_unit_fields_mfg_line";
    write_shapefile(mfgrid.get(), line_path, "line");
    expect_integer_fields(line_path + ".dbf", {"layer"});

    //refining the quadtree splits the cells of the MODFLOW grid, so it is built after the
    //MODFLOW grid is written; one cell split in four gives 27 cells
    std::unique_ptr<QuadTree3D> qtree = make_qtree(mfgrid.get(), {{5.0, 25.0, 1, 0}});
    for (const char * feature_type : {"polygon", "point"})
    {
        SCOPED_TRACE(feature_type);
        std::string path = folder + "gridgen_unit_fields_qtg_" + feature_type;
        write_shapefile(qtree.get(), path, feature_type);
        expect_integer_fields(path + ".dbf", cell_fields);
        EXPECT_EQ(integer_value(path + ".dbf", 26, "layer"), 1);
        EXPECT_EQ(integer_value(path + ".dbf", 26, "row"), 2);
        EXPECT_EQ(integer_value(path + ".dbf", 26, "col"), 3);
    }
    std::string qt_line_path = folder + "gridgen_unit_fields_qtg_line";
    write_shapefile(qtree.get(), qt_line_path, "line");
    expect_integer_fields(qt_line_path + ".dbf", {"layer"});
}

TEST(ShapefileFields, IntersectionIntegerFieldsAreIntegers)
{
    // the node number, point id, and an integer attribute copied from the input shapefile are
    // read as integers, and the attribute keeps its value (#188)
    std::string folder = ::testing::TempDir();
    std::string input = folder + "gridgen_unit_fields_wells";
    {
        SHPHandle shp = SHPCreate(input.c_str(), SHPT_POINT);
        DBFHandle dbf = DBFCreate(input.c_str());
        ASSERT_NE(shp, nullptr);
        ASSERT_NE(dbf, nullptr);
        DBFAddField(dbf, "well_id", FTInteger, 10, 0);
        const double x[] = {5.0, 35.0};
        const double y[] = {25.0, 5.0};
        const int well_id[] = {2000000000, -999999999};
        for (int i = 0; i < 2; i++)
        {
            SHPObject * obj = SHPCreateSimpleObject(SHPT_POINT, 1, &x[i], &y[i], nullptr);
            SHPWriteObject(shp, -1, obj);
            SHPDestroyObject(obj);
            DBFWriteIntegerAttribute(dbf, i, 0, well_id[i]);
        }
        SHPClose(shp);
        DBFClose(dbf);
    }

    std::istringstream dfn("BEGIN MODFLOW_GRID mfg\n"
                           "  ROTATION_ANGLE = 0\n"
                           "  X_OFFSET = 0\n"
                           "  Y_OFFSET = 0\n"
                           "  NLAY = 1\n"
                           "  NROW = 3\n"
                           "  NCOL = 4\n"
                           "  DELR = CONSTANT 10.\n"
                           "  DELC = CONSTANT 10.\n"
                           "  TOP = CONSTANT 10.\n"
                           "  BOTTOM LAYER 1 = CONSTANT 0.\n"
                           "END MODFLOW_GRID\n"
                           "BEGIN GRID_INTERSECTION wells\n"
                           "  GRID = mfg\n"
                           "  LAYER = 1\n"
                           "  SHAPEFILE = " +
                           input +
                           "\n"
                           "  FEATURE_TYPE = point\n"
                           "  OUTPUT_FILE = unused.txt\n"
                           "END GRID_INTERSECTION\n");
    DefinitionParser parser;
    ASSERT_TRUE(parser.parse(dfn));
    GridIntersectionBlock * data =
        dynamic_cast<GridIntersectionBlock *>(parser.get_block_data("wells"));
    ASSERT_NE(data, nullptr);
    for (const char * name : {"mfg", "wells"})
        parser.get_block_data(name)->initialize(parser.get_block_data());

    GridIntersection intersection(0, data->layer);
    intersection.set_save_intersection_info(true);
    intersection.build(*data);
    std::string output = folder + "gridgen_unit_fields_wells_out";
    ASSERT_TRUE(intersection.save_shp(output, data->attributes));

    std::string dbf_path = output + ".dbf";
    expect_integer_fields(dbf_path, {"nodenumber", "pointid", "well_id"});
    // the first well is in row 1, column 1 and the second in row 3, column 4
    EXPECT_EQ(integer_value(dbf_path, 0, "nodenumber"), 1);
    EXPECT_EQ(integer_value(dbf_path, 0, "well_id"), 2000000000);
    EXPECT_EQ(integer_value(dbf_path, 1, "nodenumber"), 12);
    EXPECT_EQ(integer_value(dbf_path, 1, "well_id"), -999999999);
}
