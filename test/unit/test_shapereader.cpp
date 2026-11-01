// Shapes read from shapefiles
#include <string>
#include <vector>

#include <gtest/gtest.h>

#include "shapefile/shpReader.h"

using namespace cusg;

namespace
{

struct Part
{
    int type;
    std::vector<double> x;
    std::vector<double> y;
};

// a closed square ring with its lower left corner at (x0, y0)
Part square(int type, double x0, double y0, double size)
{
    return {type,
            {x0, x0, x0 + size, x0 + size, x0},
            {y0, y0 + size, y0 + size, y0, y0}};
}

// write a multipatch shapefile with one shape for each list of parts
void write_multipatch(const std::string & path, const std::vector<std::vector<Part>> & shapes)
{
    SHPHandle shp = SHPCreate(path.c_str(), SHPT_MULTIPATCH);
    DBFHandle dbf = DBFCreate(path.c_str());
    ASSERT_NE(shp, nullptr);
    ASSERT_NE(dbf, nullptr);
    DBFAddField(dbf, "ID", FTInteger, 4, 0);
    for (size_t i = 0; i < shapes.size(); i++)
    {
        std::vector<int> starts, types;
        std::vector<double> x, y;
        for (const Part & p : shapes[i])
        {
            starts.push_back(static_cast<int>(x.size()));
            types.push_back(p.type);
            x.insert(x.end(), p.x.begin(), p.x.end());
            y.insert(y.end(), p.y.begin(), p.y.end());
        }
        SHPObject * obj = SHPCreateObject(SHPT_MULTIPATCH, -1, static_cast<int>(starts.size()),
                                          starts.data(), types.data(), static_cast<int>(x.size()),
                                          x.data(), y.data(), nullptr, nullptr);
        SHPWriteObject(shp, -1, obj);
        SHPDestroyObject(obj);
        DBFWriteIntegerAttribute(dbf, static_cast<int>(i), 0, static_cast<int>(i + 1));
    }
    SHPClose(shp);
    DBFClose(dbf);
}

} // namespace

TEST(ShapeReader, MultipatchRingsAreReadAsPolygonsWithHoles)
{
    // outer and first rings start a polygon, inner rings and a ring after
    // another ring are holes, and triangle strips are skipped (#100)
    std::string path = testing::TempDir() + "gridgen_unit_multipatch";
    write_multipatch(path, {{square(SHPP_OUTERRING, 0, 0, 10), square(SHPP_INNERRING, 3, 3, 4)},
                            {square(SHPP_FIRSTRING, 20, 0, 10), square(SHPP_RING, 23, 3, 4)},
                            {square(SHPP_RING, 40, 0, 10)},
                            {{SHPP_TRISTRIP, {60, 70, 60, 70}, {0, 0, 10, 10}},
                             square(SHPP_OUTERRING, 60, 0, 10)}});

    ShapeReader reader;
    ASSERT_TRUE(reader.read(path));
    std::vector<GIS_polygon> & polygons = reader.getPolyList();
    ASSERT_EQ(polygons.size(), 4u);

    const size_t rings[] = {2, 2, 1, 1};
    for (size_t i = 0; i < polygons.size(); i++)
    {
        EXPECT_EQ(polygons[i].m_id, static_cast<int>(i));
        ASSERT_EQ(polygons[i].size(), rings[i]) << "shape " << i;
        EXPECT_EQ(polygons[i].front().getType(), c_ply::POUT) << "shape " << i;
        if (rings[i] == 2) EXPECT_EQ(polygons[i].back().getType(), c_ply::PIN) << "shape " << i;
    }
    EXPECT_TRUE(reader.getArcList().empty());
}
