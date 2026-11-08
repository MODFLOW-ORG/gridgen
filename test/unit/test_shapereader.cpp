// Shapes read from shapefiles
#include <map>
#include <set>
#include <string>
#include <utility>
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
    return {type, {x0, x0, x0 + size, x0 + size, x0}, {y0, y0 + size, y0 + size, y0, y0}};
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

// write a polygon shapefile with one shape made of the rings
void write_polygon(const std::string & path, const std::vector<Part> & rings)
{
    SHPHandle shp = SHPCreate(path.c_str(), SHPT_POLYGON);
    DBFHandle dbf = DBFCreate(path.c_str());
    ASSERT_NE(shp, nullptr);
    ASSERT_NE(dbf, nullptr);
    DBFAddField(dbf, "ID", FTInteger, 4, 0);
    std::vector<int> starts;
    std::vector<double> x, y;
    for (const Part & p : rings)
    {
        starts.push_back(static_cast<int>(x.size()));
        x.insert(x.end(), p.x.begin(), p.x.end());
        y.insert(y.end(), p.y.begin(), p.y.end());
    }
    SHPObject * obj =
        SHPCreateObject(SHPT_POLYGON, -1, static_cast<int>(starts.size()), starts.data(), nullptr,
                        static_cast<int>(x.size()), x.data(), y.data(), nullptr, nullptr);
    SHPWriteObject(shp, -1, obj);
    SHPDestroyObject(obj);
    DBFWriteIntegerAttribute(dbf, 0, 0, 1);
    SHPClose(shp);
    DBFClose(dbf);
}

// the lower left corner of a ring
std::pair<double, double> lower_left(c_ply & ring)
{
    double x = ring[0]->getPos()[0], y = ring[0]->getPos()[1];
    for (int i = 1; i < ring.getSize(); i++)
    {
        x = std::min(x, ring[i]->getPos()[0]);
        y = std::min(y, ring[i]->getPos()[1]);
    }
    return {x, y};
}

} // namespace

TEST(ShapeReader, MultipatchRingsAreReadAsPolygonsWithHoles)
{
    // outer and first rings start a polygon, inner rings and a ring after
    // another ring are holes, and triangle strips and fans, and an inner ring
    // before any other ring, are skipped (#100)
    std::string path = testing::TempDir() + "gridgen_unit_multipatch";
    write_multipatch(path, {{square(SHPP_OUTERRING, 0, 0, 10), square(SHPP_INNERRING, 3, 3, 4)},
                            {square(SHPP_FIRSTRING, 20, 0, 10), square(SHPP_RING, 23, 3, 4)},
                            {square(SHPP_RING, 40, 0, 10)},
                            {{SHPP_TRISTRIP, {60, 70, 60, 70}, {0, 0, 10, 10}},
                             square(SHPP_OUTERRING, 60, 0, 10)},
                            {square(SHPP_OUTERRING, 80, 0, 10), square(SHPP_OUTERRING, 100, 0, 10)},
                            {{SHPP_TRIFAN, {120, 130, 130}, {0, 0, 10}},
                             square(SHPP_INNERRING, 123, 3, 4),
                             square(SHPP_RING, 120, 0, 10),
                             square(SHPP_RING, 123, 3, 4)}});

    ShapeReader reader;
    ASSERT_TRUE(reader.read(path));
    std::vector<GIS_polygon> & polygons = reader.getPolyList();
    ASSERT_EQ(polygons.size(), 7u);

    // the shape each polygon is from, and its rings
    const int shapes[] = {0, 1, 2, 3, 4, 4, 5};
    const size_t rings[] = {2, 2, 1, 1, 1, 1, 2};
    for (size_t i = 0; i < polygons.size(); i++)
    {
        EXPECT_EQ(polygons[i].m_id, shapes[i]);
        ASSERT_EQ(polygons[i].size(), rings[i]) << "polygon " << i;
        EXPECT_EQ(polygons[i].front().getType(), c_ply::POUT) << "polygon " << i;
        if (rings[i] == 2) EXPECT_EQ(polygons[i].back().getType(), c_ply::PIN) << "polygon " << i;
    }
    EXPECT_TRUE(reader.getArcList().empty());
}

TEST(ShapeReader, PolygonWithManyRingsFindsEachHole)
{
    // a 30 by 30 block of outer rings, two of every three with a hole and an
    // island in the hole, written in scrambled order, so each outer ring is a
    // polygon with its own hole and each island is a polygon, in the order
    // they were written
    const int side = 30;
    std::vector<Part> rings;
    std::map<std::pair<double, double>, std::pair<double, double>> hole_of;
    std::set<std::pair<double, double>> islands;
    for (int i = 0; i < side * side; i++)
    {
        double x0 = (i % side) * 10.0, y0 = (i / side) * 10.0;
        rings.push_back(square(SHPP_RING, x0, y0, 9));
        if (i % 3 == 0) continue;
        rings.push_back(square(SHPP_RING, x0 + 2, y0 + 2, 5));
        rings.push_back(square(SHPP_RING, x0 + 4, y0 + 4, 1));
        hole_of[{x0, y0}] = {x0 + 2, y0 + 2};
        islands.insert({x0 + 4, y0 + 4});
    }
    std::vector<Part> scrambled;
    std::vector<std::pair<double, double>> polygon_order; // outer rings and islands, as written
    for (size_t i = 0; i < rings.size(); i++)
    {
        scrambled.push_back(rings[(i * 7919) % rings.size()]);
        const Part & ring = scrambled.back();
        if (ring.x[2] - ring.x[0] != 5) polygon_order.push_back({ring.x[0], ring.y[0]});
    }
    std::string path = testing::TempDir() + "gridgen_unit_many_rings";
    write_polygon(path, scrambled);

    ShapeReader reader;
    ASSERT_TRUE(reader.read(path));
    std::vector<GIS_polygon> & polygons = reader.getPolyList();
    ASSERT_EQ(polygons.size(), static_cast<size_t>(side * side) + islands.size());
    size_t nholes = 0, nislands = 0;
    for (size_t i = 0; i < polygons.size(); i++)
    {
        GIS_polygon & polygon = polygons[i];
        auto corner = lower_left(polygon.front());
        EXPECT_EQ(corner, polygon_order[i]) << "polygon " << i;
        EXPECT_EQ(polygon.front().getType(), c_ply::POUT);
        if (islands.count(corner))
        {
            EXPECT_EQ(polygon.size(), 1u);
            nislands++;
            continue;
        }
        auto hole = hole_of.find(corner);
        if (hole == hole_of.end())
        {
            EXPECT_EQ(polygon.size(), 1u);
            continue;
        }
        ASSERT_EQ(polygon.size(), 2u);
        EXPECT_EQ(polygon.back().getType(), c_ply::PIN);
        EXPECT_EQ(lower_left(polygon.back()), hole->second);
        nholes++;
    }
    EXPECT_EQ(nholes, hole_of.size());
    EXPECT_EQ(nislands, islands.size());
}
