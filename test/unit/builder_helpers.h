// Shapefiles written for the unit tests and quadtree builders that read them
#pragma once

#include <cstdio>
#include <sstream>
#include <string>
#include <tuple>
#include <utility>
#include <vector>

#include <gtest/gtest.h>

#include "def_parser.h"
#include "QuadTree3D.h"
#include "def_struct.h"
#include "shapefile/shpReader.h"

namespace gridgen_test
{

using namespace cusg;

using Ring = std::vector<std::pair<double, double>>;

// a closed square ring with its lower left corner at (x0, y0)
inline Ring square(double x0, double y0, double size)
{
    return {{x0, y0}, {x0, y0 + size}, {x0 + size, y0 + size}, {x0 + size, y0}, {x0, y0}};
}

// write a shapefile of one shape for each ring, each with an integer field
// LEVEL of 1 unless fields is false
inline std::string write_shapefile(const std::string & name, int type,
                                   const std::vector<Ring> & rings, bool fields = true)
{
    std::string path = testing::TempDir() + name;
    SHPHandle shp = SHPCreate(path.c_str(), type);
    DBFHandle dbf = nullptr;
    if (fields)
    {
        dbf = DBFCreate(path.c_str());
        DBFAddField(dbf, "LEVEL", FTInteger, 4, 0);
    }
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
    if (fields)
        DBFClose(dbf);
    else
    {
        // shapelib passes a null buffer to fwrite for a table without fields,
        // so write the header and blank records here
        unsigned char header[32] = {0x03, 126, 11, 4};
        header[4] = static_cast<unsigned char>(rings.size());
        header[8] = 33; // header length
        header[10] = 1; // record length, the deletion flag
        FILE * f = std::fopen((path + ".dbf").c_str(), "wb");
        std::fwrite(header, 1, sizeof(header), f);
        std::fputc(0x0D, f);
        for (size_t i = 0; i < rings.size(); i++) std::fputc(' ', f);
        std::fputc(0x1A, f);
        std::fclose(f);
    }
    return path;
}

// a shape of a shapefile written by write_shapes: its parts, none for a null
// shape, and its LEVEL
struct Shape
{
    std::vector<Ring> parts;
    int level;
};

// write a shapefile of type with an integer field LEVEL; the points of a
// point or multipoint shape are its first part, and each part of a multipatch
// shape is an outer ring
inline std::string write_shapes(const std::string & name, int type,
                                const std::vector<Shape> & shapes)
{
    std::string path = testing::TempDir() + name;
    SHPHandle shp = SHPCreate(path.c_str(), type);
    DBFHandle dbf = DBFCreate(path.c_str());
    DBFAddField(dbf, "LEVEL", FTInteger, 4, 0);
    bool has_parts = type == SHPT_ARC || type == SHPT_POLYGON || type == SHPT_MULTIPATCH;
    for (size_t i = 0; i < shapes.size(); i++)
    {
        std::vector<double> x, y;
        std::vector<int> starts, part_types;
        for (const Ring & part : shapes[i].parts)
        {
            starts.push_back(static_cast<int>(x.size()));
            part_types.push_back(SHPP_OUTERRING);
            for (const auto & p : part)
            {
                x.push_back(p.first);
                y.push_back(p.second);
            }
        }
        SHPObject * obj =
            shapes[i].parts.empty()
                ? SHPCreateSimpleObject(SHPT_NULL, 0, nullptr, nullptr, nullptr)
                : SHPCreateObject(type, -1, has_parts ? static_cast<int>(starts.size()) : 0,
                                  starts.data(),
                                  type == SHPT_MULTIPATCH ? part_types.data() : nullptr,
                                  static_cast<int>(x.size()), x.data(), y.data(), nullptr, nullptr);
        SHPWriteObject(shp, -1, obj);
        SHPDestroyObject(obj);
        DBFWriteIntegerAttribute(dbf, static_cast<int>(i), 0, shapes[i].level);
    }
    SHPClose(shp);
    DBFClose(dbf);
    return path;
}

// a quadtree builder over a 3 x 3 grid of 100 m cells from (0, 0), rotated
// by rotation degrees about (0, 0), refined by a shapefile
inline std::string builder_dfn(const std::string & shapefile, const std::string & feature_type,
                               const std::string & level = "REFINEMENT_LEVEL = 1",
                               const std::string & rotation = "0")
{
    return "BEGIN MODFLOW_GRID mfg\n"
           "  ROTATION_ANGLE = " +
           rotation +
           "\n"
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

// build the quadtree of builder_dfn; nodes, if given, is set to its number of
// cells, and cells to the refinement level and center of each cell in node order
inline void build_qtree(const std::string & dfn_text, int * nodes = nullptr,
                        std::vector<std::tuple<int, double, double>> * cells = nullptr)
{
    std::istringstream dfn(dfn_text);
    DefinitionParser parser;
    ASSERT_TRUE(parser.parse(dfn));
    auto * data = dynamic_cast<quadtree_builder_raw_data *>(parser.getBlockData("qtb"));
    ASSERT_NE(data, nullptr);
    ASSERT_TRUE(data->initialize(parser.getBlockData()));
    auto * qtree = dynamic_cast<QuadTree3D *>(data->build());
    ASSERT_NE(qtree, nullptr);
    if (nodes != nullptr) *nodes = qtree->nodes();
    if (cells != nullptr)
    {
        cells->clear();
        for (int n = 0; n < qtree->nodes(); n++)
        {
            Box * box = qtree->get_nodeobj(n);
            cells->emplace_back(box->depth, box->x, box->y);
        }
    }
}

} // namespace gridgen_test
