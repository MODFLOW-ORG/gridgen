
#ifndef SHPREADER_H
#define SHPREADER_H

#include "polygon.h"
#include "polyline.h"
#include "shapelib/shapefil.h"
#include <map>

using namespace std;

namespace cusg
{
struct ShapeIntField
{
    int id;
    string name;
    map<int, int> values;
};
struct ShapeDoubleField
{
    int id;
    string name;
    map<int, double> values;
};
struct ShapeStringField
{
    int id;
    string name;
    map<int, string> values;
};
struct ShapeLogField
{
    int id;
    string name;
    map<int, string> values;
};

// this reader reads shapefile format and convert
// polygons to poly format

class ShapeReader
{
public:
    ShapeReader()
    {
        m_already_read = false;
        m_bbox[0] = m_bbox[2] = DBL_MAX;
        m_bbox[1] = m_bbox[3] = -DBL_MAX;
    }

    ShapeReader(const ShapeReader &) = delete;
    ShapeReader & operator=(const ShapeReader &) = delete;

    ~ShapeReader() { destroy(); }

    //read field attributes
    bool read_attributes(DBFHandle & dbf_handle, int shape_index);

    //read from file
    bool read(const string & filename);

    //get poly/arc/pts
    vector<GisPolygon> & get_poly_list() { return m_polygon_list; }
    vector<GisPolyline> & get_arc_list() { return m_arc_list; }
    vector<GisPoint2d> & get_point_list() { return m_pt_list; }

    //set tags for getting information from database
    void set_base_elevation_tag(const string & tag) { m_base_elevation_tag = tag; }
    void set_height_tag(const string & tag) { m_height_tag = tag; }

    void destroy();

    vector<ShapeIntField> & get_int_fields() { return m_int_fields; }
    vector<ShapeDoubleField> & get_double_fields() { return m_double_fields; }
    vector<ShapeStringField> & get_str_fields() { return m_str_fields; }
    vector<ShapeLogField> & get_logical_fields() { return m_logical_fields; }

private:
    void add_polygon(Ring & ring, int shape_index);

    void read_polygon(SHPObject * shape, int shape_index);
    void read_arc(SHPObject * shape, int shape_index);
    void read_pts(SHPObject * shape, int shape_index);
    void read_multi_patch(SHPObject * shape, int shape_index);

    vector<GisPolygon> m_polygon_list; //a list of polygons
    vector<GisPolyline> m_arc_list;    //a list of arcs
    vector<GisPoint2d> m_pt_list;      //a list of points
    string m_base_elevation_tag;
    string m_height_tag;

    double m_bbox[4];

    bool m_already_read;
    string m_filename; //for warnings

    ////add by Guilin

    //field id, names and values
    vector<ShapeIntField> m_int_fields;
    vector<ShapeDoubleField> m_double_fields;
    vector<ShapeStringField> m_str_fields;
    vector<ShapeLogField> m_logical_fields;
};

} //namespace cusg

#endif //SHPREADER_H
