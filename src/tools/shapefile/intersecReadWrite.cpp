#include "intersecReadWrite.h"
#include "RoundTrip.h"
#include "Warnings.h"
#include <sstream>
using namespace std;
using namespace cusg;

const string INTERSECTION_DELIMITER = ",";

// signed area of a ring, positive when it runs counterclockwise
static double ring_area(LinkedVertex * head)
{
    double area = 0.0;
    LinkedVertex * pv = head;
    do
    {
        const Point2d & a = pv->get_pos();
        const Point2d & b = pv->get_next()->get_pos();
        area += a[0] * b[1] - b[0] * a[1];
        pv = pv->get_next();
    } while (pv != head);
    return 0.5 * area;
}

bool check_attributes_contain(const char * attribute_name, vector<string> & attributes)
{
    for (auto & str : attributes)
    {
        if (str == attribute_name || str == "all")
            return true;
        else if (str == "none")
            return false;
    }
    return false;
}
void check_shp_attr_contain(GridIntersection * grid_intersection, vector<bool> & write_int_field,
                            vector<bool> & write_double_field, vector<bool> & write_str_field,
                            vector<bool> & write_logical_field, vector<string> & int_names,
                            vector<string> & double_names, vector<string> & str_names,
                            vector<string> & logical_names, vector<string> & attributes)
{
    vector<ShapeIntField> & int_fields = grid_intersection->m_shp_reader.get_int_fields();
    vector<ShapeDoubleField> & double_fields = grid_intersection->m_shp_reader.get_double_fields();
    vector<ShapeStringField> & str_fields = grid_intersection->m_shp_reader.get_str_fields();
    vector<ShapeLogField> & logical_fields = grid_intersection->m_shp_reader.get_logical_fields();
    for (auto & fld : int_fields)
    {
        if (check_attributes_contain(fld.name.c_str(), attributes))
        {
            write_int_field.push_back(true);
            int_names.push_back(fld.name);
        }
        else
            write_int_field.push_back(false);
    }
    for (auto & fld : double_fields)
    {
        if (check_attributes_contain(fld.name.c_str(), attributes))
        {
            write_double_field.push_back(true);
            double_names.push_back(fld.name);
        }
        else
            write_double_field.push_back(false);
    }
    for (auto & fld : str_fields)
    {
        if (check_attributes_contain(fld.name.c_str(), attributes))
        {
            write_str_field.push_back(true);
            str_names.push_back(fld.name);
        }
        else
            write_str_field.push_back(false);
    }
    for (auto & fld : logical_fields)
    {
        if (check_attributes_contain(fld.name.c_str(), attributes))
        {
            write_logical_field.push_back(true);
            logical_names.push_back(fld.name);
        }
        else
            write_logical_field.push_back(false);
    }
}
void add_shp_attr_val(int shape_index, GridIntersection * grid_intersection,
                      vector<bool> & write_int_field, vector<bool> & write_double_field,
                      vector<bool> & write_str_field, vector<bool> & write_logical_field,
                      vector<int> & int_values, vector<double> & double_values,
                      vector<string> & string_values, vector<string> & logical_values)
{
    vector<ShapeIntField> & int_fields = grid_intersection->m_shp_reader.get_int_fields();
    vector<ShapeDoubleField> & double_fields = grid_intersection->m_shp_reader.get_double_fields();
    vector<ShapeStringField> & str_fields = grid_intersection->m_shp_reader.get_str_fields();
    vector<ShapeLogField> & logical_fields = grid_intersection->m_shp_reader.get_logical_fields();
    vector<bool>::iterator ibit = write_int_field.begin();
    for (vector<ShapeIntField>::iterator sit = int_fields.begin(); sit != int_fields.end();
         ++sit, ++ibit)
    {
        if (*ibit) int_values.push_back(sit->values[shape_index]);
    }
    vector<bool>::iterator idit = write_double_field.begin();
    for (vector<ShapeDoubleField>::iterator sit = double_fields.begin(); sit != double_fields.end();
         ++sit, ++idit)
    {
        if (*idit) double_values.push_back(sit->values[shape_index]);
    }
    vector<bool>::iterator sbit = write_str_field.begin();
    for (vector<ShapeStringField>::iterator sit = str_fields.begin(); sit != str_fields.end();
         ++sit, ++sbit)
    {
        if (*sbit) string_values.push_back(sit->values[shape_index]);
    }
    vector<bool>::iterator lbit = write_logical_field.begin();
    for (vector<ShapeLogField>::iterator sit = logical_fields.begin(); sit != logical_fields.end();
         ++sit, ++lbit)
    {
        if (*lbit) logical_values.push_back(sit->values[shape_index]);
    }
}

void add_db_handle_fields(DBFHandle & dbf_handle, int int_field_count,
                          vector<string> & int_field_names, int double_field_count,
                          vector<string> & double_field_names, vector<int> & int_field_ids,
                          vector<int> & double_field_ids)
{
    for (int i = 0; i < int_field_count; i++)
    {
        int tid = DBFAddField(dbf_handle, int_field_names[i].c_str(), FTInteger, 20, 0);
        int_field_ids.push_back(tid);
    }
    for (int i = 0; i < double_field_count; i++)
    {
        int tid = DBFAddField(dbf_handle, double_field_names[i].c_str(), FTDouble, 20, 8);
        double_field_ids.push_back(tid);
    }
}
void add_db_handle_str_fields(DBFHandle & dbf_handle, int str_field_count,
                              vector<string> & str_field_names, vector<int> & str_field_ids)
{
    for (int i = 0; i < str_field_count; i++)
    {
        int tid = DBFAddField(dbf_handle, str_field_names[i].c_str(), FTString, 50, 0);
        str_field_ids.push_back(tid);
    }
}
//write int and double attributes
void write_int_db_val(DBFHandle & dbf_handle, int entity, vector<int> & int_field_ids,
                      vector<int> & int_values, vector<int> & double_field_ids,
                      vector<double> & double_values, vector<int> & str_field_ids,
                      vector<string> & str_values, vector<int> & logical_field_ids,
                      vector<string> & logical_values, ofstream & ofile, bool is_txt)
{
    int int_field_count = int_values.size();
    int double_field_count = double_values.size();
    int str_field_count = str_values.size();
    int logical_field_count = logical_values.size();
    //write integer fields
    for (int i = 0; i < int_field_count; ++i)
    {
        if (!is_txt)
            DBFWriteIntegerAttribute(dbf_handle, entity, int_field_ids[i], int_values[i]);
        else
            ofile << int_values[i] << INTERSECTION_DELIMITER;
    }
    for (int i = 0; i < double_field_count; ++i)
    {
        if (!is_txt)
            DBFWriteDoubleAttribute(dbf_handle, entity, double_field_ids[i], double_values[i]);
        else
            ofile << RoundTrip(double_values[i]) << INTERSECTION_DELIMITER;
    }
    for (int i = 0; i < str_field_count; i++)
    {
        if (!is_txt)
            DBFWriteStringAttribute(dbf_handle, entity, str_field_ids[i], str_values[i].c_str());
        else
            ofile << str_values[i].c_str() << INTERSECTION_DELIMITER;
    }
    for (int i = 0; i < logical_field_count; i++)
    {
        if (!is_txt)
            DBFWriteStringAttribute(dbf_handle, entity, logical_field_ids[i],
                                    logical_values[i].c_str());
        else
            ofile << logical_values[i].c_str() << INTERSECTION_DELIMITER;
    }
    ofile << endl;
}
void write_intersection(string file, int int_field_count, vector<vector<int>> & int_field_values,
                        vector<string> & int_field_names, int double_field_count,
                        vector<vector<double>> & double_field_values,
                        vector<string> & double_field_names, int str_field_count,
                        vector<vector<string>> & str_field_values, vector<string> & str_field_names,
                        int logical_field_count, vector<vector<string>> & logical_field_values,
                        vector<string> & logical_field_names, GridIntersection * grid_intersection,
                        WRITETYPE wtype, bool is_txt)
{
    //create a shape file
    ofstream ofile;
    SHPHandle handler = nullptr;
    if (wtype == WRITE_POINT)
    {
        if (!is_txt)
            handler = SHPCreate(file.c_str(), SHPT_POINT);
        else
        {
            ofile.open(file.c_str(), std::ofstream::out | std::ofstream::trunc);
        }
    }
    else if (wtype == WRITE_POLYLINE)
    {
        if (!is_txt)
            handler = SHPCreate(file.c_str(), SHPT_ARC);
        else
        {
            ofile.open(file.c_str(), std::ofstream::out | std::ofstream::trunc);
        }
    }
    else
    {
        if (!is_txt)
            handler = SHPCreate(file.c_str(), SHPT_POLYGON);
        else
        {
            ofile.open(file.c_str(), std::ofstream::out | std::ofstream::trunc);
        }
    }

    //create a .dbf file
    DBFHandle dbf_handle = nullptr;
    if (!is_txt)
    {
        string dbf = file;
        if (file.find('.') != string::npos)
        {
            int dot_index = file.find_last_of('.');
            dbf = file.substr(0, dot_index);
        }
        dbf += ".dbf";
        dbf_handle = DBFCreate(dbf.c_str());
    }

    //create cell field
    int node_field_id = -1;
    if (!is_txt) node_field_id = DBFAddField(dbf_handle, "nodenumber", FTInteger, 12, 0);
    vector<int> int_field_ids;     //create integer fields
    vector<int> double_field_ids;  //create double fields
    vector<int> str_field_ids;     //create string fields
    vector<int> logical_field_ids; //log fields
    if (!is_txt)
    {
        add_db_handle_fields(dbf_handle, int_field_count, int_field_names, double_field_count,
                             double_field_names, int_field_ids, double_field_ids);
        add_db_handle_str_fields(dbf_handle, str_field_count, str_field_names, str_field_ids);
        add_db_handle_str_fields(dbf_handle, logical_field_count, logical_field_names,
                                 logical_field_ids);
    }
    else
    {
        ofile << "nodenumber" << INTERSECTION_DELIMITER;
        for (int i = 0; i < int_field_count; ++i)
        {
            ofile << int_field_names[i] << INTERSECTION_DELIMITER;
        }
        for (int i = 0; i < double_field_count; ++i)
        {
            ofile << double_field_names[i] << INTERSECTION_DELIMITER;
        }
        for (int i = 0; i < str_field_count; i++)
        {
            ofile << str_field_names[i].c_str() << INTERSECTION_DELIMITER;
        }
        for (int i = 0; i < logical_field_count; i++)
        {
            ofile << logical_field_names[i].c_str() << INTERSECTION_DELIMITER;
        }
        ofile << endl;
    }

    vector<vector<int>>::iterator iit = int_field_values.begin();
    vector<vector<double>>::iterator dit = double_field_values.begin();
    vector<vector<string>>::iterator sit = str_field_values.begin();
    vector<vector<string>>::iterator lit = logical_field_values.begin();
    vector<double> xs, ys; //vertex coordinates of the shape being written
    double *x = nullptr, *y = nullptr;
    //write points
    if (wtype == WRITE_POINT)
    {
        xs.assign(1, 0.0);
        ys.assign(1, 0.0);
        x = xs.data();
        y = ys.data();
        vector<PointCellIntersectionInfo> & point_intersections =
            grid_intersection->get_point_intersection_info();
        vector<PointCellIntersectionInfo>::iterator pit = point_intersections.begin();
        for (; pit != point_intersections.end(); ++pit, ++iit, ++dit, ++sit, ++lit)
        {
            vector<int> & int_values = *iit;
            vector<double> & double_values = *dit;
            vector<string> & str_values = *sit;
            vector<string> & logical_values = *lit;
            x[0] = pit->point[0];
            y[0] = pit->point[1];

            SHPObject * shp_obj = nullptr;
            int entity = -1;
            if (!is_txt)
            {
                shp_obj = SHPCreateSimpleObject(SHPT_POINT, 1, x, y, nullptr);
                entity = SHPWriteObject(handler, -1, shp_obj);
            }

            //write cell
            if (!is_txt)
                DBFWriteIntegerAttribute(dbf_handle, entity, node_field_id,
                                         pit->box->number == -1 ? -1 : pit->box->number + 1);
            else
                ofile << (pit->box->number == -1 ? -1 : pit->box->number + 1)
                      << INTERSECTION_DELIMITER;
            write_int_db_val(dbf_handle, entity, int_field_ids, int_values, double_field_ids,
                             double_values, str_field_ids, str_values, logical_field_ids,
                             logical_values, ofile, is_txt);
            if (!is_txt) SHPDestroyObject(shp_obj);
        }
    }
    else if (wtype == WRITE_POLYLINE)
    {
        vector<PolylineCellIntersectionInfo> & polyline_intersections =
            grid_intersection->get_polyline_intersection_info();
        vector<PolylineCellIntersectionInfo>::iterator pit = polyline_intersections.begin();
        for (; pit != polyline_intersections.end(); ++pit, ++iit, ++dit, ++sit, ++lit)
        {
            PolylineCellIntersectionInfo & line_info = *pit;
            vector<int> & int_values = *iit;
            vector<double> & double_values = *dit;
            vector<string> & str_values = *sit;
            vector<string> & logical_values = *lit;

            int vertex_count = line_info.m_polyline->get_size();
            if (vertex_count == 1)
            {
                warn("The part of a line feature in node " + to_string(line_info.box->number + 1) +
                     " is a single point and is not written to " + file + ".");
                continue;
            }
            xs.assign(vertex_count, 0.0);
            ys.assign(vertex_count, 0.0);
            x = xs.data();
            y = ys.data();
            LinkedVertex * pv = line_info.m_polyline->get_head();
            x[0] = pv->get_pos()[0];
            y[0] = pv->get_pos()[1];
            int vertex_index = 0;
            do
            {
                pv = pv->get_next();
                vertex_index++;
                x[vertex_index] = pv->get_pos()[0];
                y[vertex_index] = pv->get_pos()[1];
            } while (pv != line_info.m_polyline->get_tail());

            SHPObject * shp_obj = nullptr;
            int entity = -1;
            if (!is_txt)
            {
                shp_obj = SHPCreateSimpleObject(SHPT_ARC, vertex_count, x, y, nullptr);
                entity = SHPWriteObject(handler, -1, shp_obj);
            }

            //write cell
            if (!is_txt)
                DBFWriteIntegerAttribute(dbf_handle, entity, node_field_id,
                                         pit->box->number == -1 ? -1 : pit->box->number + 1);
            else
                ofile << (pit->box->number == -1 ? -1 : pit->box->number + 1)
                      << INTERSECTION_DELIMITER;
            write_int_db_val(dbf_handle, entity, int_field_ids, int_values, double_field_ids,
                             double_values, str_field_ids, str_values, logical_field_ids,
                             logical_values, ofile, is_txt);

            if (!is_txt) SHPDestroyObject(shp_obj);
        }
    }
    else if (wtype == WRITE_POLYGON)
    {
        vector<PolygonCellIntersectionInfo> & polygon_intersections =
            grid_intersection->get_polygon_intersection_info();
        vector<PolygonCellIntersectionInfo>::iterator pit = polygon_intersections.begin();
        for (int shape_index = 0; pit != polygon_intersections.end();
             ++pit, ++iit, ++dit, ++shape_index, ++sit, ++lit)
        {
            PolygonCellIntersectionInfo & polygon_info = *pit;
            vector<int> & int_values = *iit;
            vector<double> & double_values = *dit;
            vector<string> & str_values = *sit;
            vector<string> & logical_values = *lit;

            int part_size = polygon_info.m_inter_polygon->size();
            int vsize = polygon_info.m_inter_polygon->get_size();
            if (vsize < 3)
            {
                warn("The part of a polygon feature in node " +
                     to_string(polygon_info.box->number + 1) +
                     " has fewer than three vertices and is not written to " + file + ".");
                continue;
            }
            // each ring repeats its first vertex at the end
            vsize += part_size;
            xs.assign(vsize, 0.0);
            ys.assign(vsize, 0.0);
            x = xs.data();
            y = ys.data();

            //make sure the out boundary
            list<Polygon::iterator> ordered_p_itrs;
            for (Polygon::iterator ring_it = polygon_info.m_inter_polygon->begin();
                 ring_it != polygon_info.m_inter_polygon->end(); ++ring_it)
            {
                if (ring_it->get_type() == Ring::POUT)
                    ordered_p_itrs.push_front(ring_it);
                else
                    ordered_p_itrs.push_back(ring_it);
            }

            //part number
            int part_start = 0;
            vector<int> parts(part_size, 0);
            int * part_starts = parts.data();
            int k = 0, vid = 0;
            for (list<Polygon::iterator>::iterator lit = ordered_p_itrs.begin();
                 lit != ordered_p_itrs.end(); ++lit, k++)
            {
                Polygon::iterator & ring_it = *lit;
                part_starts[k] = part_start;
                part_start += ring_it->get_size() + 1;

                // the shapefile specification orders the outer ring clockwise
                // and holes counterclockwise
                LinkedVertex * head = ring_it->get_head();
                bool forward = (ring_area(head) < 0.0) == (ring_it->get_type() == Ring::POUT);
                LinkedVertex * pv = head;
                do
                {
                    x[vid] = pv->get_pos()[0];
                    y[vid] = pv->get_pos()[1];
                    pv = forward ? pv->get_next() : pv->get_pre();
                    vid++;
                } while (pv != head);
                x[vid] = head->get_pos()[0];
                y[vid] = head->get_pos()[1];
                vid++;
            }

            //create a polygon shape here
            //SHPObject * shp_obj = SHPCreateObject(SHPT_POLYGON,vsize,x, y,NULL);
            SHPObject * shp_obj = nullptr;
            int entity = 0;
            if (!is_txt)
            {
                shp_obj = SHPCreateObject(SHPT_POLYGON, shape_index, part_size, part_starts,
                                          nullptr, vsize, x, y, nullptr, nullptr);
                entity = SHPWriteObject(handler, -1, shp_obj);
            }
            //write cell
            if (!is_txt)
                DBFWriteIntegerAttribute(dbf_handle, entity, node_field_id,
                                         pit->box->number == -1 ? -1 : pit->box->number + 1);
            else
                ofile << (pit->box->number == -1 ? -1 : pit->box->number + 1)
                      << INTERSECTION_DELIMITER;
            write_int_db_val(dbf_handle, entity, int_field_ids, int_values, double_field_ids,
                             double_values, str_field_ids, str_values, logical_field_ids,
                             logical_values, ofile, is_txt);

            if (!is_txt) SHPDestroyObject(shp_obj);
        }
    }

    if (!is_txt)
    {
        DBFClose(dbf_handle);
        SHPClose(handler);
    }
    else
        ofile.close();
}

//write the grid and point intersection
void write_pnt_intersect(GridIntersection * grid_intersection, string file,
                         vector<string> & attributes, FILETYPE filetype)
{
    vector<PointCellIntersectionInfo> & point_intersections =
        grid_intersection->get_point_intersection_info();

    vector<string> int_field_names, double_field_names, str_field_names, logical_field_names;
    vector<vector<int>> intflds;
    vector<vector<double>> dbflds;
    vector<vector<string>> strflds;
    vector<vector<string>> logflds;
    //to record whether each attribute needs to be written
    vector<bool> write_int_field, write_double_field, write_str_field, write_logical_field;

    //parse the attributes
    bool write_point_id =
        true; // check_attributes_contain("pointid", attributes);//write shapefile_pointid
    if (write_point_id) int_field_names.push_back(string("pointid"));
    //add the attributes in shape file
    check_shp_attr_contain(grid_intersection, write_int_field, write_double_field, write_str_field,
                           write_logical_field, int_field_names, double_field_names,
                           str_field_names, logical_field_names, attributes);

    for (auto & info : point_intersections)
    {
        vector<int> int_values;
        vector<double> dbs;
        vector<string> strs;
        vector<string> logs;
        if (write_point_id) int_values.push_back(info.m_point_id);

        add_shp_attr_val(info.m_point_id, grid_intersection, write_int_field, write_double_field,
                         write_str_field, write_logical_field, int_values, dbs, strs, logs);

        intflds.push_back(int_values);
        dbflds.push_back(dbs);
        strflds.push_back(strs);
        logflds.push_back(logs);
    }

    if (filetype == WRITE_SHP || filetype == WRITE_TXT)
    {
        bool is_txt = (filetype == WRITE_TXT);
        write_intersection(file, int_field_names.size(), intflds, int_field_names,
                           double_field_names.size(), dbflds, double_field_names,
                           str_field_names.size(), strflds, str_field_names,
                           logical_field_names.size(), logflds, logical_field_names,
                           grid_intersection, WRITE_POINT, is_txt);
    }
    else
    {
    }
}

//write the line and grid intersection
void write_line_intersect(GridIntersection * grid_intersection, string file,
                          vector<string> & attributes, FILETYPE filetype)
{
    vector<PolylineCellIntersectionInfo> & polyline_infos =
        grid_intersection->get_polyline_intersection_info();

    vector<string> int_field_names, double_field_names, str_field_names, logical_field_names;
    vector<vector<int>> intflds;
    vector<vector<double>> dbflds;
    vector<vector<string>> strflds;
    vector<vector<string>> logflds;
    //to record whether each attribute needs to be written
    vector<bool> write_int_field, write_double_field, write_str_field, write_logical_field;

    bool write_arc_id =
        true; // check_attributes_contain("arcid", attributes);//write shapefile_arcid
    bool write_length = true; //check_attributes_contain("length", attributes);//write length
    bool write_start_distance =
        true; //check_attributes_contain("starting_distance", attributes);//write starting_distance
    bool write_end_distance =
        true; //check_attributes_contain("ending_distance", attributes);//write ending_distance

    if (write_arc_id) int_field_names.push_back(string("arcid"));
    if (write_length) double_field_names.push_back(string("length"));
    if (write_start_distance) double_field_names.push_back(string("starting_distance"));
    if (write_end_distance) double_field_names.push_back(string("ending_distance"));
    //add the shape file attributes
    check_shp_attr_contain(grid_intersection, write_int_field, write_double_field, write_str_field,
                           write_logical_field, int_field_names, double_field_names,
                           str_field_names, logical_field_names, attributes);

    for (auto & info : polyline_infos)
    {
        vector<int> int_values;
        vector<double> dbs;
        vector<string> strs;
        vector<string> logs;

        if (write_arc_id) int_values.push_back(info.m_curve_id);
        if (write_length) dbs.push_back(info.m_len);
        if (write_start_distance) dbs.push_back(info.m_start_geo_dist);
        if (write_end_distance) dbs.push_back(info.m_end_geo_dist);

        add_shp_attr_val(info.m_curve_id, grid_intersection, write_int_field, write_double_field,
                         write_str_field, write_logical_field, int_values, dbs, strs, logs);

        intflds.push_back(int_values);
        dbflds.push_back(dbs);
        strflds.push_back(strs);
        logflds.push_back(logs);
    }

    if (filetype == WRITE_SHP || filetype == WRITE_TXT)
    {
        bool is_txt = (filetype == WRITE_TXT);
        write_intersection(file, int_field_names.size(), intflds, int_field_names,
                           double_field_names.size(), dbflds, double_field_names,
                           str_field_names.size(), strflds, str_field_names,
                           logical_field_names.size(), logflds, logical_field_names,
                           grid_intersection, WRITE_POLYLINE, is_txt);
    }
    else
    {
    }
}

//writhe the polygon and grid intersection
void write_poly_intersect(GridIntersection * grid_intersection, string file,
                          vector<string> & attributes, FILETYPE filetype)
{
    vector<PolygonCellIntersectionInfo> & sci = grid_intersection->get_polygon_intersection_info();

    vector<string> int_field_names, double_field_names, str_field_names, logical_field_names;
    vector<vector<int>> intflds;
    vector<vector<double>> dbflds;
    vector<vector<string>> strflds;
    vector<vector<string>> logflds;
    //to record whether each attribute needs to be written
    vector<bool> write_int_field, write_double_field, write_str_field, write_logical_field;

    bool write_polygon_id =
        true; // check_attributes_contain("polyid", attributes);//write shapefile_polyid
    bool write_area = true; //check_attributes_contain("totalarea", attributes);//write totalarea
    if (write_polygon_id) int_field_names.push_back(string("polyid"));
    if (write_area) double_field_names.push_back(string("totalarea"));

    check_shp_attr_contain(grid_intersection, write_int_field, write_double_field, write_str_field,
                           write_logical_field, int_field_names, double_field_names,
                           str_field_names, logical_field_names, attributes);

    for (auto & info : sci)
    {
        vector<int> int_values;
        vector<double> dbs;
        vector<string> strs;
        vector<string> logs;

        if (write_polygon_id) int_values.push_back(info.m_poly_id);
        if (write_area) dbs.push_back(info.m_area);

        add_shp_attr_val(info.m_poly_id, grid_intersection, write_int_field, write_double_field,
                         write_str_field, write_logical_field, int_values, dbs, strs, logs);

        intflds.push_back(int_values);
        dbflds.push_back(dbs);
        strflds.push_back(strs);
        logflds.push_back(logs);
    }

    if (filetype == WRITE_SHP || filetype == WRITE_TXT)
    {
        bool is_txt = (filetype == WRITE_TXT);
        write_intersection(file, int_field_names.size(), intflds, int_field_names,
                           double_field_names.size(), dbflds, double_field_names,
                           str_field_names.size(), strflds, str_field_names,
                           logical_field_names.size(), logflds, logical_field_names,
                           grid_intersection, WRITE_POLYGON, is_txt);
    }
    else
    {
    }
}
