#include "GridToShp.h"
#include "GridgenError.h"
#include "Box.h"
#include "GridVertices.h"
#include "shpGridToLine.h"
#include <fstream>
#include <sstream>
using namespace std;

//corners, clockwise from the upper left, and center of a box
inline void box_corners(const GridVertices & gv, Box * b, double * corner_x, double * corner_y,
                        double * center_x, double * center_y)
{
    gv.get_corners(b, corner_x, corner_y);
    const Point2d & center = gv.center(b);
    center_x[0] = center[0];
    center_y[0] = center[1];
}

SHPHandle create_shp_handle(string shp_name, string feature_type, int & shptype)
{
    SHPHandle handler;
    if (feature_type == "polygon")
    {
        handler = SHPCreate(shp_name.c_str(), SHPT_POLYGON);
        shptype = 1;
    }
    else if (feature_type == "point")
    {
        handler = SHPCreate(shp_name.c_str(), SHPT_POINT);
        shptype = 0;
    }
    else
    {
        throw GridgenError("FEATURE_TYPE " + feature_type + " is not line, polygon, or point");
    }
    return handler;
}

void write_db(int entity_num, DBFHandle & dbf, vector<int> & field_ids,
              vector<int> & double_field_ids, int num, int layer, int row, int col,
              string & location, double top, double bot, double dx, double dy)
{
    //write to DB
    DBFWriteIntegerAttribute(dbf, entity_num, field_ids[0], num);
    DBFWriteIntegerAttribute(dbf, entity_num, field_ids[1], layer);
    DBFWriteIntegerAttribute(dbf, entity_num, field_ids[2], row);
    DBFWriteIntegerAttribute(dbf, entity_num, field_ids[3], col);

    DBFWriteStringAttribute(dbf, entity_num, field_ids[4], location.c_str());

    DBFWriteDoubleAttribute(dbf, entity_num, double_field_ids[0], top);
    DBFWriteDoubleAttribute(dbf, entity_num, double_field_ids[1], bot);
    DBFWriteDoubleAttribute(dbf, entity_num, double_field_ids[2], dx);
    DBFWriteDoubleAttribute(dbf, entity_num, double_field_ids[3], dy);
}

//write feature and info using index
void write_feat_entity(Grid * grid, SHPHandle & handler, DBFHandle & dbf, vector<int> & field_ids,
                       vector<int> & double_field_ids, vector<Box *> & boxes, vector<int> & numbers,
                       vector<int> & layers, vector<int> & rows, vector<int> & cols,
                       vector<string> & locations, double cx, double cy, double rotation,
                       int shptype, bool one_base_number, vector<size_t> sort_idx,
                       bool close_polygons)
{
    QuadTree3D * qtree = dynamic_cast<QuadTree3D *>(grid);
    ModflowGrid * mfgrid = dynamic_cast<ModflowGrid *>(grid);
    bool is_quadtree = (qtree != nullptr);

    double corner_x[5] = {0}, corner_y[5] = {0}, center_x[1] = {0}, center_y[1] = {0};
    GridVertices gv(grid);

    int count = sort_idx.size(), id;
    for (int i = 0; i < count; i++)
    {
        id = sort_idx[i];

        box_corners(gv, boxes[id], corner_x, corner_y, center_x, center_y);

        SHPObject * shape = nullptr;
        int entity_num = -1;
        if (shptype == 1)
        {
            if (!close_polygons)
                shape = SHPCreateSimpleObject(SHPT_POLYGON, 4, corner_x, corner_y, nullptr);
            else
            {
                corner_x[4] = corner_x[0];
                corner_y[4] = corner_y[0];
                shape = SHPCreateSimpleObject(SHPT_POLYGON, 5, corner_x, corner_y, nullptr);
            }

            entity_num = SHPWriteObject(handler, -1, shape);
        }
        else if (shptype == 0)
        {
            shape = SHPCreateSimpleObject(SHPT_POINT, 1, center_x, center_y, nullptr);

            entity_num = SHPWriteObject(handler, -1, shape);
        }
        else
            assert(false);
        SHPDestroyObject(shape);

        int layer = (one_base_number ? (layers[id] - 1) : layers[id]);
        int row = (one_base_number ? (rows[id] - 1) : rows[id]);
        int col = (one_base_number ? (cols[id] - 1) : cols[id]);

        double top, bot;
        if (is_quadtree)
        {
            top = qtree->top[boxes[id]->number];
            bot = qtree->bot[boxes[id]->number];
        }
        else
        {
            top = mfgrid->bot[layer * mfgrid->nrow * mfgrid->ncol + row * mfgrid->ncol + col];
            bot = mfgrid->bot[(layer + 1) * mfgrid->nrow * mfgrid->ncol + row * mfgrid->ncol + col];
        }

        write_db(entity_num, dbf, field_ids, double_field_ids, numbers[id], layer, row, col,
                 locations[id], top, bot, boxes[id]->dx, boxes[id]->dy);
    }
}

////write feature info

void dfs_collect(Grid * grid, Box * b, int layer, int row, int col, string & location,
                 bool one_based_numbering, vector<Box *> & boxes, vector<int> & numbers,
                 vector<int> & layers, vector<int> & rows, vector<int> & cols,
                 vector<string> & locations, bool without_inactive)
{

    char quadrant[8] = {'\0'};
    if (b->is_leaf)
    {
        if (without_inactive && (!b->active)) return;

        boxes.push_back(b);
        numbers.push_back(b->number == -1 ? b->number
                                          : (one_based_numbering ? b->number + 1 : b->number));
        layers.push_back(one_based_numbering ? (layer + 1) : layer);
        rows.push_back(one_based_numbering ? (row + 1) : row);
        cols.push_back(one_based_numbering ? (col + 1) : col);
        locations.push_back(location);
    }
    else if (b->is_leaf == false) //b is not a leaf, call this function recursively
    {
        for (int i = 0; i < 4; i++)
        {
            int quadrant_index = i;
            if (i == 2) quadrant_index = 3;
            if (i == 3) quadrant_index = 2;

            sprintf(quadrant, "%d", quadrant_index + 1);

            string child_location = location;
            child_location += quadrant;

            dfs_collect(grid, b->p_children[i], layer, row, col, child_location,
                        one_based_numbering, boxes, numbers, layers, rows, cols, locations,
                        without_inactive);
        }
    }
}

void write_shapefile(QuadTree3D * qtree, string shp_name, string feature_type, bool close_polygons)
{
    if (feature_type == "line")
    {
        write_line_shapefile(qtree, shp_name);
        return;
    }

    //create a shape file
    int shptype = true; //polygon : true ; point : false
    SHPHandle handler = create_shp_handle(shp_name, feature_type, shptype);
    if (handler == nullptr) throw GridgenError("cannot write shapefile " + shp_name);

    //create a .dbf file
    string dbf_name = shp_name;

    if (shp_name.find('.') != string::npos)
    {
        int dot_index = shp_name.find_last_of('.');
        dbf_name = shp_name.substr(0, dot_index);
    }

    dbf_name += ".dbf";
    DBFHandle dbf_handle = DBFCreate(dbf_name.c_str());
    int node_field_id = DBFAddField(dbf_handle, "nodenumber", FTInteger, 12, 0);
    int layer_field_id = DBFAddField(dbf_handle, "layer", FTInteger, 20, 0);
    int row_field_id = DBFAddField(dbf_handle, "row", FTInteger, 20, 0);
    int col_field_id = DBFAddField(dbf_handle, "col", FTInteger, 20, 0);
    int location_field_id = DBFAddField(dbf_handle, "child_location", FTString, 50, 10);

    int top_field_id = DBFAddField(dbf_handle, "top", FTDouble, 20, 8);
    int bottom_field_id = DBFAddField(dbf_handle, "bottom", FTDouble, 20, 8);
    int delr_field_id = DBFAddField(dbf_handle, "delr", FTDouble, 20, 8);
    int delc_field_id = DBFAddField(dbf_handle, "delc", FTDouble, 20, 8);

    vector<int> field_ids;
    field_ids.push_back(node_field_id);
    field_ids.push_back(layer_field_id);
    field_ids.push_back(row_field_id);
    field_ids.push_back(col_field_id);
    field_ids.push_back(location_field_id);

    vector<int> double_field_ids;
    double_field_ids.push_back(top_field_id);
    double_field_ids.push_back(bottom_field_id);
    double_field_ids.push_back(delr_field_id);
    double_field_ids.push_back(delc_field_id);

    ////loop over all the leaf nodes

    ModflowGrid * mfgrid = qtree->get_modflow_grid();

    double cx, cy, rotation;
    mfgrid->get_rotation(cx, cy, rotation);

    vector<Box *> boxes;
    vector<int> numbers, layers, rows, cols;
    vector<string> locations;
    //write shapefile data for each base grid cell
    for (int l = 0; l < mfgrid->nlay; l++)
    {
        for (int r = 0; r < mfgrid->nrow; r++)
        {
            for (int c = 0; c < mfgrid->ncol; c++)
            {
                int bid = mfgrid->get_box_index(l, r, c);
                Box * box = qtree->node_group[bid];
                string location;

                dfs_collect(qtree, box, l, r, c, location, qtree->get_one_based_node_numbering(),
                            boxes, numbers, layers, rows, cols, locations, true);
            }
        }
    }

    //sort
    std::vector<size_t> sort_idx;
    sort(numbers, sort_idx);

    //write the sorted
    write_feat_entity(qtree, handler, dbf_handle, field_ids, double_field_ids, boxes, numbers,
                      layers, rows, cols, locations, cx, cy, rotation, shptype,
                      qtree->get_one_based_node_numbering(), sort_idx, close_polygons);

    //close the shape file
    DBFClose(dbf_handle);
    SHPClose(handler);
}

void write_shapefile(ModflowGrid * modflow, string grid_shp_name, string feature_type,
                     bool without_inactive, bool one_based_numbering, bool close_polygons)
{
    if (feature_type == "line")
    {
        write_line_shapefile(modflow, grid_shp_name);
        return;
    }

    //create a shape file
    int shptype = 0; //polygon : true ; point : false
    SHPHandle handler = create_shp_handle(grid_shp_name, feature_type, shptype);
    if (handler == nullptr) throw GridgenError("cannot write shapefile " + grid_shp_name);

    //create a .dbf file
    string dbf_name = grid_shp_name; /*
	int lastL = grid_shp_name.find_last_of('/') == string::npos ? grid_shp_name.find;*/
    if (grid_shp_name.find('.') != string::npos)
    {
        int dot_index = grid_shp_name.find_last_of('.');
        dbf_name = grid_shp_name.substr(0, dot_index);
    }

    dbf_name += ".dbf";
    DBFHandle dbf_handle = DBFCreate(dbf_name.c_str());
    int node_field_id = DBFAddField(dbf_handle, "nodenumber", FTInteger, 12, 0);
    int layer_field_id = DBFAddField(dbf_handle, "layer", FTInteger, 20, 0);
    int row_field_id = DBFAddField(dbf_handle, "row", FTInteger, 20, 0);
    int col_field_id = DBFAddField(dbf_handle, "col", FTInteger, 20, 0);
    int location_field_id = DBFAddField(dbf_handle, "child_location", FTString, 50, 10);

    int top_field_id = DBFAddField(dbf_handle, "top", FTDouble, 20, 8);
    int bottom_field_id = DBFAddField(dbf_handle, "bottom", FTDouble, 20, 8);
    int delr_field_id = DBFAddField(dbf_handle, "delr", FTDouble, 20, 8);
    int delc_field_id = DBFAddField(dbf_handle, "delc", FTDouble, 20, 8);

    vector<int> field_ids;
    field_ids.push_back(node_field_id);
    field_ids.push_back(layer_field_id);
    field_ids.push_back(row_field_id);
    field_ids.push_back(col_field_id);
    field_ids.push_back(location_field_id);

    vector<int> double_field_ids;
    double_field_ids.push_back(top_field_id);
    double_field_ids.push_back(bottom_field_id);
    double_field_ids.push_back(delr_field_id);
    double_field_ids.push_back(delc_field_id);

    ////loop over all the leaf nodes

    ModflowGrid * mfgrid = modflow;
    char extra[64];

    double cx, cy, rotation;
    modflow->get_rotation(cx, cy, rotation);

    vector<Box *> boxes;
    vector<int> numbers, layers, rows, cols;
    vector<string> locations;
    for (int l = 0; l < mfgrid->nlay; l++)
    {
        for (int r = 0; r < mfgrid->nrow; r++)
        {
            for (int c = 0; c < mfgrid->ncol; c++)
            {
                int bid = mfgrid->get_box_index(l, r, c);
                Box * box = mfgrid->node_group[bid];
                string extra_str = extra;
                string location;

                dfs_collect(modflow, box, l, r, c, location, one_based_numbering, boxes, numbers,
                            layers, rows, cols, locations, true);
            }
        }
    }

    //sort
    std::vector<size_t> sort_idx;
    sort(numbers, sort_idx);

    //write the sorted
    write_feat_entity(modflow, handler, dbf_handle, field_ids, double_field_ids, boxes, numbers,
                      layers, rows, cols, locations, cx, cy, rotation, shptype, one_based_numbering,
                      sort_idx, close_polygons);

    //close the shape file
    DBFClose(dbf_handle);
    SHPClose(handler);
}
