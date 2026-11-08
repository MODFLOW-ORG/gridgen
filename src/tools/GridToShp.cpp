#include "GridToShp.h"
#include "GridgenError.h"
#include "Box.h"
#include "GridVertices.h"
#include "shpGridToLine.h"
#include <fstream>
#include <sstream>
using namespace std;

//corners, clockwise from the upper left, and center of a box
inline void boxCorners(const GridVertices & gv, Box * b, double * corner_x, double * corner_y,
                       double * center_x, double * center_y)
{
    gv.get_corners(b, corner_x, corner_y);
    const Point2d & center = gv.center(b);
    center_x[0] = center[0];
    center_y[0] = center[1];
}

SHPHandle createShpHandle(string shpName, string feature_type, int & shptype)
{
    SHPHandle handler;
    if (feature_type == "polygon")
    {
        handler = SHPCreate(shpName.c_str(), SHPT_POLYGON);
        shptype = 1;
    }
    else if (feature_type == "point")
    {
        handler = SHPCreate(shpName.c_str(), SHPT_POINT);
        shptype = 0;
    }
    else
    {
        throw GridgenError("FEATURE_TYPE " + feature_type + " is not line, polygon, or point");
    }
    return handler;
}

void writeDB(int entityNum, DBFHandle & dbf, vector<int> & fieldIds, vector<int> & double_field_ids,
             int num, int layer, int row, int col, string & location, double top, double bot,
             double dx, double dy)
{
    //write to DB
    DBFWriteIntegerAttribute(dbf, entityNum, fieldIds[0], num);
    DBFWriteIntegerAttribute(dbf, entityNum, fieldIds[1], layer);
    DBFWriteIntegerAttribute(dbf, entityNum, fieldIds[2], row);
    DBFWriteIntegerAttribute(dbf, entityNum, fieldIds[3], col);

    DBFWriteStringAttribute(dbf, entityNum, fieldIds[4], location.c_str());

    DBFWriteDoubleAttribute(dbf, entityNum, double_field_ids[0], top);
    DBFWriteDoubleAttribute(dbf, entityNum, double_field_ids[1], bot);
    DBFWriteDoubleAttribute(dbf, entityNum, double_field_ids[2], dx);
    DBFWriteDoubleAttribute(dbf, entityNum, double_field_ids[3], dy);
}

//write feature and info using index
void writeFeatEntity(Grid * grid, SHPHandle & handler, DBFHandle & dbf, vector<int> & fieldIds,
                     vector<int> & double_field_ids, vector<Box *> & boxes, vector<int> & numbers,
                     vector<int> & layers, vector<int> & rows, vector<int> & cols,
                     vector<string> & locations, double cx, double cy, double rotation, int shptype,
                     bool one_base_number, vector<size_t> sort_idx, bool close_polygons)
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

        boxCorners(gv, boxes[id], corner_x, corner_y, center_x, center_y);

        SHPObject * shape = nullptr;
        int entityNum = -1;
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

            entityNum = SHPWriteObject(handler, -1, shape);
        }
        else if (shptype == 0)
        {
            shape = SHPCreateSimpleObject(SHPT_POINT, 1, center_x, center_y, nullptr);

            entityNum = SHPWriteObject(handler, -1, shape);
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

        writeDB(entityNum, dbf, fieldIds, double_field_ids, numbers[id], layer, row, col,
                locations[id], top, bot, boxes[id]->dx, boxes[id]->dy);
    }
}

////write feature info

void DFSCollect(Grid * grid, Box * b, int layer, int row, int col, string & location,
                bool one_based_numbering, vector<Box *> & boxes, vector<int> & numbers,
                vector<int> & layers, vector<int> & rows, vector<int> & cols,
                vector<string> & locations, bool without_inactive)
{

    char quadrant[8] = {'\0'};
    if (b->isLeaf)
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
    else if (b->isLeaf == false) //b is not a leaf, call this function recursively
    {
        for (int i = 0; i < 4; i++)
        {
            int quadrant_index = i;
            if (i == 2) quadrant_index = 3;
            if (i == 3) quadrant_index = 2;

            sprintf(quadrant, "%d", quadrant_index + 1);

            string child_location = location;
            child_location += quadrant;

            DFSCollect(grid, b->pChildren[i], layer, row, col, child_location, one_based_numbering,
                       boxes, numbers, layers, rows, cols, locations, without_inactive);
        }
    }
}

void write_shapefile(QuadTree3D * qtree, string shpName, string feature_type, bool close_polygons)
{
    if (feature_type == "line")
    {
        write_line_shapefile(qtree, shpName);
        return;
    }

    //create a shape file
    int shptype = true; //polygon : true ; point : false
    SHPHandle handler = createShpHandle(shpName, feature_type, shptype);
    if (handler == nullptr) throw GridgenError("cannot write shapefile " + shpName);

    //create a .dbf file
    string dbfName = shpName;

    if (shpName.find('.') != string::npos)
    {
        int dotIndex = shpName.find_last_of('.');
        dbfName = shpName.substr(0, dotIndex);
    }

    dbfName += ".dbf";
    DBFHandle dbfHandle = DBFCreate(dbfName.c_str());
    int nodeFieldId = DBFAddField(dbfHandle, "nodenumber", FTInteger, 12, 0);
    int layerFieldId = DBFAddField(dbfHandle, "layer", FTInteger, 20, 0);
    int rowFieldId = DBFAddField(dbfHandle, "row", FTInteger, 20, 0);
    int colFieldId = DBFAddField(dbfHandle, "col", FTInteger, 20, 0);
    int locationFieldId = DBFAddField(dbfHandle, "child_location", FTString, 50, 10);

    int topFieldId = DBFAddField(dbfHandle, "top", FTDouble, 20, 8);
    int bottomFieldId = DBFAddField(dbfHandle, "bottom", FTDouble, 20, 8);
    int delrFieldId = DBFAddField(dbfHandle, "delr", FTDouble, 20, 8);
    int delcFieldId = DBFAddField(dbfHandle, "delc", FTDouble, 20, 8);

    vector<int> fieldIds;
    fieldIds.push_back(nodeFieldId);
    fieldIds.push_back(layerFieldId);
    fieldIds.push_back(rowFieldId);
    fieldIds.push_back(colFieldId);
    fieldIds.push_back(locationFieldId);

    vector<int> double_field_ids;
    double_field_ids.push_back(topFieldId);
    double_field_ids.push_back(bottomFieldId);
    double_field_ids.push_back(delrFieldId);
    double_field_ids.push_back(delcFieldId);

    ////loop over all the leaf nodes

    ModflowGrid * mfgrid = qtree->getModflowGrid();

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
                Box * box = qtree->nodegroup[bid];
                string location;

                DFSCollect(qtree, box, l, r, c, location, qtree->get_one_based_node_numbering(),
                           boxes, numbers, layers, rows, cols, locations, true);
            }
        }
    }

    //sort
    std::vector<size_t> sort_idx;
    sort(numbers, sort_idx);

    //write the sorted
    writeFeatEntity(qtree, handler, dbfHandle, fieldIds, double_field_ids, boxes, numbers, layers,
                    rows, cols, locations, cx, cy, rotation, shptype,
                    qtree->get_one_based_node_numbering(), sort_idx, close_polygons);

    //close the shape file
    DBFClose(dbfHandle);
    SHPClose(handler);
}

void write_shapefile(ModflowGrid * modflow, string gridShpName, string feature_type,
                     bool without_inactive, bool one_based_numbering, bool close_polygons)
{
    if (feature_type == "line")
    {
        write_line_shapefile(modflow, gridShpName);
        return;
    }

    //create a shape file
    int shptype = 0; //polygon : true ; point : false
    SHPHandle handler = createShpHandle(gridShpName, feature_type, shptype);
    if (handler == nullptr) throw GridgenError("cannot write shapefile " + gridShpName);

    //create a .dbf file
    string dbfName = gridShpName; /*
	int lastL = gridShpName.find_last_of('/') == string::npos ? gridShpName.find;*/
    if (gridShpName.find('.') != string::npos)
    {
        int dotIndex = gridShpName.find_last_of('.');
        dbfName = gridShpName.substr(0, dotIndex);
    }

    dbfName += ".dbf";
    DBFHandle dbfHandle = DBFCreate(dbfName.c_str());
    int nodeFieldId = DBFAddField(dbfHandle, "nodenumber", FTInteger, 12, 0);
    int layerFieldId = DBFAddField(dbfHandle, "layer", FTInteger, 20, 0);
    int rowFieldId = DBFAddField(dbfHandle, "row", FTInteger, 20, 0);
    int colFieldId = DBFAddField(dbfHandle, "col", FTInteger, 20, 0);
    int locationFieldId = DBFAddField(dbfHandle, "child_location", FTString, 50, 10);

    int topFieldId = DBFAddField(dbfHandle, "top", FTDouble, 20, 8);
    int bottomFieldId = DBFAddField(dbfHandle, "bottom", FTDouble, 20, 8);
    int delrFieldId = DBFAddField(dbfHandle, "delr", FTDouble, 20, 8);
    int delcFieldId = DBFAddField(dbfHandle, "delc", FTDouble, 20, 8);

    vector<int> fieldIds;
    fieldIds.push_back(nodeFieldId);
    fieldIds.push_back(layerFieldId);
    fieldIds.push_back(rowFieldId);
    fieldIds.push_back(colFieldId);
    fieldIds.push_back(locationFieldId);

    vector<int> double_field_ids;
    double_field_ids.push_back(topFieldId);
    double_field_ids.push_back(bottomFieldId);
    double_field_ids.push_back(delrFieldId);
    double_field_ids.push_back(delcFieldId);

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
                Box * box = mfgrid->nodegroup[bid];
                string extra_str = extra;
                string location;

                DFSCollect(modflow, box, l, r, c, location, one_based_numbering, boxes, numbers,
                           layers, rows, cols, locations, true);
            }
        }
    }

    //sort
    std::vector<size_t> sort_idx;
    sort(numbers, sort_idx);

    //write the sorted
    writeFeatEntity(modflow, handler, dbfHandle, fieldIds, double_field_ids, boxes, numbers, layers,
                    rows, cols, locations, cx, cy, rotation, shptype, one_based_numbering, sort_idx,
                    close_polygons);

    //close the shape file
    DBFClose(dbfHandle);
    SHPClose(handler);
}
