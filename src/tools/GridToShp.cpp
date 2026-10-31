#include "GridToShp.h"
#include "Box.h"
#include "GridVertices.h"
#include "shpGridToLine.h"
#include <fstream>
#include <sstream>
using namespace std;

//corners, clockwise from the upper left, and center of a box
inline void boxCorners(const GridVertices & gv, Box * b, double * padfX, double * padfY,
                       double * ctrX, double * ctrY)
{
    gv.get_corners(b, padfX, padfY);
    const Point2d & ctr = gv.center(b);
    ctrX[0] = ctr[0];
    ctrY[0] = ctr[1];
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
        cerr << "! Error: only support exporting as polygon/point shape files!" << endl;
        return NULL;
    }
    return handler;
}

void writeDB(int entityNum, DBFHandle & dbf, vector<int> & fieldIds, vector<int> & dbIds, int num,
             int l, int r, int c, string & location, double top, double bot, double dx, double dy)
{
    //write to DB
    DBFWriteIntegerAttribute(dbf, entityNum, fieldIds[0], num);
    DBFWriteIntegerAttribute(dbf, entityNum, fieldIds[1], l);
    DBFWriteIntegerAttribute(dbf, entityNum, fieldIds[2], r);
    DBFWriteIntegerAttribute(dbf, entityNum, fieldIds[3], c);

    DBFWriteStringAttribute(dbf, entityNum, fieldIds[4], location.c_str());

    DBFWriteDoubleAttribute(dbf, entityNum, dbIds[0], top);
    DBFWriteDoubleAttribute(dbf, entityNum, dbIds[1], bot);
    DBFWriteDoubleAttribute(dbf, entityNum, dbIds[2], dx);
    DBFWriteDoubleAttribute(dbf, entityNum, dbIds[3], dy);
}

//write feature and info using index
void writeFeatEntity(Grid * grid, SHPHandle & handler, DBFHandle & dbf, vector<int> & fieldIds,
                     vector<int> & dbIds, vector<Box *> & b_vec, vector<int> & number_vec,
                     vector<int> & l_vec, vector<int> & r_vec, vector<int> & c_vec,
                     vector<string> & loc_vec, double cx, double cy, double rotation, int shptype,
                     bool one_base_number, vector<size_t> sort_idx, bool bConcide)
{
    QuadTree3D * qtree = dynamic_cast<QuadTree3D *>(grid);
    ModflowGrid * mfgrid = dynamic_cast<ModflowGrid *>(grid);
    bool isqtree = (qtree != NULL);

    double padfX[5] = {0}, padfY[5] = {0}, ctrX[1] = {0}, ctrY[1] = {0};
    GridVertices gv(grid);

    int n = sort_idx.size(), id;
    for (int i = 0; i < n; i++)
    {
        id = sort_idx[i];

        boxCorners(gv, b_vec[id], padfX, padfY, ctrX, ctrY);

        SHPObject * shpObj = NULL;
        int entityNum = -1;
        if (shptype == 1)
        {
            if (!bConcide)
                shpObj = SHPCreateSimpleObject(SHPT_POLYGON, 4, padfX, padfY, NULL);
            else
            {
                padfX[4] = padfX[0];
                padfY[4] = padfY[0];
                shpObj = SHPCreateSimpleObject(SHPT_POLYGON, 5, padfX, padfY, NULL);
            }

            entityNum = SHPWriteObject(handler, -1, shpObj);
        }
        else if (shptype == 0)
        {
            shpObj = SHPCreateSimpleObject(SHPT_POINT, 1, ctrX, ctrY, NULL);

            entityNum = SHPWriteObject(handler, -1, shpObj);
        }
        else
            assert(false);
        SHPDestroyObject(shpObj);

        int l = (one_base_number ? (l_vec[id] - 1) : l_vec[id]);
        int r = (one_base_number ? (r_vec[id] - 1) : r_vec[id]);
        int c = (one_base_number ? (c_vec[id] - 1) : c_vec[id]);

        double top, bot;
        if (isqtree)
        {
            top = qtree->top[b_vec[id]->number];
            bot = qtree->bot[b_vec[id]->number];
        }
        else
        {
            top = mfgrid->bot[l * mfgrid->nrow * mfgrid->ncol + r * mfgrid->ncol + c];
            bot = mfgrid->bot[(l + 1) * mfgrid->nrow * mfgrid->ncol + r * mfgrid->ncol + c];
        }

        writeDB(entityNum, dbf, fieldIds, dbIds, number_vec[id], l, r, c, loc_vec[id], top, bot,
                b_vec[id]->dx, b_vec[id]->dy);
    }
}

////write feature info

void DFSCollect(Grid * grid, Box * b, int l, int r, int c, string & location,
                bool one_based_numbering, vector<Box *> & b_vec, vector<int> & number_vec,
                vector<int> & l_vec, vector<int> & r_vec, vector<int> & c_vec,
                vector<string> & loc_vec, bool without_inactive)
{

    char tmp[8] = {'\0'};
    if (b->isLeaf)
    {
        if (without_inactive && (!b->active)) return;

        b_vec.push_back(b);
        number_vec.push_back(b->number == -1 ? b->number
                                             : (one_based_numbering ? b->number + 1 : b->number));
        l_vec.push_back(one_based_numbering ? (l + 1) : l);
        r_vec.push_back(one_based_numbering ? (r + 1) : r);
        c_vec.push_back(one_based_numbering ? (c + 1) : c);
        loc_vec.push_back(location);
    }
    else if (b->isLeaf == false) //b is not a leaf, call this function recursively
    {
        for (int i = 0; i < 4; i++)
        {
            int qid = i;
            if (i == 2) qid = 3;
            if (i == 3) qid = 2;

            sprintf(tmp, "%d", qid + 1);

            string mylocation = location;
            mylocation += tmp;

            DFSCollect(grid, b->pChildren[i], l, r, c, mylocation, one_based_numbering, b_vec,
                       number_vec, l_vec, r_vec, c_vec, loc_vec, without_inactive);
        }
    }
}

void writeQuadtree2Shp(QuadTree3D * qtree, string shpName, string feature_type, bool bConcide)
{
    if (feature_type == "line")
    {
        writeQuadtree2ShpLine(qtree, shpName);
        return;
    }

    //create a shape file
    int shptype = true; //polygon : true ; point : false
    SHPHandle handler = createShpHandle(shpName, feature_type, shptype);
    if (handler == NULL)
    {
        cerr << "! Error: failed to open shapefile: " << shpName << endl;
        return;
    }

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
    int baseLyrId = DBFAddField(dbfHandle, "layer", FTInteger, 20, 0);
    int baseRowId = DBFAddField(dbfHandle, "row", FTInteger, 20, 0);
    int baseColId = DBFAddField(dbfHandle, "col", FTInteger, 20, 0);
    int childLocId = DBFAddField(dbfHandle, "child_location", FTString, 50, 10);

    int tpId = DBFAddField(dbfHandle, "top", FTDouble, 20, 8);
    int btId = DBFAddField(dbfHandle, "bottom", FTDouble, 20, 8);
    int delrId = DBFAddField(dbfHandle, "delr", FTDouble, 20, 8);
    int delcId = DBFAddField(dbfHandle, "delc", FTDouble, 20, 8);

    vector<int> fieldIds;
    fieldIds.push_back(nodeFieldId);
    fieldIds.push_back(baseLyrId);
    fieldIds.push_back(baseRowId);
    fieldIds.push_back(baseColId);
    fieldIds.push_back(childLocId);

    vector<int> dbIds;
    dbIds.push_back(tpId);
    dbIds.push_back(btId);
    dbIds.push_back(delrId);
    dbIds.push_back(delcId);

    ////loop over all the leaf nodes

    ModflowGrid * mfgrid = qtree->getModflowGrid();

    double cx, cy, rotation;
    mfgrid->getRotatePara(cx, cy, rotation);

    vector<Box *> b_vec;
    vector<int> number_vec, l_vec, r_vec, c_vec;
    vector<string> loc_vec;
    //write shapefile data for each base grid cell
    for (int l = 0; l < mfgrid->nlay; l++)
    {
        for (int r = 0; r < mfgrid->nrow; r++)
        {
            for (int c = 0; c < mfgrid->ncol; c++)
            {
                int bid = mfgrid->get_nodeid(l, r, c);
                Box * box = qtree->nodegroup[bid];
                string location;

                DFSCollect(qtree, box, l, r, c, location, qtree->get_one_based_node_numbering(),
                           b_vec, number_vec, l_vec, r_vec, c_vec, loc_vec, true);
            }
        }
    }

    //sort
    std::vector<size_t> sort_idx;
    sort(number_vec, sort_idx);

    //write the sorted
    writeFeatEntity(qtree, handler, dbfHandle, fieldIds, dbIds, b_vec, number_vec, l_vec, r_vec,
                    c_vec, loc_vec, cx, cy, rotation, shptype,
                    qtree->get_one_based_node_numbering(), sort_idx, bConcide);

    //close the shape file
    DBFClose(dbfHandle);
    SHPClose(handler);
}

void writeModflowGrid2Shp(ModflowGrid * modflow, string gridShpName, string feature_type,
                          bool without_inactive, bool one_based_numbering, bool bConcide)
{
    if (feature_type == "line")
    {
        writeModflowGrid2ShpLine(modflow, gridShpName);
        return;
    }

    //create a shape file
    int shptype = 0; //polygon : true ; point : false
    SHPHandle handler = createShpHandle(gridShpName, feature_type, shptype);
    if (handler == NULL) return;

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
    int baseLyrId = DBFAddField(dbfHandle, "layer", FTInteger, 20, 0);
    int baseRowId = DBFAddField(dbfHandle, "row", FTInteger, 20, 0);
    int baseColId = DBFAddField(dbfHandle, "col", FTInteger, 20, 0);
    int childLocId = DBFAddField(dbfHandle, "child_location", FTString, 50, 10);

    int tpId = DBFAddField(dbfHandle, "top", FTDouble, 20, 8);
    int btId = DBFAddField(dbfHandle, "bottom", FTDouble, 20, 8);
    int delrId = DBFAddField(dbfHandle, "delr", FTDouble, 20, 8);
    int delcId = DBFAddField(dbfHandle, "delc", FTDouble, 20, 8);

    vector<int> fieldIds;
    fieldIds.push_back(nodeFieldId);
    fieldIds.push_back(baseLyrId);
    fieldIds.push_back(baseRowId);
    fieldIds.push_back(baseColId);
    fieldIds.push_back(childLocId);

    vector<int> dbIds;
    dbIds.push_back(tpId);
    dbIds.push_back(btId);
    dbIds.push_back(delrId);
    dbIds.push_back(delcId);

    ////loop over all the leaf nodes

    ModflowGrid * mfgrid = modflow;
    char extra[64];

    double cx, cy, rotation;
    modflow->getRotatePara(cx, cy, rotation);

    vector<Box *> b_vec;
    vector<int> number_vec, l_vec, r_vec, c_vec;
    vector<string> loc_vec;
    for (int l = 0; l < mfgrid->nlay; l++)
    {
        for (int r = 0; r < mfgrid->nrow; r++)
        {
            for (int c = 0; c < mfgrid->ncol; c++)
            {
                int bid = mfgrid->get_nodeid(l, r, c);
                Box * box = mfgrid->nodegroup[bid];
                string extra_str = extra;
                string location;

                DFSCollect(modflow, box, l, r, c, location, one_based_numbering, b_vec, number_vec,
                           l_vec, r_vec, c_vec, loc_vec, true);
            }
        }
    }

    //sort
    std::vector<size_t> sort_idx;
    sort(number_vec, sort_idx);

    //write the sorted
    writeFeatEntity(modflow, handler, dbfHandle, fieldIds, dbIds, b_vec, number_vec, l_vec, r_vec,
                    c_vec, loc_vec, cx, cy, rotation, shptype, one_based_numbering, sort_idx,
                    bConcide);

    //close the shape file
    DBFClose(dbfHandle);
    SHPClose(handler);
}
