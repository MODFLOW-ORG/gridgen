#include "shpGridToLine.h"
#include "GridVertices.h"
using namespace cusg;

//coordinates of a corner of a box, numbered clockwise from the upper left
void cornerPoint(const GridVertices & gv, const Box * b, int corner, double & x, double & y)
{
    const Point2d & p = gv.vertices()[gv.corners(b)[corner]];
    x = p[0];
    y = p[1];
}
void prepareWriteLine(string gridShpName, SHPHandle * shp_handle, DBFHandle * dbf,
                      int * layer_field)
{
    *shp_handle = SHPCreate(gridShpName.c_str(), SHPT_ARC);
    if ((*shp_handle) == nullptr)
    {
        cerr << "can't create shape file: " << gridShpName.c_str() << endl;
        return;
    }

    //create a .dbf file
    string dbfName = gridShpName; /*
	int lastL = gridShpName.find_last_of('/') == string::npos ? gridShpName.find;*/
    if (gridShpName.find('.') != string::npos)
    {
        int dotIndex = gridShpName.find_last_of('.');
        dbfName = gridShpName.substr(0, dotIndex);
    }
    dbfName += ".dbf";
    *dbf = DBFCreate(dbfName.c_str());

    *layer_field = DBFAddField(*dbf, "layer", FTInteger, 20, 0);
}

void writeOneLine(int nvertices, double * x, double * y, int layer_field, int layer,
                  SHPHandle handler, DBFHandle dbfHandle)
{
    SHPObject * shape = nullptr;
    //write to the shape file
    shape = SHPCreateSimpleObject(SHPT_ARC, nvertices, x, y, nullptr);
    int entityNum = SHPWriteObject(handler, -1, shape);
    //clear memory
    SHPDestroyObject(shape);

    //write the entity feature info
    DBFWriteIntegerAttribute(dbfHandle, entityNum, layer_field, layer);
}

//write the column and row edges of the base grid; grid holds the base cells
void writeOneMFGrid(cusg::Grid * grid, cusg::ModflowGrid * mfgrid, const GridVertices & gv,
                    SHPHandle & handler, DBFHandle & dbfHandle, int layer_field)
{
    double line_x[2], line_y[2];
    int nrow = mfgrid->nrow;
    int ncol = mfgrid->ncol;

    for (int k = 0; k < mfgrid->nlay; k++)
    {
        //write columns, from the top of the first row to the bottom of the last
        for (int i = 0; i <= ncol; i++)
        {
            int c = min(i, ncol - 1);
            Box * top = grid->nodegroup[mfgrid->get_nodeid(k, 0, c)];
            Box * bot = grid->nodegroup[mfgrid->get_nodeid(k, nrow - 1, c)];
            cornerPoint(gv, top, (i < ncol) ? 0 : 1, line_x[0], line_y[0]);
            cornerPoint(gv, bot, (i < ncol) ? 3 : 2, line_x[1], line_y[1]);

            writeOneLine(2, line_x, line_y, layer_field, k, handler, dbfHandle);
        }

        //write rows, from the west edge to the east edge
        for (int i = 0; i <= nrow; i++)
        {
            int r = min(i, nrow - 1);
            Box * west = grid->nodegroup[mfgrid->get_nodeid(k, r, 0)];
            Box * east = grid->nodegroup[mfgrid->get_nodeid(k, r, ncol - 1)];
            cornerPoint(gv, west, (i < nrow) ? 0 : 3, line_x[0], line_y[0]);
            cornerPoint(gv, east, (i < nrow) ? 1 : 2, line_x[1], line_y[1]);

            writeOneLine(2, line_x, line_y, layer_field, k, handler, dbfHandle);
        }
    }
}
void writeQuadtree2ShpLine(cusg::QuadTree3D * qtree, string gridShpName)
{
    SHPHandle handler;
    DBFHandle dbfHandle;
    int layer_field;
    prepareWriteLine(gridShpName, &handler, &dbfHandle, &layer_field);

    ModflowGrid * mfgrid = qtree->getModflowGrid();

    GridVertices gv(qtree);
    writeOneMFGrid(qtree, mfgrid, gv, handler, dbfHandle, layer_field);

    double line_x[2], line_y[2];
    //get through each cell
    for (auto b : qtree->nodegroup)
    {
        if (b->isLeaf == false)
        {
            //from the south edge to the north edge through the center
            cornerPoint(gv, b->pChildren[3], 2, line_x[0], line_y[0]);
            cornerPoint(gv, b->pChildren[0], 1, line_x[1], line_y[1]);

            writeOneLine(2, line_x, line_y, layer_field, b->layer, handler, dbfHandle);

            ////////////////////////////////////////
            //from the west edge to the east edge through the center
            cornerPoint(gv, b->pChildren[3], 0, line_x[0], line_y[0]);
            cornerPoint(gv, b->pChildren[2], 1, line_x[1], line_y[1]);

            writeOneLine(2, line_x, line_y, layer_field, b->layer, handler, dbfHandle);
        }
    }

    DBFClose(dbfHandle);
    SHPClose(handler);
}

void writeModflowGrid2ShpLine(cusg::ModflowGrid * mfgrid, string gridShpName)
{
    SHPHandle handler;
    DBFHandle dbfHandle;
    int layer_field;
    prepareWriteLine(gridShpName, &handler, &dbfHandle, &layer_field);

    GridVertices gv(mfgrid);
    writeOneMFGrid(mfgrid, mfgrid, gv, handler, dbfHandle, layer_field);

    DBFClose(dbfHandle);
    SHPClose(handler);
}