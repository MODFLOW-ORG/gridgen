#include "shpGridToLine.h"
#include "GridVertices.h"
#include "GridgenError.h"
using namespace cusg;

//coordinates of a corner of a box, numbered clockwise from the upper left
void corner_point(const GridVertices & gv, const Box * b, int corner, double & x, double & y)
{
    const Point2d & p = gv.vertices()[gv.corners(b)[corner]];
    x = p[0];
    y = p[1];
}
void prepare_write_line(string grid_shp_name, SHPHandle * shp_handle, DBFHandle * dbf,
                        int * layer_field)
{
    *shp_handle = SHPCreate(grid_shp_name.c_str(), SHPT_ARC);
    if ((*shp_handle) == nullptr) throw GridgenError("cannot write shapefile " + grid_shp_name);

    //create a .dbf file
    string dbf_name = grid_shp_name; /*
	int lastL = grid_shp_name.find_last_of('/') == string::npos ? grid_shp_name.find;*/
    if (grid_shp_name.find('.') != string::npos)
    {
        int dot_index = grid_shp_name.find_last_of('.');
        dbf_name = grid_shp_name.substr(0, dot_index);
    }
    dbf_name += ".dbf";
    *dbf = DBFCreate(dbf_name.c_str());

    *layer_field = DBFAddField(*dbf, "layer", FTInteger, 20, 0);
}

void write_one_line(int nvertices, double * x, double * y, int layer_field, int layer,
                    SHPHandle handler, DBFHandle dbf_handle)
{
    SHPObject * shape = nullptr;
    //write to the shape file
    shape = SHPCreateSimpleObject(SHPT_ARC, nvertices, x, y, nullptr);
    int entity_num = SHPWriteObject(handler, -1, shape);
    //clear memory
    SHPDestroyObject(shape);

    //write the entity feature info
    DBFWriteIntegerAttribute(dbf_handle, entity_num, layer_field, layer);
}

//write the column and row edges of the base grid; grid holds the base cells
void write_one_mf_grid(cusg::Grid * grid, cusg::ModflowGrid * mfgrid, const GridVertices & gv,
                       SHPHandle & handler, DBFHandle & dbf_handle, int layer_field)
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
            Box * top = grid->node_group[mfgrid->get_box_index(k, 0, c)];
            Box * bot = grid->node_group[mfgrid->get_box_index(k, nrow - 1, c)];
            corner_point(gv, top, (i < ncol) ? 0 : 1, line_x[0], line_y[0]);
            corner_point(gv, bot, (i < ncol) ? 3 : 2, line_x[1], line_y[1]);

            write_one_line(2, line_x, line_y, layer_field, k, handler, dbf_handle);
        }

        //write rows, from the west edge to the east edge
        for (int i = 0; i <= nrow; i++)
        {
            int r = min(i, nrow - 1);
            Box * west = grid->node_group[mfgrid->get_box_index(k, r, 0)];
            Box * east = grid->node_group[mfgrid->get_box_index(k, r, ncol - 1)];
            corner_point(gv, west, (i < nrow) ? 0 : 3, line_x[0], line_y[0]);
            corner_point(gv, east, (i < nrow) ? 1 : 2, line_x[1], line_y[1]);

            write_one_line(2, line_x, line_y, layer_field, k, handler, dbf_handle);
        }
    }
}
void write_line_shapefile(cusg::QuadTree3D * qtree, string grid_shp_name)
{
    SHPHandle handler;
    DBFHandle dbf_handle;
    int layer_field;
    prepare_write_line(grid_shp_name, &handler, &dbf_handle, &layer_field);

    ModflowGrid * mfgrid = qtree->get_modflow_grid();

    GridVertices gv(qtree);
    write_one_mf_grid(qtree, mfgrid, gv, handler, dbf_handle, layer_field);

    double line_x[2], line_y[2];
    //get through each cell
    for (auto b : qtree->node_group)
    {
        if (b->is_leaf == false)
        {
            //from the south edge to the north edge through the center
            corner_point(gv, b->p_children[3], 2, line_x[0], line_y[0]);
            corner_point(gv, b->p_children[0], 1, line_x[1], line_y[1]);

            write_one_line(2, line_x, line_y, layer_field, b->layer, handler, dbf_handle);

            ////////////////////////////////////////
            //from the west edge to the east edge through the center
            corner_point(gv, b->p_children[3], 0, line_x[0], line_y[0]);
            corner_point(gv, b->p_children[2], 1, line_x[1], line_y[1]);

            write_one_line(2, line_x, line_y, layer_field, b->layer, handler, dbf_handle);
        }
    }

    DBFClose(dbf_handle);
    SHPClose(handler);
}

void write_line_shapefile(cusg::ModflowGrid * mfgrid, string grid_shp_name)
{
    SHPHandle handler;
    DBFHandle dbf_handle;
    int layer_field;
    prepare_write_line(grid_shp_name, &handler, &dbf_handle, &layer_field);

    GridVertices gv(mfgrid);
    write_one_mf_grid(mfgrid, mfgrid, gv, handler, dbf_handle, layer_field);

    DBFClose(dbf_handle);
    SHPClose(handler);
}