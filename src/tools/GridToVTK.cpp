
#include "GridToVTK.h"
#include "GridVertices.h"
#include "RoundTrip.h"
#include <cfloat>
#include <fstream>
#include <sstream>
using namespace std;

static void write_vtk_separate_vertices(ModflowGrid * modflow, string vtk_name,
                                        bool without_inactive, bool one_based_numbering);
static void write_vtk_separate_vertices(QuadTree3D * qtree, string vtk_name);
static void write_vtk_shared_vertices(ModflowGrid * mfgrid, string vtk_name, bool without_inactive,
                                      bool one_based_numbering);
static void write_vtk_shared_vertices(QuadTree3D * qtree, string vtk_name);

void write_vtk(ModflowGrid * modflow, string vtk_name, bool share_vertices, bool without_inactive,
               bool one_based_numbering)
{
    if (share_vertices)
        write_vtk_shared_vertices(modflow, vtk_name, without_inactive, one_based_numbering);
    else
        write_vtk_separate_vertices(modflow, vtk_name, without_inactive, one_based_numbering);
}

void write_vtk(QuadTree3D * qtree, string vtk_name, bool share_vertices)
{
    if (share_vertices)
        write_vtk_shared_vertices(qtree, vtk_name);
    else
        write_vtk_separate_vertices(qtree, vtk_name);
}

//reverse interpolate the z coordinates for all boxes
//void reverseInterpolate(QuadTree3D* qtree);

/********************************************************/
//write the vtk file
void write_vtk_header(ofstream & ofile)
{
    ofile << "<?xml version=\"1.0\"?>\n";
    ofile << "<VTKFile type=\"UnstructuredGrid\" version=\"0.1\" byte_order=\"BigEndian\">\n";
    ofile << "<UnstructuredGrid GhostLevel=\"0\">\n";
}
void write_vtk_tail(ofstream & ofile)
{
    ofile << "</UnstructuredGrid>\n";
    ofile << "</VTKFile>\n";
}

void write_vtk_cell(ofstream & ofile, int nboxes)
{
    /***************************************************************************************/
    //write the cells
    ofile << "<Cells>\n";
    //write the connectivity
    ofile << "<DataArray type=\"Int32\" Name=\"connectivity\" format=\"ascii\">\n";
    int sid = 0, id = 0;
    for (int i = 0; i < nboxes; i++)
    {
        ofile << (sid) << "\t" << (sid + 6) << "\t" << (sid + 2) << "\t" << (sid + 4) << "\t"
              << (sid + 1) << "\t" << (sid + 7) << "\t" << (sid + 3) << "\t" << (sid + 5) << "\n";

        sid += 8;

        id++;
        if (id % 500 == 0) ofile.flush();
    }
    ofile << "</DataArray>\n";
    //write the offset
    sid = id = 0;
    sid = 8;
    ofile << "<DataArray type=\"Int32\" Name=\"offsets\" format=\"ascii\">\n";
    for (int i = 0; i < nboxes; i++)
    {
        ofile << sid << "\t";
        sid += 8;

        id++;
        if (id % 1000 == 0) ofile.flush();
    }
    ofile << "\n";
    ofile << "</DataArray>\n";
    //write the types
    id = 0;
    ofile << "<DataArray type=\"UInt8\" Name=\"types\" format=\"ascii\">\n";
    for (int i = 0; i < nboxes; i++)
    {
        ofile << "11\t";
        id++;
        if (id % 1000 == 0) ofile.flush();
    }
    ofile << "</DataArray>\n";
    ofile << "</Cells>\n";
}

template <typename Dtype>
void write_data_array(ofstream & ofile, string type, string name, const vector<Dtype> & data)
{
    ofile << "<DataArray type=\"" << type << "\" Name=\"" << name << "\" format=\"ascii\">";
    int id = 0;
    for (typename vector<Dtype>::const_iterator dit = data.begin(); dit != data.end(); dit++)
    {
        ofile << round_trip_value(*dit) << "\t";

        if (id % 1000 == 0) ofile.flush();
        id++;
    }
    ofile << "</DataArray>\n";
}

//export modflow grid to a paraview vtk file
static void write_vtk_separate_vertices(ModflowGrid * modflow, string vtk_name,
                                        bool without_inactive, bool one_based_numbering)
{
    GridVertices gv(modflow);
    //write the points, loop over each layer
    int nlay = modflow->nlay;
    int nrow = modflow->nrow;
    int ncol = modflow->ncol;

    double corner_x[4] = {0};
    double corner_y[4] = {0};

    stringstream ss0;
    ss0 << vtk_name.c_str();
    ss0 << ".vtu";
    ofstream ofile(ss0.str().c_str());
    write_vtk_header(ofile);

    /***************************************************************************************/
    int total_cell = nlay * nrow * ncol;
    ofile << "<Piece NumberOfPoints=\"" << total_cell * 8 << "\" NumberOfCells=\"" << total_cell
          << "\" >";

    vector<int> box_indices, cells, layers, rows, cols, vtkids;
    vector<double> tops, bottoms, delrs, delcs;

    //write the points
    ofile << "<Points>\n";
    ofile << "<DataArray type=\"Float64\" NumberOfComponents=\"3\" format=\"ascii\">\n";
    int id = 0, cell_num = 0;
    for (int l = 0; l < nlay; l++)
    {
        for (int r = 0; r < nrow; r++)
        {
            for (int c = 0; c < ncol; c++)
            {
                int bid = modflow->get_box_index(l, r, c);
                Box * box = modflow->node_group[bid];

                //if(box->is_leaf && box->active)
                if (box->is_leaf)
                {
                    if (without_inactive && (!box->active)) continue;

                    double top = modflow->bot[l * nrow * ncol + r * ncol + c];
                    double bot = modflow->bot[(l + 1) * nrow * ncol + r * ncol + c];

                    gv.get_corners(box, corner_x, corner_y);
                    for (int i = 0; i < 4; i++)
                    {
                        ofile << RoundTrip(corner_x[i]) << "\t" << RoundTrip(corner_y[i]) << "\t"
                              << RoundTrip(top * ZSCALE) << "\t"; //(box->z - box->dz / 2)<<"\t";
                        ofile << RoundTrip(corner_x[i]) << "\t" << RoundTrip(corner_y[i]) << "\t"
                              << RoundTrip(bot * ZSCALE) << "\t"; //(box->z + box->dz / 2)<<"\t";

                        id += 2;
                        if (id % 500 == 0) ofile.flush();
                    }

                    //store info
                    if (!one_based_numbering)
                    {
                        box_indices.push_back(box->id);
                        cells.push_back(box->number);
                        layers.push_back(l);
                        rows.push_back(r);
                        cols.push_back(c);
                    }
                    else
                    {
                        box_indices.push_back(box->id + 1);
                        cells.push_back(box->number + 1);
                        layers.push_back(l + 1);
                        rows.push_back(r + 1);
                        cols.push_back(c + 1);
                    }

                    tops.push_back(top);
                    bottoms.push_back(bot);
                    delrs.push_back(box->dx);
                    delcs.push_back(box->dy);

                    vtkids.push_back(cell_num);
                    cell_num++;
                }
            }
        }
    }
    ofile << "\n";
    ofile << "</DataArray>\n";
    ofile << "</Points>\n";

    //write cells
    write_vtk_cell(ofile, id / 8);

    //write the attributes
    ofile << "<CellData Scalars=\"scalars\">\n";
    write_data_array(ofile, "Int32", "id", box_indices);
    write_data_array(ofile, "Int32", "nodenumber", cells);
    write_data_array(ofile, "Int32", "layer", layers);
    write_data_array(ofile, "Int32", "row", rows);
    write_data_array(ofile, "Int32", "col", cols);
    write_data_array(ofile, "Int32", "vtkid", vtkids);

    write_data_array(ofile, "Float64", "top", tops);
    write_data_array(ofile, "Float64", "bottom", bottoms);
    write_data_array(ofile, "Float64", "delr", delrs);
    write_data_array(ofile, "Float64", "delc", delcs);
    ofile << "</CellData>\n";

    ofile << "</Piece>\n";
    write_vtk_tail(ofile);

    ofile.close();
}

void dfs_write_quadtree_vtk(QuadTree3D * qtree, ofstream & ofile, int & point_count, Box * b,
                            int layer, int row, int col, bool one_based_numbering,
                            const GridVertices & gv, vector<int> * int_values,
                            vector<double> * double_values)
{
    //if(b->is_leaf)
    if (b->is_leaf && b->active)
    {
        double corner_x[4] = {0};
        double corner_y[4] = {0};

        gv.get_corners(b, corner_x, corner_y);

        double top = qtree->top[b->number];
        double bot = qtree->bot[b->number];

        for (int i = 0; i < 4; i++)
        {
            ofile << RoundTrip(corner_x[i]) << "\t" << RoundTrip(corner_y[i]) << "\t"
                  << RoundTrip(top * ZSCALE) << "\t";
            ofile << RoundTrip(corner_x[i]) << "\t" << RoundTrip(corner_y[i]) << "\t"
                  << RoundTrip(bot * ZSCALE) << "\t";

            point_count += 2;
            if (point_count % 500 == 0) ofile.flush();
        }
        //store info

        if (!one_based_numbering)
        {
            int_values[0].push_back(b->id);
            int_values[1].push_back(b->number);
            int_values[2].push_back(layer);
            int_values[3].push_back(row);
            int_values[4].push_back(col);
        }
        else
        {
            int_values[0].push_back(b->id + 1);
            int_values[1].push_back(b->number + 1);
            int_values[2].push_back(layer + 1);
            int_values[3].push_back(row + 1);
            int_values[4].push_back(col + 1);
        }

        double_values[0].push_back(top);
        double_values[1].push_back(bot);
        double_values[2].push_back(b->dx);
        double_values[3].push_back(b->dy);

        int_values[5].push_back(point_count); //the vtk_id
    }
    else if (b->is_leaf == false) //b is not a leaf, call this function recursively
    {
        for (auto & i : b->p_children)
        {
            dfs_write_quadtree_vtk(qtree, ofile, point_count, i, layer, row, col,
                                   one_based_numbering, gv, int_values, double_values);
        }
    } //end b->is_leaf
}

//write the Quadtree grid
static void write_vtk_separate_vertices(QuadTree3D * qtree, string vtk_name)
{
    //make sure the number of the nodes are correct
    qtree->number_nodes();

    ModflowGrid * mfgrid = qtree->get_modflow_grid();
    GridVertices gv(qtree);

    //write the points, loop over each layer
    stringstream ss0;
    ss0 << vtk_name.c_str();
    ss0 << ".vtu";
    ofstream ofile(ss0.str().c_str());
    write_vtk_header(ofile);

    /***************************************************************************************/
    int total_cell = qtree->nodes();
    ofile << "<Piece NumberOfPoints=\"" << total_cell * 8 << "\" NumberOfCells=\"" << total_cell
          << "\" >";

    vector<string> int_names, double_names;
    vector<int> int_values[6]; // box_index, cells, layers, rows, cols, vtkid;
    int_names.push_back("nodeid");
    int_names.push_back("nodenumber");
    int_names.push_back("layer");
    int_names.push_back("row");
    int_names.push_back("col");
    int_names.push_back("vtkid");

    vector<double> double_values[4]; //tops, bottoms, delrs, delcs;
    double_names.push_back("top");
    double_names.push_back("bottom");
    double_names.push_back("delr");
    double_names.push_back("delc");

    //write the points
    ofile << "<Points>\n";
    ofile << "<DataArray type=\"Float64\" NumberOfComponents=\"3\" format=\"ascii\">\n";
    int point_count = 0;
    for (int l = 0; l < mfgrid->nlay; l++)
    {
        for (int r = 0; r < mfgrid->nrow; r++)
        {
            for (int c = 0; c < mfgrid->ncol; c++)
            {
                int bid = mfgrid->get_box_index(l, r, c);
                Box * box = qtree->node_group[bid];

                if (box->active)
                    dfs_write_quadtree_vtk(qtree, ofile, point_count, box, l, r, c,
                                           qtree->get_one_based_node_numbering(), gv, int_values,
                                           double_values);
            }
        }
    }

    ofile << "\n";
    ofile << "</DataArray>\n";
    ofile << "</Points>\n";

    //write cells
    write_vtk_cell(ofile, point_count / 8);

    //write the attributes
    ofile << "<CellData Scalars=\"scalars\">\n";
    int int_size = int_names.size();
    for (int i = 0; i < int_size; i++)
    {
        write_data_array(ofile, "Int32", int_names[i], int_values[i]);
    }

    for (int i = 0; i < 4; i++)
    {
        write_data_array(ofile, "Float64", double_names[i], double_values[i]);
    }
    ofile << "</CellData>\n";

    ofile << "</Piece>\n";
    write_vtk_tail(ofile);

    ofile.close();
}

void write_cells(ofstream & ofile, vector<Box *> & all_boxes, vector<int> * int_values,
                 vector<double> * double_values, BoxVtx * box_vertices)
{
    /***************************************************************************************/
    //write the cells
    ofile << "<Cells>\n";
    //write the connectivity
    ofile << "<DataArray type=\"Int32\" Name=\"connectivity\" format=\"ascii\">\n";

    int rid = 0;
    int cell_num = all_boxes.size();
    for (auto b : all_boxes)
    {
        ofile << box_vertices[b->id].vertices[0]->id << "\t" << box_vertices[b->id].vertices[1]->id
              << "\t" << box_vertices[b->id].vertices[3]->id << "\t"
              << box_vertices[b->id].vertices[2]->id << "\t" << box_vertices[b->id].vertices[4]->id
              << "\t" << box_vertices[b->id].vertices[5]->id << "\t"
              << box_vertices[b->id].vertices[7]->id << "\t" << box_vertices[b->id].vertices[6]->id
              << "\n";

        rid++;
        if (rid % 5000 == 0) ofile.flush();
    }
    ofile << "</DataArray>\n";

    //write the offset
    int sid, id = 0;
    sid = 8;
    ofile << "<DataArray type=\"Int32\" Name=\"offsets\" format=\"ascii\">\n";
    for (int i = 0; i < cell_num; i++)
    {
        ofile << sid << "\t";
        sid += 8;

        id++;
        if (id % 5000 == 0) ofile.flush();
    }
    ofile << "\n";
    ofile << "</DataArray>\n";
    //write the types
    id = 0;
    ofile << "<DataArray type=\"UInt8\" Name=\"types\" format=\"ascii\">\n";
    for (int i = 0; i < cell_num; i++)
    {
        ofile << "11\t";
        id++;
        if (id % 5000 == 0) ofile.flush();
    }
    ofile << "</DataArray>\n";
    ofile << "</Cells>\n";
}
/******************************************************************************************************************
* write VTK file with sharing vertices for modflow grid
*******************************************************************************************************************/
void build_sharing_vertices_mfgrid(ModflowGrid * mfgrid, vector<Vertex *> & vdata,
                                   BoxVtx * box_vertices)
{
    double z;
    GridVertices gv(mfgrid);

    int mfgrid_csize = (mfgrid->nrow) * (mfgrid->ncol);

    //get all the Xe and Ye
    vector<double> xe = mfgrid->get_local_xe_array();
    vector<double> ye = mfgrid->get_local_ye_array();

    //use the edges to locate the box corners
    int vid = 0;
    for (int k = 0; k <= mfgrid->nlay; k++)
    {
        double * zgrid = &(mfgrid->bot[mfgrid_csize * k]);
        ModflowAsciiGrid ag(mfgrid, zgrid, mfgrid->ncol, mfgrid->nrow);

        for (int r = 0; r <= mfgrid->nrow; r++)
        {
            for (int c = 0; c <= mfgrid->ncol; c++)
            {
                Point2d pos(xe[c], ye[r]);
                z = QuadTree3D::bilinear_interpolation(pos, &ag);

                //the corner of the cell to the southeast, or of the last row or column
                Box * box = mfgrid->node_group[mfgrid->get_box_index(0, min(r, mfgrid->nrow - 1),
                                                                     min(c, mfgrid->ncol - 1))];
                int corner = (r < mfgrid->nrow) ? ((c < mfgrid->ncol) ? 0 : 1)
                                                : ((c < mfgrid->ncol) ? 3 : 2);
                const Point2d & xy = gv.vertices()[gv.corners(box)[corner]];

                Vertex * nv = new Vertex(xy[0], xy[1], z, vid++);
                vdata.push_back(nv);
            }
        }
    }

    int mfgrid_esize = (mfgrid->nrow + 1) * (mfgrid->ncol + 1);
    int perrow_esize = mfgrid->ncol + 1;
    //new assign the vertex pointer to boxes
    for (int k = 0; k < mfgrid->nlay; k++)
    {
        for (int r = 0; r < mfgrid->nrow; r++)
        {
            for (int c = 0; c < mfgrid->ncol; c++)
            {
                int bid = mfgrid->get_box_index(k, r, c);
                Box * box = mfgrid->node_group[bid];
                box_vertices[box->id].vertices[0] = vdata[k * mfgrid_esize + r * perrow_esize + c];
                box_vertices[box->id].vertices[1] =
                    vdata[k * mfgrid_esize + (r + 1) * perrow_esize + c];
                box_vertices[box->id].vertices[2] =
                    vdata[k * mfgrid_esize + (r + 1) * perrow_esize + c + 1];
                box_vertices[box->id].vertices[3] =
                    vdata[k * mfgrid_esize + r * perrow_esize + c + 1];

                box_vertices[box->id].vertices[4] =
                    vdata[(k + 1) * mfgrid_esize + r * perrow_esize + c];
                box_vertices[box->id].vertices[5] =
                    vdata[(k + 1) * mfgrid_esize + (r + 1) * perrow_esize + c];
                box_vertices[box->id].vertices[6] =
                    vdata[(k + 1) * mfgrid_esize + (r + 1) * perrow_esize + c + 1];
                box_vertices[box->id].vertices[7] =
                    vdata[(k + 1) * mfgrid_esize + r * perrow_esize + c + 1];
            }
        }
    }
}
int collect_cell_share(ModflowGrid * mfgrid, vector<Box *> & all_boxes, vector<int> * int_values,
                       vector<double> * double_values, bool without_inactive,
                       bool one_based_numbering)
{
    int cell_num = 0;
    for (int k = 0; k < mfgrid->nlay; k++)
    {
        for (int r = 0; r < mfgrid->nrow; r++)
        {
            for (int c = 0; c < mfgrid->ncol; c++)
            {
                int bid = mfgrid->get_box_index(k, r, c);
                Box * b = mfgrid->node_group[bid];
                if (b->is_leaf)
                {
                    if (without_inactive && (!b->active)) continue;

                    all_boxes.push_back(b);

                    //store info
                    if (one_based_numbering)
                    {
                        int_values[0].push_back(b->id + 1);
                        int_values[1].push_back(b->number + 1);
                        int_values[2].push_back(k + 1);
                        int_values[3].push_back(r + 1);
                        int_values[4].push_back(c + 1);
                    }
                    else
                    {
                        int_values[0].push_back(b->id);
                        int_values[1].push_back(b->number);
                        int_values[2].push_back(k);
                        int_values[3].push_back(r);
                        int_values[4].push_back(c);
                    }

                    double top =
                        mfgrid->bot[k * mfgrid->nrow * mfgrid->ncol + r * mfgrid->ncol + c];
                    double bot =
                        mfgrid->bot[(k + 1) * mfgrid->nrow * mfgrid->ncol + r * mfgrid->ncol + c];

                    double_values[0].push_back(top);
                    double_values[1].push_back(bot);
                    double_values[2].push_back(b->dx);
                    double_values[3].push_back(b->dy);

                    int_values[5].push_back(cell_num++); //the vtk_id
                }
            }
        }
    }
    return cell_num;
}

static void write_vtk_shared_vertices(ModflowGrid * mfgrid, string vtk_name, bool without_inactive,
                                      bool one_based_numbering)
{
    int nboxes = mfgrid->node_group.size();
    vector<BoxVtx> box_vertex_data(nboxes);
    BoxVtx * box_vertices = box_vertex_data.data();
    vector<Vertex *> vdata;
    build_sharing_vertices_mfgrid(mfgrid, vdata, box_vertices);
    //begin to write

    /*****************************************************************************************/
    //collect all the boxes and necessary info
    vector<Box *> all_boxes;

    vector<string> int_names, double_names;
    vector<int> int_values[6]; // box_index, cells, layers, rows, cols, vtkid;
    int_names.push_back("nodeid");
    int_names.push_back("nodenumber");
    int_names.push_back("layer");
    int_names.push_back("row");
    int_names.push_back("col");
    int_names.push_back("vtkid");

    vector<double> double_values[4]; //tops, bottoms, delrs, delcs;
    double_names.push_back("top");
    double_names.push_back("bottom");
    double_names.push_back("delr");
    double_names.push_back("delc");

    collect_cell_share(mfgrid, all_boxes, int_values, double_values, without_inactive,
                       one_based_numbering);

    //write the points, loop over each layer
    stringstream ss0;
    ss0 << vtk_name.c_str();
    ss0 << ".vtu";
    ofstream ofile(ss0.str().c_str());
    write_vtk_header(ofile);

    /***************************************************************************************/
    ofile << "<Piece NumberOfPoints=\"" << vdata.size() << "\" NumberOfCells=\"" << all_boxes.size()
          << "\" >";

    //write the points
    ofile << "<Points>\n";
    ofile << "<DataArray type=\"Float64\" NumberOfComponents=\"3\" format=\"ascii\">\n";
    int flash = 0;
    for (auto nv : vdata)
    {
        ofile << RoundTrip(nv->pos[0]) << "\t" << RoundTrip(nv->pos[1]) << "\t"
              << RoundTrip(nv->pos[2] * ZSCALE) << "\t";

        if (flash % 5000 == 0)
        {
            ofile.flush();
        }
        flash++;
    }
    ofile << "\n";
    ofile << "</DataArray>\n";
    ofile << "</Points>\n";

    //write the boxes
    write_cells(ofile, all_boxes, int_values, double_values, box_vertices);

    //write the attributes
    ofile << "<CellData Scalars=\"scalars\">\n";
    int int_size = int_names.size();
    for (int i = 0; i < int_size; i++)
    {
        write_data_array(ofile, "Int32", int_names[i], int_values[i]);
    }

    for (int i = 0; i < 4; i++)
    {
        write_data_array(ofile, "Float64", double_names[i], (double_values[i]));
    }
    ofile << "</CellData>\n";

    ofile << "</Piece>\n";
    write_vtk_tail(ofile);

    ofile.close();

    //clear memory
    for (auto & vit : vdata)
    {
        delete vit;
    }
}

/******************************************************************************************************************
* write VTK file with sharing vertices for quadtree
******************************************************************************************************************/
//void DFSCollectInfo(Box* b, int l, int r, int c, int& nboxes,
void dfs_collect_cell_share(Box * b, int layer, int row, int col, QuadTree3D * qtree,
                            vector<Box *> & all_boxes, int & nboxes, vector<int> * int_values,
                            vector<double> * double_values, bool one_based_numbering)
{
    if (b->is_leaf && b->active)
    {
        all_boxes.push_back(b);

        //store info
        if (one_based_numbering)
        {
            int_values[0].push_back(b->id + 1);
            int_values[1].push_back(b->number + 1);
            int_values[2].push_back(layer + 1);
            int_values[3].push_back(row + 1);
            int_values[4].push_back(col + 1);
        }
        else
        {
            int_values[0].push_back(b->id);
            int_values[1].push_back(b->number);
            int_values[2].push_back(layer);
            int_values[3].push_back(row);
            int_values[4].push_back(col);
        }

        double top = qtree->top[b->number];
        double bot = qtree->bot[b->number];

        double_values[0].push_back(top);
        double_values[1].push_back(bot);
        double_values[2].push_back(b->dx);
        double_values[3].push_back(b->dy);

        int_values[5].push_back(nboxes++); //the vtk_id
    }
    else if (b->is_leaf == false)
    {
        for (auto & i : b->p_children)
        {
            dfs_collect_cell_share(i, layer, row, col, qtree, all_boxes, nboxes, int_values,
                                   double_values, one_based_numbering);
        }
    }
}
int collect_cell_share(QuadTree3D * qtree, vector<Box *> & all_boxes, vector<int> * int_values,
                       vector<double> * double_values, bool one_based_numbering)
{
    ModflowGrid * mfgrid = qtree->get_modflow_grid();
    int cell_num = 0;
    for (int k = 0; k < mfgrid->nlay; k++)
    {
        for (int r = 0; r < mfgrid->nrow; r++)
        {
            for (int c = 0; c < mfgrid->ncol; c++)
            {
                int bid = mfgrid->get_box_index(k, r, c);
                Box * box = qtree->node_group[bid];
                if (box->active)
                    dfs_collect_cell_share(box, k, r, c, qtree, all_boxes, cell_num, int_values,
                                           double_values, one_based_numbering);
            }
        }
    }
    return cell_num;
}

void rotate_one_pnt(double x, double y, double cx, double cy, double rotation, double & rotated_x,
                    double & rotated_y)
{
    rotated_x = cx + (x - cx) * cos(rotation) - (y - cy) * sin(rotation);
    rotated_y = cy + (y - cy) * cos(rotation) + (x - cx) * sin(rotation);
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////
//fit a plane using 3 points
double fit_z(double ** points, double x, double y)
{
    double a = (points[1][1] - points[0][1]) * (points[2][2] - points[0][2]) -
               (points[2][1] - points[0][1]) * (points[1][2] - points[0][2]);
    double b = (points[1][2] - points[0][2]) * (points[2][0] - points[0][0]) -
               (points[2][2] - points[0][2]) * (points[1][0] - points[0][0]);
    double c = (points[1][0] - points[0][0]) * (points[2][1] - points[0][1]) -
               (points[2][0] - points[0][0]) * (points[1][1] - points[0][1]);
    double d = -(a * points[0][0] + b * points[0][1] + c * points[0][2]);

    assert(fabs(c) >= 1e-10);
    double z = -(a * x + b * y + d) / c;

    return z;
}

//-------------------------------------------------------------------------------------------------

//
// for QuadTree3D sharing Vertices
//

struct VtkVertex
{
    VtkVertex()
    {
        number = -1;
        for (auto & i : box) i = nullptr;
    }

    Point3d pos;
    Box * box[8];
    int number;
};

/*		3
  _____________
		| e31
		|
0		*		2
		|
		| e13
  -------------
		1 
map: map<3,1>--> e31
*/
struct VtkEdge
{
    VtkEdge()
    {
        child_v = nullptr;
        v[0] = nullptr;
        v[1] = nullptr;
        child_e[0] = nullptr;
        child_e[1] = nullptr;
    }
    VtkEdge(VtkVertex * v1, VtkVertex * v2)
    {
        child_v = nullptr;
        v[0] = v1;
        v[1] = v2;
        child_e[0] = nullptr;
        child_e[1] = nullptr;
    }

    VtkVertex * child_v;
    VtkEdge * child_e[2];
    VtkVertex * v[2];

    //
    map<VtkEdge *, VtkEdge *> sub_e;
};

struct VtkBox
{
    VtkBox()
    {
        for (auto & i : e) i = nullptr;
    }
    VtkEdge * e[8];
};

//
//	build
//
void build_ve(QuadTree3D * qtree, vector<VtkVertex *> & vdata, vector<VtkEdge *> & edata,
              VtkBox * boxdata)
{
    ModflowGrid * mfgrid = qtree->get_modflow_grid();
    vector<double> xe = mfgrid->get_local_xe_array();
    vector<double> ye = mfgrid->get_local_ye_array();

    //////////////////////////////////////////////////////////////////////////////
    //initialize sharing v
    for (int k = 0; k <= mfgrid->nlay; k++)
    {
        for (int r = 0; r <= mfgrid->nrow; r++)
        {
            for (int c = 0; c <= mfgrid->ncol; c++)
            {
                VtkVertex * nv = new VtkVertex();
                nv->pos[0] = xe[c];
                nv->pos[1] = ye[r];
                vdata.push_back(nv);
            }
        }
    }
    //////////////////////////////////////////////////////////////////////////////
    //initializing the sharing edge
    int mfgrid_vsize = (mfgrid->nrow + 1) * (mfgrid->ncol + 1);
    int perrow_vsize = mfgrid->ncol + 1;
    for (int k = 0; k <= mfgrid->nlay; k++)
    {
        //initializing the row edges
        for (int r = 0; r <= mfgrid->nrow; r++)
        {
            for (int c = 0; c < mfgrid->ncol; c++)
            {
                VtkVertex * v1 = vdata[k * mfgrid_vsize + r * perrow_vsize + c];
                VtkVertex * v2 = vdata[k * mfgrid_vsize + r * perrow_vsize + c + 1];
                VtkEdge * e = new VtkEdge(v1, v2);
                edata.push_back(e);
            }
        }

        //initializing the col edges
        for (int c = 0; c <= mfgrid->ncol; c++)
        {
            for (int r = 0; r < mfgrid->nrow; r++)
            {
                VtkVertex * v1 = vdata[k * mfgrid_vsize + r * perrow_vsize + c];
                VtkVertex * v2 = vdata[k * mfgrid_vsize + (r + 1) * perrow_vsize + c];
                VtkEdge * e = new VtkEdge(v1, v2);
                edata.push_back(e);
            }
        }
    }
}
void assign_box_edge(Box * box, vector<VtkVertex *> & vdata, vector<VtkEdge *> & edata,
                     VtkBox * boxdata)
{
    int id[4];
    VtkBox & vtk_box = boxdata[box->id];
    for (int i = 0; i < 4; i++) id[i] = box->p_children[i]->id;

    ///////////////////////////////////////////////////////

    boxdata[id[0]].e[0] = vtk_box.e[0]->child_e[0];
    boxdata[id[0]].e[1] = vtk_box.e[0]->sub_e[vtk_box.e[2]];
    boxdata[id[0]].e[2] = vtk_box.e[3]->sub_e[vtk_box.e[1]];
    boxdata[id[0]].e[3] = vtk_box.e[3]->child_e[0];
    //----------------
    boxdata[id[0]].e[4] = vtk_box.e[4]->child_e[0];
    boxdata[id[0]].e[5] = vtk_box.e[4]->sub_e[vtk_box.e[6]];
    boxdata[id[0]].e[6] = vtk_box.e[7]->sub_e[vtk_box.e[5]];
    boxdata[id[0]].e[7] = vtk_box.e[7]->child_e[0];

    ///////////////////////////////////////////////////////

    boxdata[id[1]].e[0] = vtk_box.e[3]->sub_e[vtk_box.e[1]];
    boxdata[id[1]].e[1] = vtk_box.e[2]->sub_e[vtk_box.e[0]];
    boxdata[id[1]].e[2] = vtk_box.e[2]->child_e[0];
    boxdata[id[1]].e[3] = vtk_box.e[3]->child_e[1];
    //----
    boxdata[id[1]].e[4] = vtk_box.e[7]->sub_e[vtk_box.e[5]];
    boxdata[id[1]].e[5] = vtk_box.e[6]->sub_e[vtk_box.e[4]];
    boxdata[id[1]].e[6] = vtk_box.e[6]->child_e[0];
    boxdata[id[1]].e[7] = vtk_box.e[7]->child_e[1];

    ///////////////////////////////////////////////////////

    boxdata[id[2]].e[0] = vtk_box.e[1]->sub_e[vtk_box.e[3]];
    boxdata[id[2]].e[1] = vtk_box.e[1]->child_e[1];
    boxdata[id[2]].e[2] = vtk_box.e[2]->child_e[1];
    boxdata[id[2]].e[3] = vtk_box.e[2]->sub_e[vtk_box.e[0]];
    //------
    boxdata[id[2]].e[4] = vtk_box.e[5]->sub_e[vtk_box.e[7]];
    boxdata[id[2]].e[5] = vtk_box.e[5]->child_e[1];
    boxdata[id[2]].e[6] = vtk_box.e[6]->child_e[1];
    boxdata[id[2]].e[7] = vtk_box.e[6]->sub_e[vtk_box.e[4]];

    ///////////////////////////////////////////////////////

    boxdata[id[3]].e[0] = vtk_box.e[0]->child_e[1];
    boxdata[id[3]].e[1] = vtk_box.e[1]->child_e[0];
    boxdata[id[3]].e[2] = vtk_box.e[1]->sub_e[vtk_box.e[3]];
    boxdata[id[3]].e[3] = vtk_box.e[0]->sub_e[vtk_box.e[2]];
    //-------------
    boxdata[id[3]].e[4] = vtk_box.e[4]->child_e[1];
    boxdata[id[3]].e[5] = vtk_box.e[5]->child_e[0];
    boxdata[id[3]].e[6] = vtk_box.e[5]->sub_e[vtk_box.e[7]];
    boxdata[id[3]].e[7] = vtk_box.e[4]->sub_e[vtk_box.e[6]];
}
void dfs_build_box(Box * box, vector<VtkVertex *> & vdata, vector<VtkEdge *> & edata,
                   VtkBox * boxdata)
{
    if (box->is_leaf) return;

    VtkBox & vtk_box = boxdata[box->id];
    for (int k = 0; k < 2; k++)
    {
        for (int t = 0; t < 4; t++)
        {
            VtkEdge * edge = boxdata[box->id].e[k * 4 + t];
            if (edge->child_v == nullptr)
            {
                VtkVertex * nv = new VtkVertex();
                nv->pos[0] = (edge->v[0]->pos[0] + edge->v[1]->pos[0]) / 2;
                nv->pos[1] = (edge->v[0]->pos[1] + edge->v[1]->pos[1]) / 2;
                vdata.push_back(nv);
                edge->child_v = nv;

                VtkEdge * e1 = new VtkEdge(edge->v[0], nv);
                VtkEdge * e2 = new VtkEdge(nv, edge->v[1]);
                edge->child_e[0] = e1;
                edge->child_e[1] = e2;
                edata.push_back(e1);
                edata.push_back(e2);
            }
        }
    }

    // determine all the edges
    if (vtk_box.e[0]->sub_e.find(boxdata[box->id].e[2]) == vtk_box.e[0]->sub_e.end())
    {
        //center v
        VtkVertex * center_vertex = new VtkVertex();
        vdata.push_back(center_vertex);
        center_vertex->pos[0] = (vtk_box.e[0]->v[0]->pos[0] + vtk_box.e[2]->v[1]->pos[0]) / 2;
        center_vertex->pos[1] = (vtk_box.e[0]->v[0]->pos[1] + vtk_box.e[2]->v[1]->pos[1]) / 2;

        VtkEdge * e02 = new VtkEdge(vtk_box.e[0]->child_v, center_vertex);
        VtkEdge * e20 = new VtkEdge(center_vertex, vtk_box.e[2]->child_v);
        VtkEdge * e31 = new VtkEdge(vtk_box.e[3]->child_v, center_vertex);
        VtkEdge * e13 = new VtkEdge(center_vertex, vtk_box.e[1]->child_v);

        edata.push_back(e02);
        edata.push_back(e20);
        edata.push_back(e31);
        edata.push_back(e13);

        vtk_box.e[0]->sub_e.insert(make_pair(vtk_box.e[2], e02));
        vtk_box.e[2]->sub_e.insert(make_pair(vtk_box.e[0], e20));
        vtk_box.e[3]->sub_e.insert(make_pair(vtk_box.e[1], e31));
        vtk_box.e[1]->sub_e.insert(make_pair(vtk_box.e[3], e13));
    }

    if (vtk_box.e[4]->sub_e.find(vtk_box.e[6]) == vtk_box.e[4]->sub_e.end())
    {
        //center v
        VtkVertex * cv2 = new VtkVertex();
        vdata.push_back(cv2);
        cv2->pos[0] = (vtk_box.e[4]->v[0]->pos[0] + vtk_box.e[6]->v[1]->pos[0]) / 2;
        cv2->pos[1] = (vtk_box.e[4]->v[0]->pos[1] + vtk_box.e[6]->v[1]->pos[1]) / 2;

        VtkEdge * e02 = new VtkEdge(vtk_box.e[4]->child_v, cv2);
        VtkEdge * e20 = new VtkEdge(cv2, vtk_box.e[6]->child_v);
        VtkEdge * e31 = new VtkEdge(vtk_box.e[7]->child_v, cv2);
        VtkEdge * e13 = new VtkEdge(cv2, vtk_box.e[5]->child_v);

        edata.push_back(e02);
        edata.push_back(e20);
        edata.push_back(e31);
        edata.push_back(e13);

        vtk_box.e[4]->sub_e.insert(make_pair(vtk_box.e[6], e02));
        vtk_box.e[6]->sub_e.insert(make_pair(vtk_box.e[4], e20));
        vtk_box.e[7]->sub_e.insert(make_pair(vtk_box.e[5], e31));
        vtk_box.e[5]->sub_e.insert(make_pair(vtk_box.e[7], e13));
    }

    //assign edges
    assign_box_edge(box, vdata, edata, boxdata);

    //boxdata
    for (auto & i : box->p_children)
    {
        dfs_build_box(i, vdata, edata, boxdata);
    }
}
void build_box(QuadTree3D * qtree, vector<VtkVertex *> & vdata, vector<VtkEdge *> & edata,
               VtkBox * boxdata)
{
    ModflowGrid * mfgrid = qtree->get_modflow_grid();
    /////////////////////////////////////////////////////////////////////////////////////////
    //initialize the root box
    int mfgrid_esize = mfgrid->ncol * (mfgrid->nrow + 1) + mfgrid->nrow * (mfgrid->ncol + 1);
    for (int k = 0; k < mfgrid->nlay; k++)
    {
        for (int r = 0; r < mfgrid->nrow; r++)
        {
            for (int c = 0; c < mfgrid->ncol; c++)
            {
                int bid = mfgrid->get_box_index(k, r, c);
                Box * box = qtree->node_group[bid];

                for (int t = 0; t < 2; t++)
                {
                    int baseid = (k + t) * mfgrid_esize;
                    int id1 = baseid + (mfgrid->nrow + 1) * mfgrid->ncol + c * mfgrid->nrow + r;
                    int id2 = baseid + (r + 1) * mfgrid->ncol + c;
                    int id3 =
                        baseid + (mfgrid->nrow + 1) * mfgrid->ncol + (c + 1) * mfgrid->nrow + r;
                    int id4 = baseid + r * mfgrid->ncol + c;

                    boxdata[box->id].e[t * 4] = edata[id1];
                    boxdata[box->id].e[t * 4 + 1] = edata[id2];
                    boxdata[box->id].e[t * 4 + 2] = edata[id3];
                    boxdata[box->id].e[t * 4 + 3] = edata[id4];
                }
            }
        }
    }

    /////////////////////////////////////////////////////////////////////////////////////////
    // DFS build boxes
    for (int k = 0; k < mfgrid->nlay; k++)
    {
        for (int r = 0; r < mfgrid->nrow; r++)
        {
            for (int c = 0; c < mfgrid->ncol; c++)
            {
                int bid = mfgrid->get_box_index(k, r, c);
                Box * box = qtree->node_group[bid];

                dfs_build_box(box, vdata, edata, boxdata);
            }
        }
    }
}
void build_sharing(QuadTree3D * qtree, vector<VtkVertex *> & vdata, vector<VtkEdge *> & edata,
                   VtkBox * boxdata)
{

    // build VE
    build_ve(qtree, vdata, edata, boxdata);

    //build Box
    build_box(qtree, vdata, edata, boxdata);
}

//
// analysis
//

void determine_vto_box(QuadTree3D * qtree, vector<VtkVertex *> & vdata, vector<VtkEdge *> & edata,
                       VtkBox * boxdata)
{
    int nboxes = qtree->node_group.size();
    for (int i = 0; i < nboxes; i++)
    {
        VtkBox & vtk_box = boxdata[i];
        Box * box = qtree->node_group[i];
        if (!(box->is_leaf)) continue;

        vtk_box.e[0]->v[0]->box[6] = box;
        vtk_box.e[0]->v[1]->box[7] = box;
        vtk_box.e[2]->v[0]->box[5] = box;
        vtk_box.e[2]->v[1]->box[4] = box;

        vtk_box.e[4]->v[0]->box[2] = box;
        vtk_box.e[4]->v[1]->box[3] = box;
        vtk_box.e[6]->v[0]->box[1] = box;
        vtk_box.e[6]->v[1]->box[0] = box;

        //child
        if (vtk_box.e[0]->child_v != nullptr)
        {
            vtk_box.e[0]->child_v->box[6] = box;
            vtk_box.e[0]->child_v->box[7] = box;
        }
        if (vtk_box.e[1]->child_v != nullptr)
        {
            vtk_box.e[1]->child_v->box[4] = box;
            vtk_box.e[1]->child_v->box[7] = box;
        }
        if (vtk_box.e[2]->child_v != nullptr)
        {
            vtk_box.e[2]->child_v->box[4] = box;
            vtk_box.e[2]->child_v->box[5] = box;
        }
        if (vtk_box.e[3]->child_v != nullptr)
        {
            vtk_box.e[3]->child_v->box[5] = box;
            vtk_box.e[3]->child_v->box[6] = box;
        }

        //child
        if (vtk_box.e[4]->child_v != nullptr)
        {
            vtk_box.e[4]->child_v->box[2] = box;
            vtk_box.e[4]->child_v->box[3] = box;
        }
        if (vtk_box.e[5]->child_v != nullptr)
        {
            vtk_box.e[5]->child_v->box[0] = box;
            vtk_box.e[5]->child_v->box[3] = box;
        }
        if (vtk_box.e[6]->child_v != nullptr)
        {
            vtk_box.e[6]->child_v->box[1] = box;
            vtk_box.e[6]->child_v->box[2] = box;
        }
        if (vtk_box.e[7]->child_v != nullptr)
        {
            vtk_box.e[7]->child_v->box[1] = box;
            vtk_box.e[7]->child_v->box[2] = box;
        }
    }
}
void number_share_vtk_v(QuadTree3D * qtree, vector<VtkVertex *> & vdata, VtkBox * boxdata,
                        vector<VtkVertex *> & vtk_vertices)
{
    int next_number = 0;
    int box_count = qtree->node_group.size();

    for (int i = 0; i < box_count; i++)
    {
        Box * b = qtree->node_group[i];
        if (!(b->is_leaf && b->active)) continue;
        VtkBox & vtk_box = boxdata[b->id];
        for (int k = 0; k < 2; k++)
        {
            if (vtk_box.e[0 + 4 * k]->v[0]->number == -1)
            {
                vtk_box.e[0 + 4 * k]->v[0]->number = next_number++;
                vtk_vertices.push_back(vtk_box.e[0 + 4 * k]->v[0]);
            }
            if (vtk_box.e[0 + 4 * k]->v[1]->number == -1)
            {
                vtk_box.e[0 + 4 * k]->v[1]->number = next_number++;
                vtk_vertices.push_back(vtk_box.e[0 + 4 * k]->v[1]);
            }
            if (vtk_box.e[2 + 4 * k]->v[0]->number == -1)
            {
                vtk_box.e[2 + 4 * k]->v[0]->number = next_number++;
                vtk_vertices.push_back(vtk_box.e[2 + 4 * k]->v[0]);
            }
            if (vtk_box.e[2 + 4 * k]->v[1]->number == -1)
            {
                vtk_box.e[2 + 4 * k]->v[1]->number = next_number++;
                vtk_vertices.push_back(vtk_box.e[2 + 4 * k]->v[1]);
            }
        }
    }
}

//assume at most 1 depth difference
double smooth_interpolate(vector<double> & cell_z, vector<int> & depth, vector<Box *> & boxes,
                          double px, double py)
{
    set<int> depths;
    int cell_count = depth.size();
    for (int i = 0; i < cell_count; i++)
    {
        depths.insert(depth[i]);
    }
    //average
    if (depths.size() == 1 || cell_count == 1)
    {
        double z_sum(0);
        for (int i = 0; i < cell_count; i++) z_sum += cell_z[i];
        return z_sum / cell_count;
    }

    //interpolation
    if (cell_count == 4)
    {
        double d[2] = {0};
        for (int k = 0; k < 2; k++)
        {
            if (depth[k] == depth[k + 2])
                d[k] = (cell_z[k] + cell_z[k + 2]) / 2;
            else if (depth[k] < depth[k + 2])
                d[k] = (cell_z[k] + cell_z[k + 2] * 2) / 3;
            else if (depth[k] > depth[k + 2])
                d[k] = (cell_z[k] * 2 + cell_z[k + 2]) / 3;
            else
                assert(false);
        }
        return (d[0] + d[1]) / 2;
    }

    //do a plane fitting
    if (cell_count == 3)
    {
        double rows[3][3] = {};
        double * points[3] = {rows[0], rows[1], rows[2]};

        for (int i = 0; i < 3; i++)
        {
            points[i][0] = boxes[i]->x;
            points[i][1] = boxes[i]->y;
            points[i][2] = cell_z[i];
        }

        double z = fit_z(points, px, py);

        return z;
    }

    if (cell_count == 2)
    {
        if (depth[0] == depth[1]) return (cell_z[0] + cell_z[1]) / 2;
        if (depth[0] < depth[1]) return (cell_z[0] + 2 * cell_z[1]) / 3;
        if (depth[0] > depth[1]) return (cell_z[0] * 2 + cell_z[1]) / 3;
    }
    assert(false);
    return 0;
}

void interpolate_z(QuadTree3D * qtree, vector<VtkVertex *> & vdata, vector<VtkEdge *> & edata,
                   VtkBox * boxdata, vector<VtkVertex *> & vtk_vertices)
{
    double * ztop = qtree->top;
    double * zbot = qtree->bot;
    for (auto v : vtk_vertices)
    {
        //cells in the order found, the sets only skip duplicates
        set<Box *> upper_seen, lower_seen;
        vector<Box *> upper_boxes, lower_boxes;
        vector<double> upper_bottoms, lower_tops;
        vector<int> upper_depths, lower_depths;
        for (int i = 0; i < 4; i++)
        {
            if (v->box[i] != nullptr && v->box[i]->is_leaf && v->box[i]->active)
            {
                if (upper_seen.find(v->box[i]) == upper_seen.end())
                {
                    upper_seen.insert(v->box[i]);
                    upper_boxes.push_back(v->box[i]);
                    upper_bottoms.push_back(zbot[v->box[i]->number]);
                    upper_depths.push_back(v->box[i]->depth);
                }
            }
            if (v->box[i + 4] != nullptr && v->box[i + 4]->is_leaf && v->box[i + 4]->active)
            {
                if (lower_seen.find(v->box[i + 4]) == lower_seen.end())
                {
                    lower_seen.insert(v->box[i + 4]);
                    lower_boxes.push_back(v->box[i + 4]);
                    lower_tops.push_back(ztop[v->box[i + 4]->number]);
                    lower_depths.push_back(v->box[i + 4]->depth);
                }
            }
        }

        double z_estimates[2] = {0};
        //1
        if (upper_seen.size() > 0)
        {
            z_estimates[0] =
                smooth_interpolate(upper_bottoms, upper_depths, upper_boxes, v->pos[0], v->pos[1]);
        }
        //2
        if (lower_seen.size() > 0)
        {
            z_estimates[1] =
                smooth_interpolate(lower_tops, lower_depths, lower_boxes, v->pos[0], v->pos[1]);
        }

        if (upper_seen.size() > 0 && lower_seen.size() > 0)
            v->pos[2] = (z_estimates[0] + z_estimates[1]) / 2;
        else if (upper_seen.size() > 0)
            v->pos[2] = z_estimates[0];
        else
            v->pos[2] = z_estimates[1];
    }
}

//
// write
//
void write_cells(ofstream & ofile, vector<Box *> & all_boxes, vector<int> * int_values,
                 vector<double> * double_values, VtkBox * boxdata)
{
    /***************************************************************************************/
    //write the cells
    ofile << "<Cells>\n";
    //write the connectivity
    ofile << "<DataArray type=\"Int32\" Name=\"connectivity\" format=\"ascii\">\n";

    int rid = 0;
    VtkVertex * vids[8];
    int cell_num = all_boxes.size();
    for (auto b : all_boxes)
    {
        VtkBox & vtk_box = boxdata[b->id];
        vids[0] = vtk_box.e[0]->v[0];
        vids[1] = vtk_box.e[0]->v[1];
        vids[2] = vtk_box.e[1]->v[1];
        vids[3] = vtk_box.e[2]->v[0];
        vids[4] = vtk_box.e[4]->v[0];
        vids[5] = vtk_box.e[4]->v[1];
        vids[6] = vtk_box.e[5]->v[1];
        vids[7] = vtk_box.e[6]->v[0];

        ofile << vids[0]->number << "\t" << vids[1]->number << "\t" << vids[3]->number << "\t"
              << vids[2]->number << "\t" << vids[4]->number << "\t" << vids[5]->number << "\t"
              << vids[7]->number << "\t" << vids[6]->number << "\n";

        rid++;
        if (rid % 5000 == 0) ofile.flush();
    }
    ofile << "</DataArray>\n";

    //write the offset
    int sid, id = 0;
    sid = 8;
    ofile << "<DataArray type=\"Int32\" Name=\"offsets\" format=\"ascii\">\n";
    for (int i = 0; i < cell_num; i++)
    {
        ofile << sid << "\t";
        sid += 8;

        id++;
        if (id % 5000 == 0) ofile.flush();
    }
    ofile << "\n";
    ofile << "</DataArray>\n";
    //write the types
    id = 0;
    ofile << "<DataArray type=\"UInt8\" Name=\"types\" format=\"ascii\">\n";
    for (int i = 0; i < cell_num; i++)
    {
        ofile << "11\t";
        id++;
        if (id % 5000 == 0) ofile.flush();
    }
    ofile << "</DataArray>\n";
    ofile << "</Cells>\n";
}
void write_vtu_file(QuadTree3D * qtree, vector<VtkVertex *> & vdata, vector<VtkEdge *> & edata,
                    VtkBox * boxdata, vector<VtkVertex *> & vtk_vertices, string vtk_name)
{
    ModflowGrid * mfgrid = qtree->get_modflow_grid();
    //begin to write
    double cx, cy, rotation;
    mfgrid->get_rotation(cx, cy, rotation);

    vector<Box *> all_boxes;

    vector<string> int_names, double_names;
    vector<int> int_values[6]; // box_index, cells, layers, rows, cols, vtkid;
    int_names.push_back("nodeid");
    int_names.push_back("nodenumber");
    int_names.push_back("layer");
    int_names.push_back("row");
    int_names.push_back("col");
    int_names.push_back("vtkid");

    vector<double> double_values[4]; //tops, bottoms, delrs, delcs;
    double_names.push_back("top");
    double_names.push_back("bottom");
    double_names.push_back("delr");
    double_names.push_back("delc");
    //
    collect_cell_share(qtree, all_boxes, int_values, double_values,
                       qtree->get_one_based_node_numbering());

    //world coordinates of each point from the nearest corner of a cell it bounds
    GridVertices gv(qtree);
    map<VtkVertex *, Point2d> world;
    for (auto & box : all_boxes)
    {
        VtkBox & vtk_box = boxdata[box->id];
        VtkVertex * vids[8] = {vtk_box.e[0]->v[0], vtk_box.e[0]->v[1], vtk_box.e[1]->v[1],
                               vtk_box.e[2]->v[0], vtk_box.e[4]->v[0], vtk_box.e[4]->v[1],
                               vtk_box.e[5]->v[1], vtk_box.e[6]->v[0]};
        const int * ids = gv.corners(box);
        for (auto & vid : vids)
        {
            if (world.count(vid)) continue;
            double rotated_x, rotated_y, dmin = DBL_MAX;
            rotate_one_pnt(vid->pos[0], vid->pos[1], cx, cy, rotation, rotated_x, rotated_y);
            for (int j = 0; j < 4; j++)
            {
                const Point2d & p = gv.vertices()[ids[j]];
                double d = (p[0] - rotated_x) * (p[0] - rotated_x) +
                           (p[1] - rotated_y) * (p[1] - rotated_y);
                if (d < dmin)
                {
                    dmin = d;
                    world[vid] = p;
                }
            }
        }
    }

    //write the points, loop over each layer
    stringstream ss0;
    ss0 << vtk_name.c_str();
    ss0 << ".vtu";
    ofstream ofile(ss0.str().c_str());
    write_vtk_header(ofile);

    /***************************************************************************************/
    ofile << "<Piece NumberOfPoints=\"" << vtk_vertices.size() << "\" NumberOfCells=\""
          << all_boxes.size() << "\" >";

    //write the points
    ofile << "<Points>\n";
    ofile << "<DataArray type=\"Float64\" NumberOfComponents=\"3\" format=\"ascii\">\n";
    int flash = 0;
    for (auto sv : vtk_vertices)
    {
        //points and cells both come from the active leaves
        assert(world.count(sv));
        const Point2d & xy = world[sv];

        ofile << RoundTrip(xy[0]) << "\t" << RoundTrip(xy[1]) << "\t"
              << RoundTrip(sv->pos[2] * ZSCALE) << "\t";

        if (flash % 5000 == 0)
        {
            ofile.flush();
        }
        flash++;
    }

    ofile << "\n";
    ofile << "</DataArray>\n";
    ofile << "</Points>\n";

    //write the boxes
    write_cells(ofile, all_boxes, int_values, double_values, boxdata);

    //write the attributes
    ofile << "<CellData Scalars=\"scalars\">\n";
    int int_size = int_names.size();
    for (int i = 0; i < int_size; i++)
    {
        write_data_array(ofile, "Int32", int_names[i], int_values[i]);
    }

    for (int i = 0; i < 4; i++)
    {
        write_data_array(ofile, "Float64", double_names[i], (double_values[i]));
    }
    ofile << "</CellData>\n";

    ofile << "</Piece>\n";
    write_vtk_tail(ofile);

    ofile.close();
}

//
static void write_vtk_shared_vertices(QuadTree3D * qtree, string vtk_name)
{
    qtree->number_nodes();

    int nboxes = qtree->node_group.size();

    vector<VtkVertex *> vdata;
    vector<VtkEdge *> edata;
    vector<VtkBox> boxes(nboxes);
    VtkBox * boxdata = boxes.data();

    //build the sharing v
    build_sharing(qtree, vdata, edata, boxdata);
    //determine the neighbor boxes for v
    determine_vto_box(qtree, vdata, edata, boxdata);

    //number the vertices
    vector<VtkVertex *> vtk_vertices;
    number_share_vtk_v(qtree, vdata, boxdata, vtk_vertices);

    //export to VTK file
    interpolate_z(qtree, vdata, edata, boxdata, vtk_vertices);

    write_vtu_file(qtree, vdata, edata, boxdata, vtk_vertices, vtk_name);

    //clear memory
    for (auto & vit : vdata)
    {
        delete vit;
    }
    for (auto & eit : edata)
    {
        delete eit;
    }
}
