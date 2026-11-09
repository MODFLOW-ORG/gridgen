
#ifndef CUSG_QUADTREE_VTK_H
#define CUSG_QUADTREE_VTK_H

#include "QuadTree3D.h"
#include "shapelib/shapefil.h"
#include "ascii_grid.h"
#include "GridVertices.h"
#include <map>
#include <string>
#include <vector>
using namespace std;
using namespace cusg;

#define ZSCALE 1

/*!\class Vertex
*  Store the information of vertex
*/
class Vertex
{
public:
    Vertex(double p1, double p2, double p3, int d)
    {
        pos.set(p1, p2, p3);
        id = d;
    }
    Point3d pos;
    int id;
};
class BoxVtx
{
public:
    BoxVtx()
    {
        for (auto & v : vertices) v = nullptr;
    }
    Vertex * vertices[8];
};

//write a MODFLOW grid to a VTK file, with the cells sharing their vertices or
//each cell with its own
void write_vtk(ModflowGrid * modflow, string vtk_name, bool share_vertices,
               bool without_inactive = true, bool one_based_numbering = true);

//write a quadtree grid to a VTK file, with the cells sharing their vertices or
//each cell with its own
void write_vtk(QuadTree3D * qtree, string vtk_name, bool share_vertices);

//elevation of each grid vertex at each layer interface (0 is the top of layer 1)
//used by an active cell, as the VTK file with shared vertices gives it
map<pair<int, int>, double> vertex_elevations(Grid * grid, const GridVertices & gv);

//build sharing vertices data here
void build_sharing_vertices(QuadTree3D * qtree, vector<Vertex *> & vdata, BoxVtx * box_vertices);

void rotate_one_pnt(double x, double y, double cx, double cy, double rotation, double & rotated_x,
                    double & rotated_y);

/////////////////////////////////////////////////////////////////////////////////////////////////////////////
//util

//collect share
int collect_cell_share(QuadTree3D * qtree, vector<Box *> & all_boxes, vector<int> * int_values,
                       vector<double> * double_values, bool one_based_numbering);

void write_vtk_header(ofstream & ofile);

void write_vtk_tail(ofstream & ofile);

void write_vtk_cell(ofstream & ofile, int nboxes);

template <typename Dtype>
void write_data_array(ofstream & ofile, string type, string name, const vector<Dtype> & data);

void write_cells(ofstream & ofile, vector<Box *> & all_boxes, vector<int> * int_values,
                 vector<double> * double_values, BoxVtx * box_vertices);
//fit a plane using 3 points
double fit_z(double ** points, double x, double y);

#endif