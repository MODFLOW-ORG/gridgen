
#ifndef CUSG_QUADTREE_VTK_H
#define CUSG_QUADTREE_VTK_H

#include "QuadTree3D.h"
#include "shapelib/shapefil.h"
#include "ascii_grid.h"
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
void write_vtk(ModflowGrid * modflow, string vtkName, bool share_vertices,
               bool without_inactive = true, bool one_based_numbering = true);

//write a quadtree grid to a VTK file, with the cells sharing their vertices or
//each cell with its own
void write_vtk(QuadTree3D * qtree, string vtkName, bool share_vertices);

//build sharing vertices data here
void buildSharingVertices(QuadTree3D * qtree, vector<Vertex *> & vdata, BoxVtx * box_vertices);

void rotateOnePnt(double x, double y, double cx, double cy, double rotation, double & rotated_x,
                  double & rotated_y);

/////////////////////////////////////////////////////////////////////////////////////////////////////////////
//util

//collect share
int collectCell_share(QuadTree3D * qtree, vector<Box *> & allBoxes, vector<int> * int_values,
                      vector<double> * double_values, bool one_based_numbering);

void writeVTKHeader(ofstream & ofile);

void writeVTKTail(ofstream & ofile);

void writeVTKCell(ofstream & ofile, int nboxes);

template <typename Dtype>
void writeDataArray(ofstream & ofile, string type, string name, const vector<Dtype> & data);

void writeCells(ofstream & ofile, vector<Box *> & allBoxes, vector<int> * int_values,
                vector<double> * double_values, BoxVtx * box_vertices);
//fit a plane using 3 points
double fitZ(double ** points, double x, double y);

#endif