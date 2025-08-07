
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
		for(int i = 0; i < 8; i++)
			vertices[i] = NULL;
	}
	Vertex* vertices[8];
};

//export modflow grid to a paraview vtk file
void writeModflowGrid2VTK(ModflowGrid* modflow, string vtkName, bool without_inactive = true, bool one_based_numbering = true);

void writeModflowGrid2VTK_shareV(ModflowGrid* mfgrid, string vtkName, bool without_inactive = true, bool one_based_numbering = true);

//write the Quadtree grid
void writeQuadtreeGrid2VTK(QuadTree3D* qtree, string vtkName);

void writeQuadtreeGrid2VTK_shareV(QuadTree3D* qtree, string vtkName);


//build sharing vertices data here
void buildSharingVertices(QuadTree3D* qtree, vector<Vertex*>& vdata, BoxVtx* boxVts);

void rotateOnePnt(double x, double y, double cx, double cy, double rotation,double& resx, double& resy);


/////////////////////////////////////////////////////////////////////////////////////////////////////////////
//util

//collect share
int collectCell_share(QuadTree3D* qtree,  vector<Box*>& allBoxes, vector<int>* intvals, vector<float>* dbvals, bool one_based_numbering);

void writeVTKHeader(ofstream& ofile);

void writeVTKTail(ofstream& ofile);

void writeVTKCell(ofstream& ofile, int boxNum);

template<typename Dtype>
void writeDataArray(ofstream& ofile, string type, string name, const vector<Dtype>& data);

void writeCells(ofstream& ofile,  vector<Box*>& allBoxes, vector<int>* intvals, vector<float>* dbvals, BoxVtx* boxVts);
//fit a plane using 3 points
double fitZ(double** vv, double x, double y);

#endif