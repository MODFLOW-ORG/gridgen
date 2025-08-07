#include "QuadTree3D.h"
#include "shapelib/shapefil.h"

//class Grid2LineWriter
//{
//private:
//	SHPHandle handler;
//	DBFHandle dbfHandle;
//
//};


void writeQuadtree2ShpLine(cusg::QuadTree3D* qtree, string gridShpName);

void writeModflowGrid2ShpLine(cusg::ModflowGrid* modflow, string gridShpName);