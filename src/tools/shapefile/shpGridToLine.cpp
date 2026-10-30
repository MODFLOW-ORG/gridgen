#include "shpGridToLine.h"
#include "GridVertices.h"
using namespace cusg;

//coordinates of a corner of a box, numbered clockwise from the upper left
void cornerPoint(const GridVertices& gv, const Box* b, int corner, double& x, double& y)
{
	const Point2d& p = gv.vertices()[gv.corners(b)[corner]];
	x = p[0];
	y = p[1];
}
void prepareWriteLine(string gridShpName, SHPHandle* shd, DBFHandle* dbf, int* lyrId)
{
	*shd = SHPCreate(gridShpName.c_str(), SHPT_ARC);
	if((*shd) == NULL)
	{
		cerr<<"can't create shape file: "<<gridShpName.c_str()<<endl;
		return;
	}

		//create a .dbf file
	string dbfName = gridShpName;/*
	int lastL = gridShpName.find_last_of('/') == string::npos ? gridShpName.find;*/
	if( gridShpName.find('.')!=string::npos )
	{
		int dotIndex = gridShpName.find_last_of('.');
		dbfName = gridShpName.substr(0,dotIndex);
	}
	dbfName += ".dbf";
	*dbf = DBFCreate(dbfName.c_str());

	*lyrId = DBFAddField(*dbf, "layer", FTInteger, 20, 0);
}



void writeOneLine(int vNum, double *tmpX, double* tmpY, int lyrID, int layer, SHPHandle handler, DBFHandle dbfHandle)
{
	SHPObject* shpObj = NULL;
	//write to the shape file
	shpObj = SHPCreateSimpleObject(SHPT_ARC, vNum, tmpX, tmpY, NULL);
	int entityNum = SHPWriteObject(handler, -1, shpObj);
	//clear memory
	SHPDestroyObject(shpObj);

	//write the entity feature info
	int success = DBFWriteIntegerAttribute(dbfHandle, entityNum, lyrID, layer);
}

//write the column and row edges of the base grid; grid holds the base cells
void writeOneMFGrid(cusg::Grid* grid, cusg::ModflowGrid* mfgrid, const GridVertices& gv, SHPHandle& handler, DBFHandle& dbfHandle, int lyrId)
{
	double resX[2], resY[2];
	int nrow = mfgrid->nrow;
	int ncol = mfgrid->ncol;

	for(int k = 0; k < mfgrid->nlay; k++)
	{
		//write columns, from the top of the first row to the bottom of the last
		for(int i = 0; i <= ncol; i++)
		{
			int c = min(i, ncol - 1);
			Box* top = grid->nodegroup[mfgrid->get_nodeid(k, 0, c)];
			Box* bot = grid->nodegroup[mfgrid->get_nodeid(k, nrow - 1, c)];
			cornerPoint(gv, top, (i < ncol) ? 0 : 1, resX[0], resY[0]);
			cornerPoint(gv, bot, (i < ncol) ? 3 : 2, resX[1], resY[1]);

			writeOneLine(2, resX, resY, lyrId, k, handler, dbfHandle);
		}

		//write rows, from the west edge to the east edge
		for(int i = 0; i <= nrow; i++)
		{
				int r = min(i, nrow - 1);
				Box* west = grid->nodegroup[mfgrid->get_nodeid(k, r, 0)];
				Box* east = grid->nodegroup[mfgrid->get_nodeid(k, r, ncol - 1)];
				cornerPoint(gv, west, (i < nrow) ? 0 : 3, resX[0], resY[0]);
				cornerPoint(gv, east, (i < nrow) ? 1 : 2, resX[1], resY[1]);

				writeOneLine(2, resX, resY, lyrId, k, handler, dbfHandle);
		}

	}

}
void writeQuadtree2ShpLine(cusg::QuadTree3D* qtree, string gridShpName)
{
	SHPHandle handler;
	DBFHandle dbfHandle;
	int lyrId;
	prepareWriteLine(gridShpName, &handler, &dbfHandle, &lyrId);

	ModflowGrid* mfgrid = qtree->getModflowGrid();
	unsigned int id_offset = 0;
	unsigned int id = id_offset;
	//char extra[64];

	GridVertices gv(qtree);
	writeOneMFGrid(qtree, mfgrid, gv, handler, dbfHandle, lyrId);

	double resX[2], resY[2];
	SHPObject* shpObj = NULL;
	//get through each cell
	for(NodeGroup::iterator nit = qtree->nodegroup.begin(); nit != qtree->nodegroup.end(); ++nit)
	{
		Box* b = *nit;
		if(b->isLeaf == false)
		{
			//from the south edge to the north edge through the center
			cornerPoint(gv, b->pChildren[3], 2, resX[0], resY[0]);
			cornerPoint(gv, b->pChildren[0], 1, resX[1], resY[1]);

			writeOneLine(2, resX, resY, lyrId, b->layer, handler, dbfHandle);

			////////////////////////////////////////
			//from the west edge to the east edge through the center
			cornerPoint(gv, b->pChildren[3], 0, resX[0], resY[0]);
			cornerPoint(gv, b->pChildren[2], 1, resX[1], resY[1]);

			writeOneLine(2, resX, resY, lyrId, b->layer, handler, dbfHandle);
		}

	}

	DBFClose(dbfHandle);
	SHPClose(handler);
}

void writeModflowGrid2ShpLine(cusg::ModflowGrid* mfgrid, string gridShpName)
{
	SHPHandle handler;
	DBFHandle dbfHandle;
	int lyrId;
	prepareWriteLine(gridShpName, &handler, &dbfHandle, &lyrId);

	GridVertices gv(mfgrid);
	writeOneMFGrid(mfgrid, mfgrid, gv, handler, dbfHandle, lyrId);

	DBFClose(dbfHandle);
	SHPClose(handler);
}