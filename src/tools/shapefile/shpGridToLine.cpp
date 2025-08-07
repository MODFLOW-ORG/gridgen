#include "shpGridToLine.h"
using namespace cusg;

void  rotate(double cx, double cy, double rotation, double sx, double sy, double& ex, double & ey)
{
	ex = cx + (sx - cx) * cos(rotation) - (sy - cy)*sin(rotation);
	ey = cy + (sy - cy)* cos(rotation) + (sx - cx)*sin(rotation);

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

void writeOneMFGrid(cusg::ModflowGrid* mfgrid, SHPHandle& handler, DBFHandle& dbfHandle, int lyrId)
{
	double cx, cy, rotation;
	mfgrid->getRotatePara(cx, cy, rotation);

	double tmpX[2], tmpY[2];
	double resX[2], resY[2];

	for(int k = 0; k < mfgrid->nlay; k++)
	{
		//write columns
		double * Xe = mfgrid->get_local_Xe_array();
		double * Ye = mfgrid->get_local_Ye_array();
		for(int i = 0; i <= mfgrid->ncol; i++)
		{
			tmpX[0] = mfgrid->Xe[i];
			tmpX[1] = mfgrid->Xe[i];
			tmpY[0] = mfgrid->Ye[0];
			tmpY[1] = mfgrid->Ye[mfgrid->nrow];

			rotate(cx, cy, rotation, tmpX[0], tmpY[0], resX[0], resY[0]);
			rotate(cx, cy, rotation, tmpX[1], tmpY[1], resX[1], resY[1]);

			writeOneLine(2, resX, resY, lyrId, k, handler, dbfHandle);
		}

		for(int i = 0; i <= mfgrid->nrow; i++)
		{
				tmpY[0]=  mfgrid->Ye[i];
				tmpY[1] = mfgrid->Ye[i];
				tmpX[0] = mfgrid->Xe[0];
				tmpX[1] = mfgrid->Xe[mfgrid->ncol];

				rotate(cx, cy, rotation, tmpX[0], tmpY[0], resX[0], resY[0]);
				rotate(cx, cy, rotation, tmpX[1], tmpY[1], resX[1], resY[1]);

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

	double cx, cy, rotation;
	mfgrid->getRotatePara(cx, cy, rotation);

	writeOneMFGrid(mfgrid, handler, dbfHandle, lyrId);

	double tmpX[2], tmpY[2];
	double resX[2], resY[2];
	SHPObject* shpObj = NULL;
	//get through each cell
	for(NodeGroup::iterator nit = qtree->nodegroup.begin(); nit != qtree->nodegroup.end(); ++nit)
	{
		Box* b = *nit;
		if(b->isLeaf == false)
		{
			tmpX[0] = b->x;
			tmpX[1] = b->x;
			tmpY[0] = b->y - b->dy/2;
			tmpY[1] = b->y + b->dy/2;

			rotate(cx, cy, rotation, tmpX[0], tmpY[0], resX[0], resY[0]);
			rotate(cx, cy, rotation, tmpX[1], tmpY[1], resX[1], resY[1]);

			writeOneLine(2, resX, resY, lyrId, b->layer, handler, dbfHandle);

			////////////////////////////////////////
			tmpX[0] = b->x - b->dx/2;
			tmpX[1] = b->x + b->dx/2;
			tmpY[0] = b->y;
			tmpY[1] = b->y;

			rotate(cx, cy, rotation, tmpX[0], tmpY[0], resX[0], resY[0]);
			rotate(cx, cy, rotation, tmpX[1], tmpY[1], resX[1], resY[1]);

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

	writeOneMFGrid(mfgrid, handler, dbfHandle, lyrId);

	DBFClose(dbfHandle);
	SHPClose(handler);
}