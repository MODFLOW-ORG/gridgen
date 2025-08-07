
#include "GridToVTK.h"
#include <fstream>
#include <sstream>
using namespace std;

//int debug_num[4] = {0};
static int ShareVID = 0;

//reverse interpolate the z coordinates for all boxes
//void reverseInterpolate(QuadTree3D* qtree);

inline void rotateOneBox(Box* b, double cx, double cy, double rotation, double* padfX, double* padfY, double* padfZ, double* padfM, double* ctrX, double* ctrY)
{
	double ldX, ldY, luX, luY, rdX, rdY, ruX, ruY, ctrx, ctry;
	b->rotate(cx, cy, rotation, ldX, ldY, luX, luY, rdX, rdY, ruX, ruY);
	b->rotateCtr(cx, cy, rotation, ctrx, ctry);

	padfX[0] = luX; padfX[1] = ruX;
	padfX[2] = rdX; padfX[3] = ldX;
	padfY[0] = luY; padfY[1] = ruY;
	padfY[2] = rdY; padfY[3] = ldY;
	ctrX[0] = ctrx; ctrY[0] = ctry;
}

/********************************************************/
//write the vtk file
void writeVTKHeader(ofstream& ofile)
{
	ofile<<"<?xml version=\"1.0\"?>\n";
	ofile<<"<VTKFile type=\"UnstructuredGrid\" version=\"0.1\" byte_order=\"BigEndian\">\n";
	ofile<<"<UnstructuredGrid GhostLevel=\"0\">\n";
}
void writeVTKTail(ofstream& ofile)
{
	ofile<<"</UnstructuredGrid>\n";
	ofile<<"</VTKFile>\n";
}

void writeVTKCell(ofstream& ofile, int boxNum)
{
	/***************************************************************************************/
	//write the cells
	ofile<<"<Cells>\n";
	//write the connectivity
	ofile<<"<DataArray type=\"Int32\" Name=\"connectivity\" format=\"ascii\">\n";
	int sid = 0, id = 0;
	for(int i = 0; i < boxNum; i++)
	{
		ofile<<(sid)<<"\t"<<(sid+6)<<"\t"<<(sid+2)<<"\t"<<(sid+4)<<"\t"
		<<(sid+1)<<"\t"<<(sid+7)<<"\t"<<(sid+3)<<"\t"<<(sid+5)<<"\n";

			//ofile<<(sid+6)<<"\t"<<sid<<"\t"<<(sid+4)<<"\t"<<(sid+2)<<"\t";
			//	ofile<<(sid+7)<<"\t"<<(sid+1)<<"\t"<<(sid+5)<<"\t"<<(sid+3)<<"\n";
				sid+=8;

				id++;
				if(id%500==0)
					ofile.flush();
	}
	ofile<<"</DataArray>\n";
	//write the offset
	sid = id = 0;
	sid = 8;
	ofile<<"<DataArray type=\"Int32\" Name=\"offsets\" format=\"ascii\">\n";
	for(int i =0; i < boxNum; i++)
	{
				ofile<<sid<<"\t";
				sid+=8;

				id++;
				if(id%1000==0)
					ofile.flush();
	}
	ofile<<"\n";
	ofile<<"</DataArray>\n";
	//write the types
	id = 0;
	ofile<<"<DataArray type=\"UInt8\" Name=\"types\" format=\"ascii\">\n";
	for(int i = 0; i < boxNum; i++)
	{
		ofile<<"11\t";
		id++;
		if(id%1000==0)
			ofile.flush();
	}
	ofile<<"</DataArray>\n";
	ofile<<"</Cells>\n";

}

template<typename Dtype>
void writeDataArray(ofstream& ofile, string type, string name, const vector<Dtype>& data)
{
	//typename vector<Dtype>::const_iterator DIT;
	ofile<<"<DataArray type=\""<<type<<"\" Name=\""<<name<<"\" format=\"ascii\">";
	int id = 0;
	for(typename vector<Dtype>::const_iterator dit = data.begin(); dit != data.end(); dit++)
	{
		ofile<<(*dit)<<"\t";
		
		if(id % 1000==0)
			ofile.flush();
		id++;
	}
	ofile<<"</DataArray>\n";
}

//export modflow grid to a paraview vtk file
void writeModflowGrid2VTK(ModflowGrid* modflow, string vtkName, bool without_inactive, bool one_based_numbering)
{
	double cx, cy, rotation;
	modflow->getRotatePara(cx, cy, rotation);
	//write the points, loop over each layer
	int nlay = modflow->nlay;
	int nrow = modflow->nrow;
	int ncol = modflow->ncol;

	double padfX[4] = {0};
	double padfY[4] = {0};
	double padfZ[4] = {0};
	double padfM[4] = {0};
	double ctrX[1] = {0};
	double ctrY[1] = {0};

	stringstream ss0;
	ss0<<vtkName.c_str();
	ss0<<".vtu";
	ofstream ofile(ss0.str().c_str());
	writeVTKHeader(ofile);

	/***************************************************************************************/
	int totalCell = nlay * nrow * ncol;
	ofile<<"<Piece NumberOfPoints=\""<<totalCell*8<<"\" NumberOfCells=\""<<totalCell<<"\" >";

	vector<int> nodeIds, nodeNums, layers, rows, cols,vtkids;
	vector<double> tops, bottoms, delrs, delcs;

	//write the points
	ofile<<"<Points>\n";
	ofile<<"<DataArray type=\"Float64\" NumberOfComponents=\"3\" format=\"ascii\">\n";
	int id = 0, cellNum = 0;
	for(int l=0;l<nlay;l++)
	{
		for(int r=0;r<nrow;r++)
		{
			for(int c=0;c<ncol;c++)
			{
				int bid=modflow->get_nodeid(l,r,c);
				Box * box=modflow->nodegroup[bid];

				//if(box->isLeaf && box->active)
				if(box->isLeaf)
				{
					if(without_inactive && (!box->active))
						continue;

					double top= modflow->bot[ l * nrow * ncol + r * ncol + c];
					double bot= modflow->bot[(l + 1) * nrow * ncol + r * ncol + c] ;

					rotateOneBox(box, cx, cy, rotation, padfX, padfY, padfZ, padfM, ctrX, ctrY);
					for(int i = 0; i < 4; i++)
					{
						ofile<<padfX[i]<<"\t"<<padfY[i]<<"\t"<<top * ZSCALE<<"\t" ;		//(box->z - box->dz / 2)<<"\t";
						ofile<<padfX[i]<<"\t"<<padfY[i]<<"\t"<<bot * ZSCALE<<"\t";		//(box->z + box->dz / 2)<<"\t";

						id+=2;
						if(id%500==0)
							ofile.flush();
					}
					
					//store info
					if(!one_based_numbering)
					{
						nodeIds.push_back(box->id);
						nodeNums.push_back(box->number); layers.push_back(l); rows.push_back(r); cols.push_back(c);
					}
					else
					{
						nodeIds.push_back(box->id + 1);
						nodeNums.push_back(box->number + 1); layers.push_back(l + 1); rows.push_back(r + 1); cols.push_back(c + 1);
					}

					tops.push_back(top); bottoms.push_back(bot); delrs.push_back(box->dx); delcs.push_back(box->dy);
				

					vtkids.push_back(cellNum);
					cellNum++;
				}
			}
		}
	}
	ofile<<"\n";
	ofile<<"</DataArray>\n";
	ofile<<"</Points>\n";

	//write cells
	writeVTKCell(ofile, id/8);

	//write the attributes
	ofile<<"<CellData Scalars=\"scalars\">\n";
	writeDataArray(ofile, "Int32", "id", nodeIds);
	writeDataArray(ofile, "Int32", "nodenumber", nodeNums);
	writeDataArray(ofile, "Int32", "layer", layers);
	writeDataArray(ofile, "Int32", "row", rows);
	writeDataArray(ofile, "Int32", "col", cols);
	writeDataArray(ofile, "Int32", "vtkid", vtkids);

	writeDataArray(ofile, "Float64", "top", tops);
	writeDataArray(ofile, "Float64", "bottom", bottoms);
	writeDataArray(ofile, "Float64", "delr", delrs);
	writeDataArray(ofile, "Float64", "delc", delcs);
	ofile<<"</CellData>\n";

	ofile<<"</Piece>\n";
	writeVTKTail(ofile);

	ofile.close();

}

void DFSWriteQDT2VTK(QuadTree3D* qtree, ofstream& ofile,  int& pntNum, Box* b, int l,int r, int c,bool one_based_numbering,double cx, double cy, double rotation,
	vector<int>* intvals, vector<double>* dbvals)
{
	char tmp[8]={'\0'};
	//if(b->isLeaf)
	if(b->isLeaf && b->active)
	{
		double padfX[4] = {0};
		double padfY[4] = {0};
		double padfZ[4] = {0};
		double padfM[4] = {0};
		double ctrX[1] = {0};
		double ctrY[1] = {0};

		rotateOneBox(b, cx, cy, rotation, padfX, padfY, padfZ, padfM, ctrX, ctrY);

		double top=qtree->top[b->number];
		double bot=qtree->bot[b->number];

		for(int i = 0; i < 4; i++)
		{
				ofile<<padfX[i]<<"\t"<<padfY[i]<<"\t"<<top * ZSCALE<<"\t";
				ofile<<padfX[i]<<"\t"<<padfY[i]<<"\t"<<bot * ZSCALE<<"\t";

				pntNum+=2;
				if(pntNum%500==0)
					ofile.flush();
		}
		//store info

		if(!one_based_numbering)
		{
			intvals[0].push_back(b->id); 
			intvals[1].push_back(b->number); 
			intvals[2].push_back(l); 
			intvals[3].push_back(r); 
			intvals[4].push_back(c);
		}
		else
		{
			intvals[0].push_back(b->id + 1); intvals[1].push_back(b->number + 1); 
			intvals[2].push_back(l + 1); intvals[3].push_back(r + 1); intvals[4].push_back(c + 1);
		}

		dbvals[0].push_back(top); 
		dbvals[1].push_back(bot); 
		dbvals[2].push_back(b->dx); 
		dbvals[3].push_back(b->dy);

		intvals[5].push_back(pntNum);//the vtk_id
	}
	else if(b->isLeaf==false)//b is not a leaf, call this function recursively
	{
		for (int i = 0; i < 4; i++)
		{
			//if(one_based_numbering) sprintf(tmp,"%d",i+1);
			//else sprintf(tmp,"%d",i);
			DFSWriteQDT2VTK(qtree, ofile,pntNum,b->pChildren[i],l,r,c,one_based_numbering, cx, cy, rotation, intvals, dbvals);
		}
	}//end b->isLeaf
}


//write the Quadtree grid
void writeQuadtreeGrid2VTK(QuadTree3D* qtree, string vtkName)
{
	//make sure the number of the nodes are correct
	qtree->number_nodes();


	ModflowGrid * mfgrid=qtree->getModflowGrid();

	double cx, cy, rotation;
	mfgrid->getRotatePara(cx, cy, rotation);

	//write the points, loop over each layer
	stringstream ss0;
	ss0<<vtkName.c_str();
	ss0<<".vtu";
	ofstream ofile(ss0.str().c_str());
	writeVTKHeader(ofile);
	
	/***************************************************************************************/
	int totalCell = qtree->nodes();
	ofile<<"<Piece NumberOfPoints=\""<<totalCell*8<<"\" NumberOfCells=\""<<totalCell<<"\" >";

	vector<string> intNames, dbNames;
	vector<int> intvals[6];// nodeid, nodeNums, layers, rows, cols, vtkid;
	intNames.push_back("nodeid");
	intNames.push_back("nodenumber"); 
	intNames.push_back("layer"); 
	intNames.push_back("row"); 
	intNames.push_back("col");
	intNames.push_back("vtkid");

	vector<double> dbvals[4];//tops, bottoms, delrs, delcs;
	dbNames.push_back("top"); 
	dbNames.push_back("bottom"); 
	dbNames.push_back("delr"); 
	dbNames.push_back("delc");

	//write the points
	ofile<<"<Points>\n";
	ofile<<"<DataArray type=\"Float64\" NumberOfComponents=\"3\" format=\"ascii\">\n";
	int pntNum = 0;
	for(int l=0;l<mfgrid->nlay;l++)
	{
		for(int r=0;r<mfgrid->nrow;r++)
		{
			for(int c=0;c<mfgrid->ncol;c++)
			{
				int bid=mfgrid->get_nodeid(l,r,c);
				Box * box=qtree->nodegroup[bid];

				if(box->active)
					DFSWriteQDT2VTK(qtree, ofile, pntNum,  box,  l, r,  c,qtree->get_one_based_node_numbering(), cx,  cy,  rotation, intvals, dbvals);
			}
		}
	}

	ofile<<"\n";
	ofile<<"</DataArray>\n";
	ofile<<"</Points>\n";

	//write cells
	writeVTKCell(ofile, pntNum/8);

	//write the attributes
	ofile<<"<CellData Scalars=\"scalars\">\n";
	int intSize = intNames.size();
	for(int i = 0; i < intSize; i++)
	{
		writeDataArray(ofile, "Int32", intNames[i], intvals[i]);
	}

	int dbSize = dbNames.size();
	for(int i = 0; i < 4; i++)
	{
		writeDataArray(ofile, "Float64", dbNames[i], dbvals[i]);
	}
	ofile<<"</CellData>\n";

	ofile<<"</Piece>\n";
	writeVTKTail(ofile);

	ofile.close();
}

void writeCells(ofstream& ofile,  vector<Box*>& allBoxes, vector<int>* intvals, vector<float>* dbvals, BoxVtx* boxVts)
{
	//ModflowGrid* mfgrid = qtree->getModflowGrid();
	/***************************************************************************************/
	//write the cells
	ofile<<"<Cells>\n";
	//write the connectivity
	ofile<<"<DataArray type=\"Int32\" Name=\"connectivity\" format=\"ascii\">\n";

	int rid = 0;
	int cellNum = allBoxes.size();
	for(vector<Box*>::iterator bit = allBoxes.begin(); bit != allBoxes.end(); ++bit)
	{
		Box* b = *bit;
		ofile<<boxVts[b->id].vertices[0]->id<<"\t"<<boxVts[b->id].vertices[1]->id<<"\t"<<boxVts[b->id].vertices[3]->id<<"\t"<<boxVts[b->id].vertices[2]->id<<"\t"
			<<boxVts[b->id].vertices[4]->id<<"\t"<<boxVts[b->id].vertices[5]->id<<"\t"<<boxVts[b->id].vertices[7]->id<<"\t"<<boxVts[b->id].vertices[6]->id<<"\n";

		rid++;
		if(rid % 5000 == 0)
			ofile.flush();
	}
	ofile<<"</DataArray>\n";

	//write the offset
	int sid, id = 0;
	sid = 8;
	ofile<<"<DataArray type=\"Int32\" Name=\"offsets\" format=\"ascii\">\n";
	for(int i =0; i < cellNum; i++)
	{
		ofile<<sid<<"\t";
		sid+=8;

		id++;
		if(id%5000==0) ofile.flush();
	}
	ofile<<"\n";
	ofile<<"</DataArray>\n";
	//write the types
	id = 0;
	ofile<<"<DataArray type=\"UInt8\" Name=\"types\" format=\"ascii\">\n";
	for(int i = 0; i < cellNum; i++)
	{
		ofile<<"11\t";
		id++;
		if(id%5000==0)
			ofile.flush();
	}
	ofile<<"</DataArray>\n";
	ofile<<"</Cells>\n";
}
/******************************************************************************************************************
* write VTK file with sharing vertices for modflow grid
*******************************************************************************************************************/
void buildSharingVertices_mfgrid(ModflowGrid* mfgrid, vector<Vertex*>& vdata, BoxVtx* boxVts)
{
	double val;//, resx, resy;
	double cx, cy, rotation;
	mfgrid->getRotatePara(cx, cy, rotation);

	int mfgrid_csize=(mfgrid->nrow) * (mfgrid->ncol);
	
	//get all the Xe and Ye
	double* Xe = mfgrid->get_local_Xe_array();
	double* Ye = mfgrid->get_local_Ye_array();

	//use the edges to locate the box corners
	int vid = 0;
	for(int k = 0; k <= mfgrid->nlay; k++)
	{
		double* zgrid = &(mfgrid->bot[mfgrid_csize * k]);
		mf_ascii_grid ag(mfgrid, zgrid, mfgrid->ncol, mfgrid->nrow);

		for(int r=0; r<=mfgrid->nrow;r++)
		{
			for(int c=0; c<=mfgrid->ncol;c++)
			{
				Point2d pos(Xe[c], Ye[r]);
				val = QuadTree3D::bilinear_interpolation(pos, &ag);

				//rotateOnePnt(pos[0], pos[1], cx, cy, rotation, resx, resy);

				Vertex* nv = new Vertex(pos[0], pos[1], val, vid++);
				vdata.push_back(nv);
			}
		}
	}

	int mfgrid_esize = (mfgrid->nrow + 1) * (mfgrid->ncol + 1);
	int perrow_esize = mfgrid->ncol + 1;
	//new assign the vertex pointer to boxes
	for(int k = 0; k < mfgrid->nlay; k++)
	{
		for(int r = 0; r <	mfgrid->nrow; r++)
		{
			for(int c = 0; c < mfgrid->ncol; c++)
			{
				int bid=mfgrid->get_nodeid(k,r,c);
				Box * box=mfgrid->nodegroup[bid];
				boxVts[box->id].vertices[0] = vdata[ k * mfgrid_esize + r * perrow_esize + c];
				boxVts[box->id].vertices[1] = vdata[ k * mfgrid_esize + (r + 1) * perrow_esize + c];
				boxVts[box->id].vertices[2] = vdata[ k * mfgrid_esize + (r + 1) * perrow_esize + c + 1];
				boxVts[box->id].vertices[3] = vdata[ k * mfgrid_esize + r * perrow_esize + c + 1];

				boxVts[box->id].vertices[4] = vdata[ (k + 1) * mfgrid_esize + r * perrow_esize + c];
				boxVts[box->id].vertices[5] = vdata[ (k + 1) * mfgrid_esize + (r + 1) * perrow_esize + c];
				boxVts[box->id].vertices[6] = vdata[ (k + 1) * mfgrid_esize + (r + 1) * perrow_esize + c + 1];
				boxVts[box->id].vertices[7] = vdata[ (k + 1) * mfgrid_esize + r * perrow_esize + c + 1];
			}
		}
	}

}
int collectCell_share(ModflowGrid* mfgrid,  vector<Box*>& allBoxes, vector<int>* intvals, vector<float>* dbvals, bool without_inactive, bool one_based_numbering)
{
	int cellNum = 0;
	for(int k=0;k<mfgrid->nlay;k++)
	{
		for(int r=0;r<mfgrid->nrow;r++)
		{
			for(int c=0;c<mfgrid->ncol;c++)
			{
				int bid=mfgrid->get_nodeid(k,r,c);
				Box * b=mfgrid->nodegroup[bid];
				if(b->isLeaf)
				{
					if(without_inactive && (!b->active))
						continue;

					allBoxes.push_back(b);

					//store info
					if(one_based_numbering)
					{
						intvals[0].push_back(b->id + 1); intvals[1].push_back(b->number + 1); 
						intvals[2].push_back(k + 1); intvals[3].push_back(r + 1); intvals[4].push_back(c + 1);
					}
					else
					{
						intvals[0].push_back(b->id); intvals[1].push_back(b->number); 
						intvals[2].push_back(k); intvals[3].push_back(r); intvals[4].push_back(c);
					}

					double top=mfgrid->bot[k * mfgrid->nrow * mfgrid->ncol + r * mfgrid->ncol + c];
					double bot=mfgrid->bot[ (k+1) * mfgrid->nrow * mfgrid->ncol + r * mfgrid->ncol + c];

					dbvals[0].push_back((float)top); 
					dbvals[1].push_back((float)bot); 
					dbvals[2].push_back((float)(b->dx)); 
					dbvals[3].push_back((float)(b->dy));

					intvals[5].push_back(cellNum++);//the vtk_id
				}
			}
		}
	}
	return cellNum;
}

void writeModflowGrid2VTK_shareV(ModflowGrid* mfgrid, string vtkName, bool without_inactive, bool one_based_numbering)
{
	int boxNum = mfgrid->nodegroup.size();
	BoxVtx* boxVts = new BoxVtx[boxNum];
	vector<Vertex*> vdata;
	buildSharingVertices_mfgrid(mfgrid, vdata, boxVts);
		//begin to write
	double cx, cy, rotation;
	mfgrid->getRotatePara(cx, cy, rotation);

	/*****************************************************************************************/
	//collect all the boxes and necessary info
	vector<Box*> allBoxes;
	
	vector<string> intNames, dbNames;
	vector<int> intvals[6];// nodeid, nodeNums, layers, rows, cols, vtkid;
	intNames.push_back("nodeid");
	intNames.push_back("nodenumber"); 
	intNames.push_back("layer"); 
	intNames.push_back("row"); 
	intNames.push_back("col");
	intNames.push_back("vtkid");

	vector<float> dbvals[4];//tops, bottoms, delrs, delcs;
	dbNames.push_back("top"); 
	dbNames.push_back("bottom"); 
	dbNames.push_back("delr"); 
	dbNames.push_back("delc");

	int cellNum  = collectCell_share(mfgrid,  allBoxes, intvals, dbvals, without_inactive, one_based_numbering);

	//write the points, loop over each layer
	stringstream ss0;
	ss0<<vtkName.c_str();
	ss0<<".vtu";
	ofstream ofile(ss0.str().c_str());
	writeVTKHeader(ofile);
	
	/***************************************************************************************/
	ofile<<"<Piece NumberOfPoints=\""<<vdata.size()<<"\" NumberOfCells=\""<<allBoxes.size()<<"\" >";

	//write the points
	ofile<<"<Points>\n";
	ofile<<"<DataArray type=\"Float32\" NumberOfComponents=\"3\" format=\"ascii\">\n";
	int flash = 0;
	for(vector<Vertex*>::iterator vit = vdata.begin(); vit != vdata.end(); ++vit)
	{
		Vertex* nv = *vit;
		double resx, resy;
		rotateOnePnt(nv->pos[0], nv->pos[1], cx, cy, rotation,resx, resy);
		ofile<<(float)(resx)<<"\t"<<(float)(resy)<<"\t"<<(float)(nv->pos[2] * ZSCALE)<<"\t";

		if(flash % 5000 == 0)
		{
			ofile.flush();
		}
		flash++;
	}
	ofile<<"\n";
	ofile<<"</DataArray>\n";
	ofile<<"</Points>\n";

	//write the boxes
	writeCells(ofile,  allBoxes, intvals, dbvals, boxVts);

	//write the attributes
	ofile<<"<CellData Scalars=\"scalars\">\n";
	int intSize = intNames.size();
	for(int i = 0; i < intSize; i++)
	{
		writeDataArray(ofile, "Int32", intNames[i], intvals[i]);
	}

	int dbSize = dbNames.size();
	for(int i = 0; i < 4; i++)
	{
		writeDataArray(ofile, "Float32", dbNames[i], (dbvals[i]));
	}
	ofile<<"</CellData>\n";

	ofile<<"</Piece>\n";
	writeVTKTail(ofile);

	ofile.close();

	//clear memory
	for(vector<Vertex*>::iterator vit = vdata.begin(); vit != vdata.end(); ++vit)
	{
		delete (*vit);
	}
	delete[] boxVts;
}

/******************************************************************************************************************
* write VTK file with sharing vertices for quadtree
******************************************************************************************************************/
//void DFSCollectInfo(Box* b, int l, int r, int c, int& boxNum, 
void DFScollectCell_share(Box* b, int k, int r, int c, QuadTree3D* qtree, vector<Box*>& allBoxes, int& boxNum,  vector<int>* intvals, vector<float>* dbvals
	, bool one_based_numbering)
{
	if(b->isLeaf && b->active)
	{
		allBoxes.push_back(b);

		//store info
		if(one_based_numbering)
		{
			intvals[0].push_back(b->id + 1); intvals[1].push_back(b->number + 1); 
			intvals[2].push_back(k + 1); intvals[3].push_back(r + 1);intvals[4].push_back(c + 1);
		}
		else
		{
			intvals[0].push_back(b->id); intvals[1].push_back(b->number); 
			intvals[2].push_back(k); intvals[3].push_back(r);intvals[4].push_back(c);
		}

		double top=qtree->top[b->number];
		double bot=qtree->bot[b->number];

		dbvals[0].push_back((float)top); 
		dbvals[1].push_back((float)bot); 
		dbvals[2].push_back((float)(b->dx)); 
		dbvals[3].push_back((float)(b->dy));

		intvals[5].push_back(boxNum++);//the vtk_id
	}
	else if(b->isLeaf==false)
	{
		for(int i = 0; i < 4; i++)
		{
			DFScollectCell_share(b->pChildren[i], k, r, c, qtree, allBoxes, boxNum, intvals, dbvals, one_based_numbering);
		}
	}
}
int collectCell_share(QuadTree3D* qtree,  vector<Box*>& allBoxes, vector<int>* intvals, vector<float>* dbvals, bool one_based_numbering)
{
	ModflowGrid* mfgrid = qtree->getModflowGrid();
	int cellNum = 0;
	for(int k=0;k<mfgrid->nlay;k++)
	{
		for(int r=0;r<mfgrid->nrow;r++)
		{
			for(int c=0;c<mfgrid->ncol;c++)
			{
				int bid=mfgrid->get_nodeid(k,r,c);
				Box * box=qtree->nodegroup[bid];
				if(box->active)
					DFScollectCell_share(box, k, r, c, qtree, allBoxes, cellNum, intvals, dbvals, one_based_numbering);
			}
		}
	}
	return cellNum;
}

Vertex* getMid(Vertex* v1, Vertex* v2, vector<Vertex*>& vdata, int& vid)
{
	Vertex* nv = new Vertex( (v1->pos[0] + v2->pos[0])/2, (v1->pos[1] +v2->pos[1])/2, (v1->pos[2] + v2->pos[2])/2, vid++);
	vdata.push_back(nv);

	return nv;
}

void rotateOnePnt(double x, double y, double cx, double cy, double rotation,double& resx, double& resy)
{
	resx= cx +  (x - cx) * cos(rotation) - (y - cy)* sin(rotation);
	resy = cy + (y - cy)*cos(rotation) + (x - cx)*sin(rotation);
}


//////////////////////////////////////////////////////////////////////////////////////////////////////////
//fit a plane using 3 points
double fitZ(double** vv, double x, double y)
{
	double a = (vv[1][1] - vv[0][1]) * (vv[2][2] - vv[0][2]) - (vv[2][1] - vv[0][1]) * (vv[1][2] - vv[0][2]);
	double b = (vv[1][2] - vv[0][2]) * (vv[2][0] - vv[0][0]) - (vv[2][2] - vv[0][2]) * (vv[1][0] - vv[0][0]);
	double c = (vv[1][0] - vv[0][0]) * (vv[2][1] - vv[0][1]) - (vv[2][0] - vv[0][0]) * (vv[1][1] - vv[0][1]);
	double d = -(a * vv[0][0] + b * vv[0][1] + c * vv[0][2]);
	
	assert(fabs(c)>=1e-10);
	double z = - (a * x + b * y + d) / c;
	
	return z;
}

int getDepth(int n[4])
{
	if(n[0] == n[1] && n[1] == n[2] && n[2] == n[3])
		return -1;
	for(int i = 0; i < 4; i++)
	{
		if(n[i] > n[(i+1)%4])
			return i;
	}
	assert(false);
	return -1;
}
bool no_collinear_add(set<Box*>& bset, Box* b)
{
	if(bset.size() == 2)
	{
		set<Box*>::iterator bit = bset.begin();

		Box* b1 = *(bset.begin());
		Box* b2 = *(++bit);
		if(bset.find(b) == bset.end())
		{
			double co_eff = fabs((b2->x - b1->x) * (b->y - b1->y) - (b->x - b1->x) * (b2->y - b1->y));
			if( co_eff < 1e-6)
				return false;
		}
	}
	bset.insert(b);
	return true;
}

//-------------------------------------------------------------------------------------------------

//build the mfgrid ascii grid
void build_mf_ascii(QuadTree3D* qtree, vector<mf_ascii_grid*>& top_ags, vector<mf_ascii_grid*>& bot_ags)
{
	ModflowGrid* mfgrid = qtree->getModflowGrid();

	int csize = mfgrid->nrow * mfgrid->ncol;
	for(int k = 0; k < mfgrid->nlay; k++)
	{
		double* top_zgrid = &(mfgrid->top[csize * k]);
		double* bot_zgrid = &(mfgrid->bot[csize * k]);
		mf_ascii_grid* top_ag = new mf_ascii_grid(mfgrid, top_zgrid, mfgrid->ncol, mfgrid->nrow); 
		mf_ascii_grid* bot_ag = new mf_ascii_grid(mfgrid, bot_zgrid, mfgrid->ncol, mfgrid->nrow);
		top_ags.push_back(top_ag);
		bot_ags.push_back(bot_ag);
	}
}



//
// for QuadTree3D sharing Vertices
//


struct VTK_V
{
	VTK_V()
	{
		num = -1;
		for(int i = 0; i < 8; i++)
			box[i] = NULL;
	}

	Point3d pos;
	Box* box[8];
	int num;
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
struct VTK_E
{
	VTK_E()
	{
		child_v = NULL; 
		v[0] = NULL; v[1] = NULL;
		child_e[0] = NULL; child_e[1] = NULL;
	}
	VTK_E(VTK_V* v1, VTK_V* v2)
	{
		child_v = NULL; 
		v[0] = v1; v[1] = v2;
		child_e[0] = NULL; child_e[1] = NULL;
	}

	VTK_V* child_v;
	VTK_E* child_e[2];
	VTK_V* v[2];

	//
	map<VTK_E*, VTK_E*> sub_e;
};

struct VTK_Box
{
	VTK_Box()
	{
		for(int i = 0; i < 8; i++)
			e[i] = NULL;
	}
	VTK_E* e[8];
};

//
//	build 
//
void buildVE(QuadTree3D* qtree, vector<VTK_V*>& vdata, vector<VTK_E*>& edata, VTK_Box* boxdata)
{
	ModflowGrid* mfgrid = qtree->getModflowGrid();
	double* Xe = mfgrid->get_local_Xe_array();
	double* Ye = mfgrid->get_local_Ye_array();

	//////////////////////////////////////////////////////////////////////////////
	//initialize sharing v
	for(int k = 0; k <= mfgrid->nlay; k++)
	{
		for(int r = 0; r <= mfgrid->nrow; r++)
		{
			for(int c = 0; c <= mfgrid->ncol; c++)
			{
				VTK_V* nv = new VTK_V();
				nv->pos[0] = Xe[c]; nv->pos[1] = Ye[r];
				vdata.push_back(nv);
			}
		}
	}
	//////////////////////////////////////////////////////////////////////////////
	//initializing the sharing edge
	int mfgrid_vsize = (mfgrid->nrow + 1) * (mfgrid->ncol + 1);
	int perrow_vsize = mfgrid->ncol + 1;
	for(int k = 0; k <= mfgrid->nlay; k++)
	{
		//initializing the row edges
		for(int r = 0; r <= mfgrid->nrow; r++)
		{
			for(int c = 0; c < mfgrid->ncol; c++)
			{
				VTK_V* v1 = vdata[ k * mfgrid_vsize + r * perrow_vsize + c];
				VTK_V* v2 = vdata[ k * mfgrid_vsize + r * perrow_vsize + c + 1];
				VTK_E* e = new VTK_E(v1, v2);
				edata.push_back(e);
			}
		}

		//initializing the col edges
		for(int c = 0; c <= mfgrid->ncol; c++)
		{
			for(int r = 0; r < mfgrid->nrow; r++)
			{
				VTK_V* v1 = vdata[ k * mfgrid_vsize + r * perrow_vsize + c];
				VTK_V* v2 = vdata[ k * mfgrid_vsize + (r + 1) * perrow_vsize + c];
				VTK_E* e = new VTK_E(v1, v2);
				edata.push_back(e);
			}
		}
	}
}
void assignBoxEdge(Box* box, vector<VTK_V*>& vdata, vector<VTK_E*>& edata, VTK_Box* boxdata)
{
	int id[4];
	VTK_Box& vtkBox = boxdata[box->id];
	for(int i = 0; i < 4; i++)
		id[i] = box->pChildren[i]->id;

	///////////////////////////////////////////////////////

	boxdata[id[0]].e[0] = vtkBox.e[0]->child_e[0];
	boxdata[id[0]].e[1] = vtkBox.e[0]->sub_e[vtkBox.e[2]];
	boxdata[id[0]].e[2] = vtkBox.e[3]->sub_e[vtkBox.e[1]];
	boxdata[id[0]].e[3] = vtkBox.e[3]->child_e[0];
	//----------------
	boxdata[id[0]].e[4] = vtkBox.e[4]->child_e[0];
	boxdata[id[0]].e[5] = vtkBox.e[4]->sub_e[vtkBox.e[6]];
	boxdata[id[0]].e[6] = vtkBox.e[7]->sub_e[vtkBox.e[5]];
	boxdata[id[0]].e[7] = vtkBox.e[7]->child_e[0];


	///////////////////////////////////////////////////////


	boxdata[id[1]].e[0] = vtkBox.e[3]->sub_e[vtkBox.e[1]];
	boxdata[id[1]].e[1] = vtkBox.e[2]->sub_e[vtkBox.e[0]];
	boxdata[id[1]].e[2] = vtkBox.e[2]->child_e[0];
	boxdata[id[1]].e[3] = vtkBox.e[3]->child_e[1];
	//----
	boxdata[id[1]].e[4] = vtkBox.e[7]->sub_e[vtkBox.e[5]];
	boxdata[id[1]].e[5] = vtkBox.e[6]->sub_e[vtkBox.e[4]];
	boxdata[id[1]].e[6] = vtkBox.e[6]->child_e[0];
	boxdata[id[1]].e[7] = vtkBox.e[7]->child_e[1];

	///////////////////////////////////////////////////////

	boxdata[id[2]].e[0] = vtkBox.e[1]->sub_e[vtkBox.e[3]];
	boxdata[id[2]].e[1] = vtkBox.e[1]->child_e[1];
	boxdata[id[2]].e[2] = vtkBox.e[2]->child_e[1];
	boxdata[id[2]].e[3] = vtkBox.e[2]->sub_e[vtkBox.e[0]];
	//------
	boxdata[id[2]].e[4] = vtkBox.e[5]->sub_e[vtkBox.e[7]];
	boxdata[id[2]].e[5] = vtkBox.e[5]->child_e[1];
	boxdata[id[2]].e[6] = vtkBox.e[6]->child_e[1];
	boxdata[id[2]].e[7] = vtkBox.e[6]->sub_e[vtkBox.e[4]];

	
	///////////////////////////////////////////////////////

	boxdata[id[3]].e[0] = vtkBox.e[0]->child_e[1];
	boxdata[id[3]].e[1] = vtkBox.e[1]->child_e[0];
	boxdata[id[3]].e[2] = vtkBox.e[1]->sub_e[vtkBox.e[3]];
	boxdata[id[3]].e[3] = vtkBox.e[0]->sub_e[vtkBox.e[2]];
	//-------------
	boxdata[id[3]].e[4] = vtkBox.e[4]->child_e[1];
	boxdata[id[3]].e[5] = vtkBox.e[5]->child_e[0];
	boxdata[id[3]].e[6] = vtkBox.e[5]->sub_e[vtkBox.e[7]];
	boxdata[id[3]].e[7] = vtkBox.e[4]->sub_e[vtkBox.e[6]];

}
void DFSbuildBox(Box* box, vector<VTK_V*>& vdata, vector<VTK_E*>& edata, VTK_Box* boxdata)
{
	if(box->isLeaf) return;

	VTK_Box& vtkBox = boxdata[box->id];
	for(int k = 0; k < 2; k++)
	{
		for(int t = 0; t < 4; t++)
		{
			VTK_E* te = boxdata[box->id].e[k * 4 + t];
			if(te->child_v == NULL)
			{
				VTK_V* nv = new VTK_V();
				nv->pos[0] = (te->v[0]->pos[0] + te->v[1]->pos[0])/2;
				nv->pos[1] = (te->v[0]->pos[1] + te->v[1]->pos[1])/2;
				vdata.push_back(nv);
				te->child_v = nv;

				VTK_E* e1 = new VTK_E(te->v[0], nv);
				VTK_E* e2 = new VTK_E(nv, te->v[1]);
				te->child_e[0] = e1; te->child_e[1] = e2;
				edata.push_back(e1); edata.push_back(e2);
			}
		}
	}

	// determine all the edges
	if(vtkBox.e[0]->sub_e.find(boxdata[box->id].e[2]) ==vtkBox.e[0]->sub_e.end())
	{
		//center v
		VTK_V* cv1 = new VTK_V();
		cv1->pos[0] = (vtkBox.e[0]->v[0]->pos[0] + vtkBox.e[2]->v[1]->pos[0]) / 2;
		cv1->pos[1] = (vtkBox.e[0]->v[0]->pos[1] + vtkBox.e[2]->v[1]->pos[1]) / 2;

		VTK_E* e02 = new VTK_E(vtkBox.e[0]->child_v, cv1);
		VTK_E* e20 = new VTK_E(cv1, vtkBox.e[2]->child_v);
		VTK_E* e31 = new VTK_E(vtkBox.e[3]->child_v, cv1);
		VTK_E* e13 = new VTK_E(cv1, vtkBox.e[1]->child_v);

		edata.push_back(e02); edata.push_back(e20); edata.push_back(e31); edata.push_back(e13);

		vtkBox.e[0]->sub_e.insert(make_pair(vtkBox.e[2], e02));
		vtkBox.e[2]->sub_e.insert(make_pair(vtkBox.e[0], e20));
		vtkBox.e[3]->sub_e.insert(make_pair(vtkBox.e[1], e31));
		vtkBox.e[1]->sub_e.insert(make_pair(vtkBox.e[3], e13));
	}

	if(vtkBox.e[4]->sub_e.find(vtkBox.e[6]) == vtkBox.e[4]->sub_e.end())
	{
		//center v
		VTK_V* cv2 = new VTK_V();
		cv2->pos[0] = (vtkBox.e[4]->v[0]->pos[0] + vtkBox.e[6]->v[1]->pos[0]) / 2;
		cv2->pos[1] = (vtkBox.e[4]->v[0]->pos[1] + vtkBox.e[6]->v[1]->pos[1]) / 2;

		VTK_E* e02 = new VTK_E(vtkBox.e[4]->child_v, cv2);
		VTK_E* e20 = new VTK_E(cv2, vtkBox.e[6]->child_v);
		VTK_E* e31 = new VTK_E(vtkBox.e[7]->child_v, cv2);
		VTK_E* e13 = new VTK_E(cv2, vtkBox.e[5]->child_v);
		
		edata.push_back(e02); edata.push_back(e20); edata.push_back(e31); edata.push_back(e13);

		vtkBox.e[4]->sub_e.insert(make_pair(vtkBox.e[6], e02));
		vtkBox.e[6]->sub_e.insert(make_pair(vtkBox.e[4], e20));
		vtkBox.e[7]->sub_e.insert(make_pair(vtkBox.e[5], e31));
		vtkBox.e[5]->sub_e.insert(make_pair(vtkBox.e[7], e13));
	}

	//assign edges
	assignBoxEdge(box, vdata, edata, boxdata);

	//boxdata
	for(int i = 0; i < 4; i++)
	{
		DFSbuildBox(box->pChildren[i], vdata, edata, boxdata);
	}
}
void buildBox(QuadTree3D* qtree, vector<VTK_V*>& vdata, vector<VTK_E*>& edata, VTK_Box* boxdata)
{
	ModflowGrid* mfgrid = qtree->getModflowGrid();
	/////////////////////////////////////////////////////////////////////////////////////////
	//initialize the root box
	int mfgrid_esize = mfgrid->ncol * ( mfgrid->nrow + 1) + mfgrid->nrow * (mfgrid->ncol + 1);
	for(int k = 0; k < mfgrid->nlay; k++)
	{
		for(int r = 0; r < mfgrid->nrow; r++)
		{
			for(int c = 0; c < mfgrid->ncol; c++)
			{
				int bid = mfgrid->get_nodeid(k, r, c);
				Box* box = qtree->nodegroup[bid];
				
				for(int t = 0; t < 2; t++)
				{
					int baseid = (k + t) * mfgrid_esize;
					int id1 = baseid + (mfgrid->nrow + 1) * mfgrid->ncol + c * mfgrid->nrow + r;
					int id2 = baseid + (r + 1) * mfgrid->ncol + c;
					int id3 = baseid + (mfgrid->nrow + 1) * mfgrid->ncol + (c + 1) * mfgrid->nrow + r;
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
	for(int k = 0; k < mfgrid->nlay; k++)
	{
		for(int r = 0; r < mfgrid->nrow; r++)
		{
			for(int c = 0; c < mfgrid->ncol; c++)
			{
				int bid = mfgrid->get_nodeid(k, r, c);
				Box* box = qtree->nodegroup[bid];

				DFSbuildBox(box, vdata, edata, boxdata);
			}
		}
	}

}
void buildSharing(QuadTree3D* qtree, vector<VTK_V*>& vdata, vector<VTK_E*>& edata, VTK_Box* boxdata)
{
	int bsize = qtree->nodegroup.size();

	// build VE
	buildVE(qtree, vdata, edata, boxdata);

	//build Box
	buildBox(qtree, vdata, edata, boxdata);
}

//
// analysis
//

void determineVtoBox(QuadTree3D* qtree, vector<VTK_V*>& vdata, vector<VTK_E*>& edata, VTK_Box* boxdata)
{
	int boxNum = qtree->nodegroup.size();
	for(int i = 0; i < boxNum; i++)
	{
		VTK_Box& vb = boxdata[i];
		Box* tbox = qtree->nodegroup[i];
		if(!(tbox->isLeaf))
			continue;

		vb.e[0]->v[0]->box[6] = tbox;
		vb.e[0]->v[1]->box[7] = tbox;
		vb.e[2]->v[0]->box[5] = tbox;
		vb.e[2]->v[1]->box[4] = tbox;
		
		vb.e[4]->v[0]->box[2] = tbox;
		vb.e[4]->v[1]->box[3] = tbox;
		vb.e[6]->v[0]->box[1] = tbox;
		vb.e[6]->v[1]->box[0] = tbox;

		//child 
		if(vb.e[0]->child_v != NULL)
		{
			vb.e[0]->child_v->box[6] = tbox;
			vb.e[0]->child_v->box[7] = tbox;
		}
		if(vb.e[1]->child_v != NULL)
		{
			vb.e[1]->child_v->box[4] = tbox;
			vb.e[1]->child_v->box[7] = tbox;
		}
		if(vb.e[2]->child_v != NULL)
		{
			vb.e[2]->child_v->box[4] = tbox;
			vb.e[2]->child_v->box[5] = tbox;
		}
		if(vb.e[3]->child_v != NULL)
		{
			vb.e[3]->child_v->box[5] = tbox;
			vb.e[3]->child_v->box[6] = tbox;
		}

		//child
		if(vb.e[4]->child_v != NULL)
		{
			vb.e[4]->child_v->box[2] = tbox;
			vb.e[4]->child_v->box[3] = tbox;
		}
		if(vb.e[5]->child_v != NULL)
		{
			vb.e[5]->child_v->box[0] = tbox;
			vb.e[5]->child_v->box[3] = tbox;
		}
		if(vb.e[6]->child_v != NULL)
		{
			vb.e[6]->child_v->box[1] = tbox;
			vb.e[6]->child_v->box[2] = tbox;
		}
		if(vb.e[7]->child_v != NULL)
		{
			vb.e[7]->child_v->box[1] = tbox;
			vb.e[7]->child_v->box[2] = tbox;
		}
	}
}
void number_share_vtk_v(QuadTree3D* qtree, vector<VTK_V*>& vdata, VTK_Box* boxdata, vector<VTK_V*>& vtk_vertices)
{
	int num = 0;
	int bsize = qtree->nodegroup.size();

	for(int i = 0; i < bsize; i++)
	{
		Box* b = qtree->nodegroup[i];
		if(!(b->isLeaf && b->active))
			continue;
		VTK_Box& vb = boxdata[b->id];
		for(int k = 0; k < 2; k++)
		{
			if(vb.e[0 + 4 * k]->v[0]->num == -1)
			{
				vb.e[0 + 4 * k]->v[0]->num = num++;
				vtk_vertices.push_back(vb.e[0 + 4 * k]->v[0]);
			}
			if(vb.e[0 + 4 * k]->v[1]->num == -1)
			{
				vb.e[0 + 4 * k]->v[1]->num = num++;
				vtk_vertices.push_back(vb.e[0 + 4 * k]->v[1]);
			}
			if(vb.e[2 + 4 * k]->v[0]->num == -1)
			{
				vb.e[2 + 4 * k]->v[0]->num = num++;
				vtk_vertices.push_back(vb.e[2 + 4 * k]->v[0]);
			}
			if(vb.e[2 + 4 * k]->v[1]->num == -1)
			{
				vb.e[2 + 4 * k]->v[1]->num = num++;
				vtk_vertices.push_back(vb.e[2 + 4 * k]->v[1]);
			}
		}
	}
}

//assume at most 1 depth difference
double smooth_interpolate(vector<double>& dvals, vector<int>& depth, vector<Box*>& bs, double px, double py)
{
	set<int> dset;
	int dsize = depth.size();
	for(int i = 0; i < dsize; i++)
	{
		dset.insert(depth[i]);
	}
	//average
	if(dset.size() == 1 || dsize == 1)
	{
		double dsum(0);
		for(int i = 0; i < dsize; i++)
			dsum += dvals[i];
		return dsum/dsize;
	}

	//interpolation
	if(dsize == 4)
	{
		double d[2] = {0};
		for(int k = 0; k < 2; k++)
		{
			if(depth[k] == depth[k + 2])
				d[k] = (dvals[k] + dvals[k + 2])/2;
			else if(depth[k] < depth[k + 2])
				d[k] = (dvals[k] + dvals[k + 2]*2)/3;
			else if(depth[k] > depth[k + 2])
				d[k] = (dvals[k]*2 + dvals[k+2])/3;
			else 
				assert(false);
		}
		return (d[0] + d[1])/2;
	}
	
	//do a plane fitting
	if(dsize == 3)
	{
		double** vv = new double*[3];
		for(int i = 0; i < 3; i++)
		{
			vv[i] = new double[3];
			memset(vv[i], 0, sizeof(double) * 3);
		}

		for(int i = 0; i < 3; i++)
		{
			vv[i][0] = bs[i]->x; vv[i][1] = bs[i]->y; vv[i][2] = dvals[i];
		}

		double val = fitZ(vv, px, py);

		for(int i = 0; i < 3; i++)
			delete[] vv[i];
		delete[] vv;

		return val;
	}

	if(dsize == 2)
	{
		if(depth[0] == depth[1])
			return (dvals[0] + dvals[1])/2;
		if(depth[0] < depth[1])
			return (dvals[0] + 2*dvals[1])/3;
		if(depth[0] > depth[1])
			return (dvals[0]*2 + dvals[1])/3;
	}
	assert(false);
	return 0;
}

void interpolate_z(QuadTree3D* qtree, vector<VTK_V*>& vdata, vector<VTK_E*>& edata, VTK_Box* boxdata, vector<VTK_V*>& vtk_vertices)
{
	double* ztop = qtree->top;
	double* zbot = qtree->bot;
	for(vector<VTK_V*>::iterator vit = vtk_vertices.begin(); vit != vtk_vertices.end(); ++vit)
	{
		VTK_V* v = *vit;
		set<Box*> bset1, bset2;
		vector<double> dvals1, dvals2;
		vector<int> depth1, depth2;
		//vector<double> zs;
		for(int i = 0; i < 4; i++)
		{
			if(v->box[i] != NULL && v->box[i]->isLeaf && v->box[i]->active)
			{
				if(bset1.find(v->box[i]) == bset1.end())
				{
					bset1.insert(v->box[i]); dvals1.push_back(zbot[v->box[i]->number]);
					depth1.push_back(v->box[i]->depth);
				}
				//zs.push_back(zbot[v->box[i]->number]);
			}
			if(v->box[i + 4] != NULL && v->box[i + 4]->isLeaf && v->box[i + 4]->active)
			{
				if(bset2.find(v->box[i + 4]) == bset2.end())
				{
					bset2.insert(v->box[i + 4]); dvals2.push_back(ztop[v->box[i + 4]->number]);
					depth2.push_back(v->box[i + 4]->depth);
				}
				//zs.push_back(ztop[v->box[i + 4]->number]);
			}
		}

		double dt[2] = {0};
		//1
		if(bset1.size() > 0)
		{
			vector<Box*> vb1(bset1.begin(), bset1.end());
			dt[0] = smooth_interpolate(dvals1, depth1, vb1, v->pos[0], v->pos[1]);
		}
		//2
		if(bset2.size() > 0)
		{
			vector<Box*> vb2(bset2.begin(), bset2.end());
			dt[1] = smooth_interpolate(dvals2, depth2, vb2, v->pos[0], v->pos[1]);
		}

		//double zsum(0);
		//for(vector<double>::iterator dit = zs.begin(); dit != zs.end(); ++dit)
		//{
		//	zsum += (*dit);
		//}
		if(bset1.size() > 0 && bset2.size() > 0)
			v->pos[2] = (dt[0] + dt[1])/2;
		else if(bset1.size() > 0)
			v->pos[2] = dt[0];
		else 
			v->pos[2] = dt[1];
	}
}

//
// write
//
void writeCells(ofstream& ofile,  vector<Box*>& allBoxes, vector<int>* intvals, vector<float>* dbvals, VTK_Box* boxdata)
{
	//ModflowGrid* mfgrid = qtree->getModflowGrid();
	/***************************************************************************************/
	//write the cells
	ofile<<"<Cells>\n";
	//write the connectivity
	ofile<<"<DataArray type=\"Int32\" Name=\"connectivity\" format=\"ascii\">\n";

	int rid = 0;
	VTK_V* vids[8];
	int cellNum = allBoxes.size();
	for(vector<Box*>::iterator bit = allBoxes.begin(); bit != allBoxes.end(); ++bit)
	{
		Box* b = *bit;
		VTK_Box& s_box = boxdata[b->id];
		vids[0] = s_box.e[0]->v[0]; vids[1] = s_box.e[0]->v[1]; vids[2] = s_box.e[1]->v[1]; vids[3] = s_box.e[2]->v[0];
		vids[4] = s_box.e[4]->v[0]; vids[5] = s_box.e[4]->v[1]; vids[6] = s_box.e[5]->v[1]; vids[7] = s_box.e[6]->v[0];

		ofile<<vids[0]->num<<"\t"<<vids[1]->num<<"\t"<<vids[3]->num<<"\t"<<vids[2]->num<<"\t"
			<<vids[4]->num<<"\t"<<vids[5]->num<<"\t"<<vids[7]->num<<"\t"<<vids[6]->num<<"\n";
		//ofile<<vids[0]->num<<"\t"<<vids[6]->num<<"\t"<<vids[2]->num<<"\t"<<vids[4]->num<<"\t"
		//	<<vids[1]->num<<"\t"<<vids[7]->num<<"\t"<<vids[3]->num<<"\t"<<vids[5]->num<<"\n";

		rid++;
		if(rid % 5000 == 0)
			ofile.flush();
	}
	ofile<<"</DataArray>\n";

	//write the offset
	int sid, id = 0;
	sid = 8;
	ofile<<"<DataArray type=\"Int32\" Name=\"offsets\" format=\"ascii\">\n";
	for(int i =0; i < cellNum; i++)
	{
		ofile<<sid<<"\t";
		sid+=8;

		id++;
		if(id%5000==0) ofile.flush();
	}
	ofile<<"\n";
	ofile<<"</DataArray>\n";
	//write the types
	id = 0;
	ofile<<"<DataArray type=\"UInt8\" Name=\"types\" format=\"ascii\">\n";
	for(int i = 0; i < cellNum; i++)
	{
		ofile<<"11\t";
		id++;
		if(id%5000==0)
			ofile.flush();
	}
	ofile<<"</DataArray>\n";
	ofile<<"</Cells>\n";
}
void write_vtu_file(QuadTree3D* qtree, vector<VTK_V*>& vdata, vector<VTK_E*>& edata, VTK_Box* boxdata, vector<VTK_V*>& vtk_vertices, string vtkName)
{
	ModflowGrid* mfgrid = qtree->getModflowGrid();
	//begin to write
	double cx, cy, rotation;
	mfgrid->getRotatePara(cx, cy, rotation);

	vector<Box*> allBoxes;
	
	vector<string> intNames, dbNames;
	vector<int> intvals[6];// nodeid, nodeNums, layers, rows, cols, vtkid;
	intNames.push_back("nodeid");
	intNames.push_back("nodenumber"); 
	intNames.push_back("layer"); 
	intNames.push_back("row"); 
	intNames.push_back("col");
	intNames.push_back("vtkid");

	vector<float> dbvals[4];//tops, bottoms, delrs, delcs;
	dbNames.push_back("top"); 
	dbNames.push_back("bottom"); 
	dbNames.push_back("delr"); 
	dbNames.push_back("delc");
	//
	int cellNum = collectCell_share(qtree, allBoxes, intvals, dbvals, qtree->get_one_based_node_numbering());

		//write the points, loop over each layer
	stringstream ss0;
	ss0<<vtkName.c_str();
	ss0<<".vtu";
	ofstream ofile(ss0.str().c_str());
	writeVTKHeader(ofile);
	
	/***************************************************************************************/
	ofile<<"<Piece NumberOfPoints=\""<<vtk_vertices.size()<<"\" NumberOfCells=\""<<allBoxes.size()<<"\" >";

	//write the points
	ofile<<"<Points>\n";
	ofile<<"<DataArray type=\"Float32\" NumberOfComponents=\"3\" format=\"ascii\">\n";
	int flash = 0;
	for(vector<VTK_V*>::iterator vit = vtk_vertices.begin(); vit != vtk_vertices.end(); ++vit)
	{
		VTK_V* sv = *vit;
		double resx, resy;
		rotateOnePnt(sv->pos[0], sv->pos[1], cx, cy, rotation,resx, resy);

		ofile<<(float)(resx)<<"\t"<<(float)(resy)<<"\t"<<(float)(sv->pos[2] * ZSCALE)<<"\t";

		if(flash % 5000 == 0)
		{
			ofile.flush();
		}
		flash++;
	}

	ofile<<"\n";
	ofile<<"</DataArray>\n";
	ofile<<"</Points>\n";

	//write the boxes
	writeCells(ofile,  allBoxes, intvals, dbvals, boxdata);

	//write the attributes
	ofile<<"<CellData Scalars=\"scalars\">\n";
	int intSize = intNames.size();
	for(int i = 0; i < intSize; i++)
	{
		writeDataArray(ofile, "Int32", intNames[i], intvals[i]);
	}

	int dbSize = dbNames.size();
	for(int i = 0; i < 4; i++)
	{
		writeDataArray(ofile, "Float32", dbNames[i], (dbvals[i]));
	}
	ofile<<"</CellData>\n";

	ofile<<"</Piece>\n";
	writeVTKTail(ofile);

	ofile.close();
}

//
void writeQuadtreeGrid2VTK_shareV(QuadTree3D* qtree, string vtkName)
{
	qtree->number_nodes();
	
	int boxNum = qtree->nodegroup.size();

	vector<VTK_V*> vdata;
	vector<VTK_E*> edata;
	VTK_Box* boxdata = new VTK_Box[boxNum];

	//build the sharing v
	buildSharing(qtree, vdata, edata, boxdata);
	//determine the neighbor boxes for v
	determineVtoBox(qtree, vdata, edata, boxdata);

	//number the vertices
	vector<VTK_V*> vtk_vertices;
	number_share_vtk_v(qtree, vdata,	boxdata, vtk_vertices);

	//export to VTK file
	interpolate_z(qtree, vdata, edata, boxdata, vtk_vertices);
	
	write_vtu_file(qtree, vdata, edata, boxdata, vtk_vertices, vtkName);
	
	//clear memory
	for(vector<VTK_V*>::iterator vit = vdata.begin(); vit != vdata.end(); ++vit)
	{
		delete (*vit);
	}
	for(vector<VTK_E*>::iterator eit = edata.begin(); eit != edata.end(); ++eit)
	{
		delete (*eit);
	}
	delete[] boxdata; 

}
