#include "intersecReadWrite.h"
#include <sstream>
using namespace std;
using namespace cusg;

const string _IFO_DEL=",";

bool checkAttributesContain(const char* attriName, vector<string>& attributes)
{
	for (vector<string>::iterator sit = attributes.begin(); sit != attributes.end(); ++sit)
	{
		string & str = *sit;
		if(str == attriName || str =="all")
			return true;
		else if(str == "none")
			return false;
	}
	return false;
}
void checkShpAttrContain(GridIntersection* gridInter, 	vector<bool>& bIntFld, vector<bool>& bDBFld,vector<bool>& bStrFld, vector<bool>& bLogFld,
	vector<string>& iName,vector<string>&dName, vector<string>&sName, vector<string>&lName, vector<string>& attributes)
{
	vector<ShapeIntField>& m_intFlds = gridInter->m_shp_reader.getIntFields();
	vector<ShapeDoubleField>& m_dbFlds = gridInter->m_shp_reader.getDBFields();
	vector<ShapeStringField>& m_strFlds = gridInter->m_shp_reader.getStrFields();
	vector<ShapeLogField>& m_logFlds = gridInter->m_shp_reader.getLogicalFields();
	for(vector<ShapeIntField>::iterator sit = m_intFlds.begin(); sit !=m_intFlds.end(); ++sit)
	{
		if(checkAttributesContain(sit->name.c_str(), attributes))
		{
			bIntFld.push_back(true); iName.push_back(sit->name);
		}
		else
			bIntFld.push_back(false);
	}
	for(vector<ShapeDoubleField>::iterator sit = m_dbFlds.begin(); sit !=m_dbFlds.end(); ++sit)
	{
		if(checkAttributesContain(sit->name.c_str(), attributes))
		{
			bDBFld.push_back(true);dName.push_back(sit->name);
		}
		else
			bDBFld.push_back(false);
	}
	for(vector<ShapeStringField>::iterator sit = m_strFlds.begin(); sit !=m_strFlds.end(); ++sit)
	{
		if(checkAttributesContain(sit->name.c_str(), attributes))
		{
			bStrFld.push_back(true); sName.push_back(sit->name);
		}
		else
			bStrFld.push_back(false);
	}
	for(vector<ShapeLogField>::iterator sit = m_logFlds.begin(); sit !=m_logFlds.end(); ++sit)
	{
		if(checkAttributesContain(sit->name.c_str(), attributes))
		{
			bLogFld.push_back(true); lName.push_back(sit->name);
		}
		else
			bLogFld.push_back(false);
	}
}
void addShpAttrVal( int iShape, GridIntersection* gridInter, vector<bool>& bIntFld, vector<bool>& bDBFld,vector<bool>& bStrFld, vector<bool>& bLogFld,
	vector<int>& ivals, vector<double>& dvals, vector<string>& svals, vector<string>& lvals)
{
	vector<ShapeIntField>& m_intFlds = gridInter->m_shp_reader.getIntFields();
	vector<ShapeDoubleField>& m_dbFlds = gridInter->m_shp_reader.getDBFields();
	vector<ShapeStringField>& m_strFlds = gridInter->m_shp_reader.getStrFields();
	vector<ShapeLogField>& m_logFlds = gridInter->m_shp_reader.getLogicalFields();
	vector<bool>::iterator ibit = bIntFld.begin();
	for(vector<ShapeIntField>::iterator sit = m_intFlds.begin(); sit !=m_intFlds.end(); ++sit,++ibit)
	{
		if(*ibit)
			ivals.push_back(sit->vals[ iShape ]);
	}
	vector<bool>::iterator idit = bDBFld.begin();
	for(vector<ShapeDoubleField>::iterator sit = m_dbFlds.begin(); sit !=m_dbFlds.end(); ++sit, ++idit)
	{
		if(*idit)
			dvals.push_back(sit->vals[ iShape ]);
	}
	vector<bool>::iterator sbit = bStrFld.begin();
	for(vector<ShapeStringField>::iterator sit = m_strFlds.begin(); sit !=m_strFlds.end(); ++sit, ++sbit)
	{
		if(*sbit)
			svals.push_back(sit->vals[ iShape]);
	}
	vector<bool>::iterator lbit = bLogFld.begin();
	for(vector<ShapeLogField>::iterator sit = m_logFlds.begin(); sit !=m_logFlds.end(); ++sit, ++lbit)
	{
		if(*lbit)
			lvals.push_back(sit->vals[iShape]);
	}
}

//void  addShpAttrVal(GridIntersection* gridInter, 	vector<bool>& bIntFld, vector<bool>& bDBFld,vector<bool>& bStrFld, vector<bool>& bLogFld,
//	vector<int>& ivals, vector<double>& dvals, vector<string>& svals, vector<string>& lvals)
//{
//
//}


void addDBHandleFields(DBFHandle& dbfhdl, int intFldNum, vector<string>& intFldNames, 
	int doubleFldNum,vector<string>& dbFldNames,vector<int>& intFldIDs,vector<int>& dbFldIDs)
{
	for(int i = 0; i < intFldNum; i++)
	{
		int tid = DBFAddField(dbfhdl, intFldNames[i].c_str(), FTInteger,20,0);
		intFldIDs.push_back(tid);
	}
	for(int i = 0; i < doubleFldNum; i++)
	{
		int tid = DBFAddField(dbfhdl, dbFldNames[i].c_str(), FTDouble,20,8);
		dbFldIDs.push_back(tid);
	}
}
void addDBHandleStrFields(DBFHandle& dbfhdl, int strFldNum, vector<string>& strFldNames, vector<int>& strFldIDs)
{
	for(int i = 0; i < strFldNum; i++)
	{
		int tid = DBFAddField(dbfhdl, strFldNames[i].c_str(), FTString,50,0);
		strFldIDs.push_back(tid);
	}
}
//write int and double attributes
void writeIntDBVal(DBFHandle& dbfhdl, int entity, vector<int>& intFldIDs, vector<int>& intVals, vector<int>& dbFldIDs, vector<double>& dbVals,
	vector<int>& strFldIDs, vector<string>& strVals,vector<int>& logFldIDs, vector<string>& logVals, ofstream& ofile, bool isTxt)
{
	int intFldNum = intVals.size();
	int dbFldNum = dbVals.size();
	int strFldNum = strVals.size();
	int logFldNum = logVals.size();
	//write integer fields
	int success = 0;
	for(int i = 0; i < intFldNum; ++i)
	{
		if(!isTxt)
			success = DBFWriteIntegerAttribute(dbfhdl, entity, intFldIDs[i], intVals[i]);
		else
			ofile<<intVals[i]<<_IFO_DEL;
	}
	for (int i = 0; i < dbFldNum; ++i)
	{
		if(!isTxt)
			success = DBFWriteDoubleAttribute(dbfhdl,entity,dbFldIDs[i], dbVals[i]);
		else
			ofile<<dbVals[i]<<_IFO_DEL;
	}
	for (int i = 0; i < strFldNum; i++)
	{
		if(!isTxt)
			success = DBFWriteStringAttribute(dbfhdl,entity,strFldIDs[i], strVals[i].c_str());
		else
			ofile<<strVals[i].c_str()<<_IFO_DEL;
	}
	for (int i = 0; i < logFldNum; i++)
	{
		if(!isTxt)
			success = DBFWriteStringAttribute(dbfhdl,entity,logFldIDs[i], logVals[i].c_str());
		else
			ofile<<logVals[i].c_str()<<_IFO_DEL;
	}
	ofile<<endl;
}
void writeIntersection(string file,int intFldNum, vector<vector<int> >& intfldVals, vector<string>& intFldNames,
					   int doubleFldNum, vector<vector<double> >& dbfldVals, vector<string>& dbFldNames,
					   int strFldNum, vector<vector<string> >& strfldVals, vector<string>& strFldNames,
					   int logFldNum, vector<vector<string> >& logfldVals, vector<string>& logFldNames,
					   GridIntersection* gridInter, WRITETYPE wtype, bool isTxt)
{
	//create a shape file
	ofstream ofile;
	SHPHandle handler = NULL;
	if(wtype == WRITE_POINT)
	{
		if(!isTxt)
			handler= SHPCreate(file.c_str(), SHPT_POINT );
		else
		{
			//ofile.open(file.c_str(),std::ofstream::out | std::ofstream::trunc); 
			ofile.open(file.c_str(), std::ofstream::out |std::ofstream::trunc); 
		}
	}
	else if(wtype == WRITE_POLYLINE)
	{
		if(!isTxt)
			handler = SHPCreate(file.c_str(), SHPT_ARC);
		else
		{
			ofile.open(file.c_str(),std::ofstream::out | std::ofstream::trunc); 
		}
	}
	else
	{
		if(!isTxt)
			handler = SHPCreate(file.c_str(), SHPT_POLYGON);
		else
		{
			ofile.open(file.c_str(),std::ofstream::out | std::ofstream::trunc); 
		}
	}

	//create a .dbf file
	DBFHandle dbfhdl = NULL;
	if(!isTxt)
	{
			string dbf = file;
			if( file.find('.')!=string::npos )
			{
				int dotidx = file.find_last_of('.');
				dbf = file.substr(0,dotidx);
			}
			dbf += ".dbf";
			dbfhdl = DBFCreate(dbf.c_str());
	}

	//create nodenumber field
	int ndfldId = -1;
	if(!isTxt)
		ndfldId = DBFAddField(dbfhdl, "nodenumber", FTInteger,12,0);
	vector<int> intFldIDs;//create integer fields
	vector<int> dbFldIDs;//create double fields
	vector<int> strFldIDs;//create string fields
	vector<int> logFldIDs;//log fields
	if(!isTxt)
	{
		addDBHandleFields(dbfhdl, intFldNum, intFldNames, doubleFldNum, dbFldNames, intFldIDs, dbFldIDs);
		addDBHandleStrFields(dbfhdl, strFldNum, strFldNames, strFldIDs);
		addDBHandleStrFields(dbfhdl, logFldNum, logFldNames, logFldIDs);
	}
	else
	{
			int success = 0;
			ofile<<"nodenumber"<<_IFO_DEL;
			for(int i = 0; i < intFldNum; ++i)
			{
					ofile<<intFldNames[i]<<_IFO_DEL;
			}
			for (int i = 0; i < doubleFldNum; ++i)
			{
					ofile<<dbFldNames[i]<<_IFO_DEL;
			}
			for (int i = 0; i < strFldNum; i++)
			{
					ofile<<strFldNames[i].c_str()<<_IFO_DEL;
			}
			for (int i = 0; i < logFldNum; i++)
			{
					ofile<<logFldNames[i].c_str()<<_IFO_DEL;
			}
			ofile<<endl;
	}

	vector<vector<int> >::iterator iit = intfldVals.begin();
	vector<vector<double> >::iterator dit = dbfldVals.begin();
	vector<vector<string> >::iterator sit = strfldVals.begin();
	vector<vector<string> >::iterator lit = logfldVals.begin();
	double* x = NULL, *y = NULL;
	//write points 
	if(wtype == WRITE_POINT)
	{
		x = new double[1]; y = new double[1];
		vector<PointCellIntersectionInfo>& pntInters = gridInter->getPointIntersectionInfo();
		vector<PointCellIntersectionInfo>::iterator pit = pntInters.begin();
		for(; pit != pntInters.end(); ++pit, ++iit, ++dit, ++sit, ++lit)
		{
			vector<int>& intVals = *iit;
			vector<double>& dbVals = *dit;
			vector<string>& strVals = *sit;
			vector<string>& logVals = *lit;
			x[0] = pit->pnt[0];
			y[0] = pit->pnt[1];
			
			SHPObject * shpObj = NULL;
			int entity = -1;
			if(!isTxt)
			{
				shpObj = SHPCreateSimpleObject(SHPT_POINT,1,x, y,NULL);
				entity = SHPWriteObject(handler, -1, shpObj);
			}

			//write nodenumber
			int success=  0;
			if(!isTxt)
				success = DBFWriteIntegerAttribute(dbfhdl, entity, ndfldId, pit->box->number == -1? -1 : pit->box->number + 1);
			else
				ofile<<(pit->box->number == -1? -1 : pit->box->number + 1)<<_IFO_DEL;
			writeIntDBVal(dbfhdl, entity, intFldIDs, intVals, dbFldIDs, dbVals, strFldIDs, strVals, logFldIDs, logVals, ofile, isTxt);
			if(!isTxt)
				SHPDestroyObject(shpObj);
		}
		delete[] x;
		delete[] y;
	}
	else if(wtype == WRITE_POLYLINE)
	{
		vector<PolylineCellIntersectionInfo>& plineInters = gridInter->getPolylineIntersectionInfo();
		vector<PolylineCellIntersectionInfo>::iterator pit = plineInters.begin();
		for(; pit != plineInters.end(); ++pit, ++iit, ++dit,++sit,++lit)
		{
			PolylineCellIntersectionInfo& lineInfo = *pit;
			vector<int>& intVals = *iit;
			vector<double>& dbVals = *dit;
			vector<string>& strVals = *sit;
			vector<string>& logVals = *lit;

			int vsize = lineInfo.m_plyline->getSize();
			if(vsize == 1)
			{
				cerr<<"! Warning: bad polyline, ignored!"<<endl;
				continue;
			}
			x = new double[vsize]; y = new double[vsize];
			ply_vertex* pv = lineInfo.m_plyline->getHead();
			x[0] = pv->getPos()[0]; y[0] = pv->getPos()[1];
			int idx = 0;
			do
			{
				pv = pv->getNext();
				idx++;
				x[idx] = pv->getPos()[0]; 
				y[idx] = pv->getPos()[1];
			}while(pv != lineInfo.m_plyline->getTail());
			
			SHPObject * shpObj = NULL;
			int entity = -1;
			if(!isTxt)
			{
					shpObj = SHPCreateSimpleObject(SHPT_ARC,vsize,x, y,NULL);
					entity = SHPWriteObject(handler, -1, shpObj);
			}

			//write nodenumber
			int success = 0;
			if(!isTxt)
				DBFWriteIntegerAttribute(dbfhdl, entity, ndfldId, pit->box->number == -1? -1 : pit->box->number + 1);
			else
				ofile<<(pit->box->number == -1? -1 : pit->box->number + 1)<<_IFO_DEL;
			writeIntDBVal(dbfhdl, entity, intFldIDs, intVals, dbFldIDs, dbVals, strFldIDs, strVals, logFldIDs, logVals, ofile, isTxt);

			if(!isTxt)
				SHPDestroyObject(shpObj);
			delete[] x;
			delete[] y;
		}
	}
	else if(wtype == WRITE_POLYGON)
	{
		vector<PolygonCellIntersectionInfo>& pgonInters = gridInter->getPolygonIntersectionInfo();
		vector<PolygonCellIntersectionInfo>::iterator pit = pgonInters.begin();
		for(int iShape = 0; pit != pgonInters.end(); ++pit, ++iit, ++dit, ++iShape,++sit,++lit)
		{
			PolygonCellIntersectionInfo& pgonInfo = *pit;
			vector<int>& intVals = *iit;
			vector<double>& dbVals = *dit;
			vector<string>& strVals = *sit;
			vector<string>& logVals = *lit;

			int partSize = pgonInfo.m_interPolygon->size();
			int vsize = pgonInfo.m_interPolygon->getSize();
			if(vsize <3)
			{
				cerr<<"! Warning: bad polygon, ignored!"<<endl;
				continue;
			}
			x = new double[vsize]; y = new double[vsize];

			//if(pgonInfo.m_interPolygon->size() > 1)
			//	cerr<<pgonInfo.m_interPolygon->size()<<_IFO_DEL;

			//make sure the out boundary 
			list<c_polygon::iterator> orderedPItrs;
			for(c_polygon::iterator plyit = pgonInfo.m_interPolygon->begin(); plyit != pgonInfo.m_interPolygon->end(); ++plyit)
			{
				if(plyit->getType() == c_ply::POUT)
					orderedPItrs.push_front(plyit);
				else
					orderedPItrs.push_back(plyit);
			}

			//part number
			int startPartIDAcumulate = 0;
			int* nParts = new int[partSize];
			memset(nParts, 0, sizeof(int)*partSize);
		    int k = 0, vid = 0;
			for(list<c_polygon::iterator>::iterator lit = orderedPItrs.begin(); lit != orderedPItrs.end(); ++lit, k++)
			{
				c_polygon::iterator& plyit = *lit;
				nParts[k] = startPartIDAcumulate;
				startPartIDAcumulate += plyit->getSize();

				ply_vertex* hv = plyit->getHead();
				ply_vertex* pv = hv;
				do
				{
					x[vid] = pv->getPos()[0];
					y[vid] = pv->getPos()[1];
					pv = pv->getNext();
					vid++;
				}while(pv != hv);
			}

			//create a polygon shape here
			//SHPObject * shpObj = SHPCreateObject(SHPT_POLYGON,vsize,x, y,NULL);
			SHPObject* shpObj = NULL;
			int entity = 0;
			if(!isTxt)
			{
					shpObj = SHPCreateObject(SHPT_POLYGON, iShape, partSize, nParts, NULL, vsize, x ,y, NULL, NULL); 
					entity = SHPWriteObject(handler, -1, shpObj);
			}
			//write nodenumber
			int success = 0;
			if(!isTxt)
				DBFWriteIntegerAttribute(dbfhdl, entity, ndfldId, pit->box->number == -1? -1 : pit->box->number + 1);
			else
				ofile<<(pit->box->number == -1? -1 : pit->box->number + 1)<<_IFO_DEL;
			writeIntDBVal(dbfhdl, entity, intFldIDs, intVals, dbFldIDs, dbVals, strFldIDs, strVals, logFldIDs, logVals, ofile, isTxt);

			if(!isTxt)
				SHPDestroyObject(shpObj);
			delete[] x;
			delete[] y;
			delete[] nParts;
		}
	}

	if(!isTxt)
	{
		DBFClose(dbfhdl);
		SHPClose(handler);
	}
	else
		ofile.close();
}

//write the grid and point intersection
void writePntIntersect(GridIntersection* gridInter, string file, vector<string>& attributes, FILETYPE filetype)
{
	//ofstream fout;
	//if(isTxt)
	//	fout.open(file.c_str(), std::ofstream::out | std::ofstream::trunc);
	vector<PointCellIntersectionInfo>& pntInters = gridInter->getPointIntersectionInfo();

	vector<string> intFldName, dbFldName, strFldName, logFldName;
	vector<vector<int> > intflds;
	vector<vector< double> > dbflds;
	vector<vector<string> > strflds;
	vector<vector<string> > logflds;
	//to record whether each attribute needs to be written
	vector<bool> bIntFld, bDBFld, bStrFld, bLogFld;

	//parse the attributes
	bool bpt = true;// checkAttributesContain("pointid", attributes);//write shapefile_pointid
	if(bpt)
		intFldName.push_back(string("pointid"));
	//add the attributes in shape file
	checkShpAttrContain(gridInter, bIntFld, bDBFld, bStrFld, bLogFld, intFldName,dbFldName,strFldName,logFldName,attributes);
	
	for(vector<PointCellIntersectionInfo>::iterator pit = pntInters.begin(); pit != pntInters.end(); ++pit)
	{
		PointCellIntersectionInfo& pci = *pit;
		vector<int> pints;
		vector<double> dbs;
		vector<string> strs;
		vector<string> logs;
		if(bpt)
			pints.push_back(pci.m_ptID);

		addShpAttrVal(pci.m_ptID, gridInter, bIntFld, bDBFld, bStrFld, bLogFld, pints, dbs, strs, logs);
		
		intflds.push_back(pints);
		dbflds.push_back(dbs);
		strflds.push_back(strs);
		logflds.push_back(logs);
	}

	if(filetype == WRITE_SHP || filetype == WRITE_TXT)
	{
		bool isTxt = (filetype == WRITE_TXT);
		writeIntersection(file, intFldName.size(), intflds, intFldName,dbFldName.size(), dbflds, dbFldName, 
			strFldName.size(), strflds,strFldName, logFldName.size(), logflds, logFldName,
			gridInter, WRITE_POINT, isTxt);
	}
	else
	{
	}
}

//write the line and grid intersection
void writeLineIntersect(GridIntersection* gridInter, string file, vector<string>& attributes, FILETYPE filetype)
{
	vector<PolylineCellIntersectionInfo>& pcti = gridInter->getPolylineIntersectionInfo();

	vector<string> intFldName, dbFldName, strFldName, logFldName;
	vector<vector<int> > intflds;
	vector<vector< double> > dbflds;
	vector<vector<string> > strflds;
	vector<vector<string> > logflds;
	//to record whether each attribute needs to be written
	vector<bool> bIntFld, bDBFld, bStrFld, bLogFld;

	bool barc = true;// checkAttributesContain("arcid", attributes);//write shapefile_arcid
	bool blen = true;//checkAttributesContain("length", attributes);//write length
	bool bstdist = true;//checkAttributesContain("starting_distance", attributes);//write starting_distance
	bool beddist = true;//checkAttributesContain("ending_distance", attributes);//write ending_distance

	if(barc)
		intFldName.push_back(string("arcid"));
	if(blen)
		dbFldName.push_back(string("length"));
	if(bstdist)
		dbFldName.push_back(string("starting_distance"));
	if(beddist)
		dbFldName.push_back(string("ending_distance"));
	//add the shape file attributes
	checkShpAttrContain(gridInter, bIntFld, bDBFld, bStrFld, bLogFld, intFldName,dbFldName,strFldName,logFldName,attributes);

	for (vector<PolylineCellIntersectionInfo>::iterator sit = pcti.begin(); sit != pcti.end(); ++sit)
	{
		PolylineCellIntersectionInfo& sio = *sit;
		vector<int> pints;
		vector<double> dbs;
		vector<string> strs;
		vector<string> logs;

		if(barc)
			pints.push_back(sio.m_curveID);
		if(blen)
			dbs.push_back(sio.m_len);
		if(bstdist)
			dbs.push_back(sio.m_startGeoDist);
		if(beddist)
			dbs.push_back(sio.m_endGeoDist);

		addShpAttrVal(sio.m_curveID, gridInter, bIntFld, bDBFld, bStrFld, bLogFld, pints, dbs, strs, logs);
		
		intflds.push_back(pints);
		dbflds.push_back(dbs);
		strflds.push_back(strs);
		logflds.push_back(logs);
	}

	if(filetype == WRITE_SHP || filetype == WRITE_TXT)
	{
		bool isTxt  = (filetype == WRITE_TXT);
		writeIntersection( file, intFldName.size(), intflds, intFldName,dbFldName.size(), dbflds, dbFldName, 
			strFldName.size(), strflds,strFldName, logFldName.size(), logflds, logFldName,  gridInter, WRITE_POLYLINE, isTxt);
	}
	else
	{
	}
}

//writhe the polygon and grid intersection
void writePolyIntersect(GridIntersection* gridInter, string file, vector<string>& attributes, FILETYPE filetype)
{
	vector<PolygonCellIntersectionInfo>& sci = gridInter->getPolygonIntersectionInfo();

	vector<string> intFldName, dbFldName, strFldName, logFldName;
	vector<vector<int> > intflds;
	vector<vector< double> > dbflds;
	vector<vector<string> > strflds;
	vector<vector<string> > logflds;
	//to record whether each attribute needs to be written
	vector<bool> bIntFld, bDBFld, bStrFld, bLogFld;

	bool bpoly = true;// checkAttributesContain("polyid", attributes);//write shapefile_polyid
	bool barea= true;//checkAttributesContain("totalarea", attributes);//write totalarea
	if(bpoly)
		intFldName.push_back(string("polyid"));
	if(barea)
		dbFldName.push_back(string("totalarea"));
	
	checkShpAttrContain(gridInter, bIntFld, bDBFld, bStrFld, bLogFld, intFldName,dbFldName,strFldName,logFldName,attributes);


	for (vector<PolygonCellIntersectionInfo>::iterator pit = sci.begin(); pit != sci.end(); ++pit)
	{
		PolygonCellIntersectionInfo& pio = *pit;
		vector<int> pints;
		vector<double> dbs;
		vector<string> strs;
		vector<string> logs;

		if(bpoly)
			pints.push_back(pio.m_polyID);
		if(barea)
			dbs.push_back(pio.m_area);

		addShpAttrVal(pio.m_polyID, gridInter, bIntFld, bDBFld, bStrFld, bLogFld, pints, dbs, strs, logs);
		
		intflds.push_back(pints);
		dbflds.push_back(dbs);
		strflds.push_back(strs);
		logflds.push_back(logs);
	}
	
	if(filetype == WRITE_SHP || filetype == WRITE_TXT)
	{
		bool isTxt  = (filetype == WRITE_TXT);
		writeIntersection(file, intFldName.size(), intflds, intFldName,dbFldName.size(), dbflds, dbFldName, 
			strFldName.size(), strflds,strFldName, logFldName.size(), logflds, logFldName, gridInter, WRITE_POLYGON, isTxt);
	}
	else
	{
	}
}

//gather the boxes
void gatherBoxes(GridIntersection* gridInter, vector<Box*>& boxes, WRITETYPE wtype)
{
	//gather the boxes
	vector<PointCellIntersectionInfo>& pntInters = gridInter->getPointIntersectionInfo();
	vector<PolylineCellIntersectionInfo>& pcti = gridInter->getPolylineIntersectionInfo();
	vector<PolygonCellIntersectionInfo>& sci = gridInter->getPolygonIntersectionInfo();
	if(wtype == WRITE_POINT)
	{
			for(vector<PointCellIntersectionInfo>::iterator pit = pntInters.begin(); pit != pntInters.end(); ++pit)
			{
				PointCellIntersectionInfo& pci = *pit;
				boxes.push_back(pci.box);
			}
	}
	else if(wtype == WRITE_POLYLINE)
	{
			for (vector<PolylineCellIntersectionInfo>::iterator sit = pcti.begin(); sit != pcti.end(); ++sit)
			{
				PolylineCellIntersectionInfo& sio = *sit;
				boxes.push_back(sio.box);
			}
	}
	else if(wtype == WRITE_POLYGON)
	{
			for (vector<PolygonCellIntersectionInfo>::iterator pit = sci.begin(); pit != sci.end(); ++pit)
			{
				PolygonCellIntersectionInfo& pio = *pit;
				boxes.push_back(pio.box);
			}
	}
}


//write the polygon and grid 
//this function omits the string type and log type fields
void writeIntersection2VTK(string vtkfile,int intFldNum, vector<vector<int> >& intfldVals, vector<string>& intFldNames,
					   int doubleFldNum, vector<vector<double> >& dbfldVals, vector<string>& dbFldNames,
					   GridIntersection* gridInter, WRITETYPE wtype)
{
	std::cerr<<"Warning: exporting VTK file would omit string type and logic type fields!\n";
	//write the points, loop over each layer
	stringstream ss0;
	ss0<<vtkfile.c_str();
	ss0<<".vtu";
	ofstream ofile(ss0.str().c_str());

	//gather boxes
	vector<Box*> boxes;
	gatherBoxes(gridInter, boxes,wtype);

	//todo: export the grid boxes to vtk file

	std::cerr<<"Exporting intersection to VTK file is not complete!\n;";

}



