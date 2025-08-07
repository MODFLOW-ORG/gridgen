
#include "shpReader.h"

#include <cassert>
using namespace std;

namespace cusg
{

typedef unsigned int uint;

void ShapeReader::destroy()
{
	m_ply_list.clear();
	m_arc_list.clear();
	m_ply_list.clear();

	m_intFlds.clear();
	m_dbFlds.clear();
	m_strFlds.clear();
	m_logFlds.clear();
}

//
//read field attributes from DBF file according to id
bool ShapeReader::readAttributes(DBFHandle& dbfHandle,int iShape)
{
	for(vector<ShapeIntField>::iterator iit = m_intFlds.begin(); iit != m_intFlds.end();  ++iit)
	{
		ShapeIntField& sifld= *iit;
		sifld.vals.insert(make_pair(iShape, DBFReadIntegerAttribute(dbfHandle, iShape, sifld.id)));
	}
	for(vector<ShapeDoubleField>::iterator iit = m_dbFlds.begin(); iit != m_dbFlds.end();  ++iit)
	{
		ShapeDoubleField& sifld= *iit;
		sifld.vals.insert(make_pair(iShape, DBFReadDoubleAttribute(dbfHandle, iShape, sifld.id)));
	}
	for(vector<ShapeStringField>::iterator iit = m_strFlds.begin(); iit != m_strFlds.end();  ++iit)
	{
		ShapeStringField& sifld= *iit;
		sifld.vals.insert(make_pair(iShape, string(DBFReadStringAttribute(dbfHandle, iShape, sifld.id))));
	}
	for(vector<ShapeLogField>::iterator iit = m_logFlds.begin(); iit != m_logFlds.end();  ++iit)
	{
		ShapeLogField& sifld= *iit;
		sifld.vals.insert(make_pair(iShape, string(DBFReadLogicalAttribute(dbfHandle, iShape, sifld.id))));
	}
	return true;
}


//
// read from files (both shp and db files)
//
bool ShapeReader::read(const string& filename)
{

    if(already_read) return true;

    //start to read
	SHPHandle	hSHP;
	hSHP = SHPOpen( filename.c_str(), "rb" );

	//something is wrong
    if( hSHP == NULL )
    {
		cout<<"- Error: Unable to open: "<<filename<<endl;
		return false;
    }



	/************************************************************/
#if 1	
	DBFHandle hDBF = DBFOpen( filename.c_str(), "rb" );

	//check
	if( hDBF == NULL )
	{
		cerr<<"! Error: DBFOpen("<< filename<<") failed."<<endl;
		return false;
	}

	//check
	if( DBFGetFieldCount(hDBF) == 0 )
	{
		cerr<<"! Warning: There are no fields in this table!"<<endl;
		return true;
	}

	////check more
	//if( (uint)DBFGetRecordCount(hDBF)!=m_ply_list.size()){
	//	cerr<<"! Error: Number of rows ("<<DBFGetRecordCount(hDBF)
	//		<<") in the database does not match to the number of buildings ("
	//		<<m_ply_list.size()<<")"<<endl;
	//}

	//used for reading from db
	int nWidth, nDecimals;

	//position of these fields in database...
	uint base_ele_pos=UINT_MAX;
	uint heigh_pos=UINT_MAX;
	for(int i = 0; i < DBFGetFieldCount(hDBF); i++ )
	{
		char    szTitle[12];
		DBFFieldType eType = DBFGetFieldInfo( hDBF, i, szTitle, &nWidth, &nDecimals );
		if(eType!=FTInteger && eType!=FTDouble && eType!=FTString && eType != FTLogical) continue; //not good...
		if(base_elevation_tag==szTitle) base_ele_pos=i;
		if(height_tag==szTitle) heigh_pos=i;

		//collect the field name and type info
				if(eType == FTInteger)
				{
					ShapeIntField sintfld;
					sintfld.id = i;
					sintfld.name  = szTitle;
					m_intFlds.push_back(sintfld);
				}
				else if(eType == FTDouble)
				{
					ShapeDoubleField sdbfld;
					sdbfld.id = i;
					sdbfld.name  = szTitle;
					m_dbFlds.push_back(sdbfld);
				}
				else if(eType == FTString)
				{
					ShapeStringField sstrfld;
					sstrfld.id = i;
					sstrfld.name  = szTitle;
					m_strFlds.push_back(sstrfld);
				}
				else if(eType == FTLogical)
				{
					ShapeLogField slogfld;
					slogfld.id = i;
					slogfld.name  = szTitle;
					m_logFlds.push_back(slogfld);
				}
		}

	////check if we can get the field ids for the information we wanted
	//if(base_ele_pos==UINT_MAX)
	//	cerr<<"! Warning: no base elevation found in the database"<<endl;

	//if(heigh_pos==UINT_MAX)
	//	cerr<<"! Warning: no height found in the database"<<endl;

#endif
	/************************************************************/
    
    //Get info
    int		nShapeType, nEntities;
    double 	adfMinBound[4], adfMaxBound[4];
    SHPGetInfo( hSHP, &nEntities, &nShapeType, adfMinBound, adfMaxBound );

    //reserve space
    m_ply_list.reserve(nEntities);

	//	Skim over the list of shapes, printing all the vertices
    for(uint i=0; i<(uint)nEntities; i++ )
    {
		
        //
        SHPObject * psShape = SHPReadObject( hSHP, i );

        //make sure psShape is not null
        assert(psShape);

		//some format check..
        if( psShape->nParts > 0 && psShape->panPartStart[0] != 0 )
        {
            cerr<<"panPartStart[0] = "<<psShape->panPartStart[0]
                <<", not zero as expected.\n";
        }

		//cerr<<"reading shape type : "<<psShape->nSHPType<<endl;
        switch(psShape->nSHPType)
        {
			case SHPT_POLYGON: 
			case SHPT_POLYGONZ:
			case SHPT_POLYGONM:
				{
					readPly(psShape, i); break;
				}
            case SHPT_ARC: 
			case SHPT_ARCZ:
			case SHPT_ARCM:
				{
					readArc(psShape, i); break;
				}
            case SHPT_POINT: 
			case SHPT_POINTM:
			case SHPT_POINTZ:
			case SHPT_MULTIPOINT:
			case SHPT_MULTIPOINTM:
			case SHPT_MULTIPOINTZ:
				{
					//cerr<<"reading shapefile with m/z \n";
					readPts(psShape, i); break;
				}
			case SHPT_MULTIPATCH:
				{
					readMultiPatch(psShape, i);
					break;
				}
            default:
				{
					cerr<<"! Error: ShapeReader::read: type unsupported: "<<psShape->nSHPType<<endl;
					break;
				}
        }

		//read the field values
		readAttributes(hDBF, i);

		 /*
        if( bValidate )
        {
            int nAltered = SHPRewindObject( hSHP, psShape );

            if( nAltered > 0 )
            {
                printf( "  %d rings wound in the wrong direction.\n",
                        nAltered );
                nInvalidCount++;
            }
        }
        */
        
        SHPDestroyObject( psShape );
    }

    SHPClose( hSHP );

#if 0
    cout.precision(20);
    cout<<"bbox[0]="<<bbox[0]<<" bbox[1]="<<bbox[1]<<" bbox[2]="<<bbox[2]<<" bbox[3]="<<bbox[3]<<endl;
    cout<<"width="<<bbox[1]-bbox[0]<<" height="<<bbox[3]-bbox[2]<<endl;
#endif

    //done readh shp file

    //
    //reading from database now.
    //
   /*
    DBFHandle hDBF = DBFOpen( filename.c_str(), "rb" );

    //check
    if( hDBF == NULL )
    {
        cerr<<"! Error: DBFOpen("<< filename<<") failed."<<endl;
        return false;
    }

    //check
    if( DBFGetFieldCount(hDBF) == 0 )
    {
        cerr<<"! Warning: There are no fields in this table!"<<endl;
        return true;
    }

    //check more
    if( (uint)DBFGetRecordCount(hDBF)!=m_ply_list.size()){
        cerr<<"! Error: Number of rows ("<<DBFGetRecordCount(hDBF)
            <<") in the database does not match to the number of buildings ("
            <<m_ply_list.size()<<")"<<endl;
    }

    //used for reading from db
    int nWidth, nDecimals;

    //position of these fields in database...
    uint base_ele_pos=UINT_MAX;
    uint heigh_pos=UINT_MAX;
    for(int i = 0; i < DBFGetFieldCount(hDBF); i++ )
    {
        char    szTitle[12];
        DBFFieldType eType = DBFGetFieldInfo( hDBF, i, szTitle, &nWidth, &nDecimals );
        if(eType!=FTInteger && eType!=FTDouble) continue; //not good...
        if(base_elevation_tag==szTitle) base_ele_pos=i;
        if(height_tag==szTitle) heigh_pos=i;
    }

    //check if we can get the field ids for the information we wanted
    if(base_ele_pos==UINT_MAX)
        cerr<<"! Warning: no base elevation found in the database"<<endl;

    if(heigh_pos==UINT_MAX)
        cerr<<"! Warning: no height found in the database"<<endl;

    DBFClose( hDBF );
    */

    already_read=true;
    return true;
}

void ShapeReader::readPly(SHPObject * psShape, int iShape)
{
    //read each part
   int iPart = 1;

   //
   c_ply ply(c_ply::POUT);
   ply.beginPoly();

   /*
   ply_extra_info& extra=ply.getExtraInfo();
   extra.box[0]=psShape->dfXMin;
   extra.box[1]=psShape->dfXMax;
   extra.box[2]=psShape->dfYMin;
   extra.box[3]=psShape->dfYMax;
   */

   for( int j = 0; j < psShape->nVertices; j++ )
   {
       string pszPartType;

       if( j == 0 && psShape->nParts > 0 )
           pszPartType = SHPPartTypeName( psShape->panPartType[0] );

       //this defines a new loop...
       if( iPart < psShape->nParts && psShape->panPartStart[iPart] == j )
       {
           pszPartType = SHPPartTypeName( psShape->panPartType[iPart] );
           iPart++;
           ply.endPoly();
           add_ply(ply, iShape);

           //create a new ply
           ply=c_ply(c_ply::PIN);
           ply.beginPoly();
       }

       //Points may duplicate....
       double& x=psShape->padfX[j];
       double& y=psShape->padfY[j];
       ply.addVertex(x,y);

       if(x<bbox[0]) bbox[0]=x;
       if(x>bbox[1]) bbox[1]=x;
       if(y<bbox[2]) bbox[2]=y;
       if(y>bbox[3]) bbox[3]=y;
   }

   ply.endPoly();

   add_ply(ply,iShape);
}
void ShapeReader::readMultiPatch(SHPObject* psShape, int iShape)
{
	for(int i = 0; i < psShape->nParts; i++)
	{
		int stId = psShape->panPartStart[i];
		int endId = psShape->nVertices - 1;
		if(i < psShape->nParts - 1)
		{
			endId = psShape->panPartStart[i + 1] - 1;
		}
		int type = psShape->panPartType[i];
		//cerr<<"current shape type : "<<type<<endl;
		if(type == SHPT_POLYGON || type == SHPT_POLYGONM || type == SHPT_POLYGONZ)
		{
			c_ply ply(c_ply::POUT);
			ply.beginPoly();
			for(int j = stId; j <= endId; j++)
			{
				//Points may duplicate....
				double& x=psShape->padfX[j];
				double& y=psShape->padfY[j];
				ply.addVertex(x,y);

				if(x<bbox[0]) bbox[0]=x;
				if(x>bbox[1]) bbox[1]=x;
				if(y<bbox[2]) bbox[2]=y;
				if(y>bbox[3]) bbox[3]=y;
			}
			ply.endPoly();
			 add_ply(ply, iShape);
		}
		else if(type == SHPT_ARC || type == SHPT_ARCZ || type == SHPT_ARCM)
		{
			GIS_plyline arc;
			arc.beginPoly();
			for(int j = stId; j <= endId; j++)
			{
				//Points may duplicate....
				double& x=psShape->padfX[j];
				double& y=psShape->padfY[j];
				arc.addVertex(x,y);

				if(x<bbox[0]) bbox[0]=x;
				if(x>bbox[1]) bbox[1]=x;
				if(y<bbox[2]) bbox[2]=y;
				if(y>bbox[3]) bbox[3]=y;
			}
			arc.endPoly();
			arc.m_id = iShape;
			//arc.m_id=m_arc_list.size();
			arc.computeGeodesicDist();
			m_arc_list.push_back(arc);
		}
		else if(type == SHPT_MULTIPOINT || type == SHPT_MULTIPOINTM || type == SHPT_MULTIPOINTZ)
		{
			for(int j = stId; j <= endId; j++)
			{
				//Points may duplicate....
				double& x=psShape->padfX[j];
				double& y=psShape->padfY[j];
				if(x<bbox[0]) bbox[0]=x;
				if(x>bbox[1]) bbox[1]=x;
				if(y<bbox[2]) bbox[2]=y;
				if(y>bbox[3]) bbox[3]=y;

				GIS_Point2d pt(x,y);
				//pt.m_id=m_pt_list.size();
				pt.m_id = iShape;
				m_pt_list.push_back(pt);
			}
		}
		else
		{
			cerr<<"error : unsupported type id: "<<type<<endl;
		}
	}

}

void ShapeReader::readArc(SHPObject * psShape, int iShape)
{
    //read each part
   int iPart = 1;

   //
   GIS_plyline arc;
   arc.beginPoly();

   for( int j = 0; j < psShape->nVertices; j++ )
   {
       string pszPartType;

       if( j == 0 && psShape->nParts > 0 )
           pszPartType = SHPPartTypeName( psShape->panPartType[0] );

       //this defines a new loop...
       if( iPart < psShape->nParts && psShape->panPartStart[iPart] == j )
       {
           pszPartType = SHPPartTypeName( psShape->panPartType[iPart] );
           iPart++;
           arc.endPoly();
           //arc.m_id=m_arc_list.size();
           arc.m_id = iShape;
		   arc.computeGeodesicDist();
           m_arc_list.push_back(arc);

           //create a new ply
           arc=GIS_plyline();
           arc.beginPoly();
       }

       //Points may duplicate....
       double& x=psShape->padfX[j];
       double& y=psShape->padfY[j];
       arc.addVertex(x,y);


       if(x<bbox[0]) bbox[0]=x;
       if(x>bbox[1]) bbox[1]=x;
       if(y<bbox[2]) bbox[2]=y;
       if(y>bbox[3]) bbox[3]=y;
   }

   arc.endPoly();
   arc.m_id = iShape;
   //arc.m_id=m_arc_list.size();
   arc.computeGeodesicDist();
   m_arc_list.push_back(arc);
}

void ShapeReader::readPts(SHPObject * psShape, int iShape)
{
   m_pt_list.reserve(psShape->nVertices);
   for( int j = 0; j < psShape->nVertices; j++ )
   {
       //Points may duplicate....
       double& x=psShape->padfX[j];
       double& y=psShape->padfY[j];
       if(x<bbox[0]) bbox[0]=x;
       if(x>bbox[1]) bbox[1]=x;
       if(y<bbox[2]) bbox[2]=y;
       if(y>bbox[3]) bbox[3]=y;

       GIS_Point2d pt(x,y);
       //pt.m_id=m_pt_list.size();
	   pt.m_id = iShape;
       m_pt_list.push_back(pt);
   }
}

void ShapeReader::add_ply(c_ply& ply, int iShape)
{
    if(ply.getType()==c_ply::POUT){
        //create a new polygon
        GIS_polygon polygon;
        //polygon.m_id=m_ply_list.size();
        polygon.m_id = iShape;
		polygon.push_back(ply);
        m_ply_list.push_back(polygon);
    }
    else{ //c_ply::PIN
        //add to the last polygon
        m_ply_list.back().push_back(ply);
    }
}

}//end namespace cusg
