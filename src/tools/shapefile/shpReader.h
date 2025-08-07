
#ifndef _SHPREADER_H_
#define _SHPREADER_H_

#include "polygon.h"
#include "polyline.h"
#include "shapelib/shapefil.h"
#include <map>

using namespace std;

namespace cusg
{
struct ShapeIntField
{
	int id;
	string name;
	map<int, int> vals;
};
struct ShapeDoubleField
{
	int id;
	string name;
	map<int,double>vals; 
};
struct ShapeStringField
{
	int id;
	string name;
	map<int, string>vals; 
};
struct ShapeLogField
{
	int id;
	string name;
	map<int, string> vals;
};

// this reader reads shapefile format and convert 
// polygons to poly format

class ShapeReader 
{
public:

    ShapeReader(){ already_read=false; bbox[0]=bbox[2]=DBL_MAX; bbox[1]=bbox[3]=-DBL_MAX; }

	//read field attributes
	bool readAttributes(DBFHandle& dbfHanld,int iShape);

    //read from file
    bool read(const string& filename);

    //get poly/arc/pts
    vector<GIS_polygon>& getPolyList(){ return m_ply_list; }
    vector<GIS_plyline>& getArcList() { return m_arc_list; }
    vector<GIS_Point2d>& getPointList() { return m_pt_list; }

    //set tags for getting information from database
    void setBaseElevationTag(const string& tag) { base_elevation_tag=tag; }
    void setHeightTag(const string& tag) { height_tag=tag; }

	void destroy();

	vector<ShapeIntField>& getIntFields(){	return m_intFlds; } 
	vector<ShapeDoubleField>& getDBFields() { return m_dbFlds;	}
	vector<ShapeStringField>& getStrFields(){ return m_strFlds;	}
	vector<ShapeLogField>& getLogicalFields(){	return m_logFlds;	}

private:

    void add_ply(c_ply& plys, int iShape);

    void readPly(SHPObject * psShape, int iShape);
    void readArc(SHPObject * psShape, int iShape);
    void readPts(SHPObject * psShape, int iShape);
	void readMultiPatch(SHPObject* psShape, int iShape);

    vector<GIS_polygon> m_ply_list; //a list of polygons
    vector<GIS_plyline> m_arc_list; //a list of arcs
    vector<GIS_Point2d> m_pt_list;  //a list of points
    string base_elevation_tag;
    string height_tag;

    double bbox[4];

    bool already_read;

	////add by Guilin
	//vector<int> m_pt_id;
	//vector<int> m_arc_id;
	//vector<int> m_ply_id;

	//field id, names and values
	vector<ShapeIntField> m_intFlds;
	vector<ShapeDoubleField> m_dbFlds;
	vector<ShapeStringField> m_strFlds;
	vector<ShapeLogField> m_logFlds;
};

}//namespace cusg

#endif //_SHPREADER_H_
