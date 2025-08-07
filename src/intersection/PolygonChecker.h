
#include "polygon.h"
#include <vector>
#include "BoxPlyClip.h"
using namespace cusg;
using namespace std;

//check the validility of the polygon, whether the holes are inside out-boundary
//return true: if the polygon is good
//return false: if 
bool CPlyEnclose(const c_ply& cply, const  Point2d& p2d);

bool checkHole(GIS_polygon& plygon, vector<GIS_polygon*>& plygons);

bool CPolygonEnclose(const c_polygon& cply, const Point2d& p2d);

////return true if only one polygon is in plygon
//bool preProPolygon(GIS_polygon& plygon, vector<GIS_polygon*>& gplygons);