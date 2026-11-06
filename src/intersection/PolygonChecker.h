
#include "polygon.h"
#include <vector>
#include "BoxPlyClip.h"
using namespace cusg;
using namespace std;

//check the validility of the polygon, whether the holes are inside out-boundary
//return true: if the polygon is good
//return false: if
bool CPlyEnclose(const c_ply & ply, const Point2d & point);

bool checkHole(GIS_polygon & polygon, vector<GIS_polygon *> & polygons);

////return true if only one polygon is in plygon