
#include "polygon.h"
#include <vector>
#include "RingClipper.h"
using namespace cusg;
using namespace std;

//check the validity of the polygon, whether the holes are inside out-boundary
//return true: if the polygon is good
//return false: if
bool ring_encloses(const Ring & ring, const Point2d & point);

bool check_hole(GisPolygon & polygon, vector<GisPolygon *> & polygons);

////return true if only one polygon is in polygon