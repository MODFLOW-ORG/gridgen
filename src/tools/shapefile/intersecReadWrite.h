
#ifndef _INTERSECT_READ_WRITE_H_
#define _INTERSECT_READ_WRITE_H_

#include "Box.h"
#include "polygon.h"
#include "polyline.h"
#include "shapelib/shapefil.h"
#include "GridIntersection.h"
#include "shpReader.h"
#include <string.h>
#include <vector>
using namespace cusg;
using namespace std;

enum WRITETYPE
{
    WRITE_POINT,
    WRITE_POLYLINE,
    WRITE_POLYGON
};

enum FILETYPE
{
    WRITE_TXT,
    WRITE_SHP,
    WRITE_VTK
};

//check whether the array attributes contain the given attribute
bool checkAttributesContain(const char * attribute_name, vector<string> & attributes);

//write the grid and point intersection
void writePntIntersect(GridIntersection * grid_intersection, string file,
                       vector<string> & attributes, FILETYPE filetype);

//write the line and grid intersection
void writeLineIntersect(GridIntersection * grid_intersection, string file,
                        vector<string> & attributes, FILETYPE filetype);

//writhe the polygon and grid intersection
void writePolyIntersect(GridIntersection * grid_intersection, string file,
                        vector<string> & attributes, FILETYPE filetype);

#endif