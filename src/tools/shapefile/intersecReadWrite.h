
#ifndef INTERSECT_READ_WRITE_H
#define INTERSECT_READ_WRITE_H

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
bool check_attributes_contain(const char * attribute_name, vector<string> & attributes);

//write the grid and point intersection
void write_pnt_intersect(GridIntersection * grid_intersection, string file,
                         vector<string> & attributes, FILETYPE filetype);

//write the line and grid intersection
void write_line_intersect(GridIntersection * grid_intersection, string file,
                          vector<string> & attributes, FILETYPE filetype);

//writhe the polygon and grid intersection
void write_poly_intersect(GridIntersection * grid_intersection, string file,
                          vector<string> & attributes, FILETYPE filetype);

#endif