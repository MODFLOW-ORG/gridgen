#pragma once
#ifndef DBF_FIELDS_H
#define DBF_FIELDS_H

//width of the integer fields of the .dbf files gridgen writes: a numeric field with no decimals
//that is wider than 10 is read as a double by shapelib and ArcGIS (#188), and 10 digits hold
//any node number, layer, row, column, or feature id
constexpr int DBF_INTEGER_WIDTH = 10;

#endif
