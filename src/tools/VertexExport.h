/*! \file VertexExport.h
\brief files that describe the cells of a grid by their shared vertices
*/
#pragma once

#include <string>

#include "Grid.h"

namespace cusg
{

//write a MODFLOW-USG grid specification file with the corners of each active
//cell at its top and bottom, and the vertices of finer neighbors on its edges
//if hanging_vertices is true
void write_gsf(Grid * grid, const std::string & filename, bool hanging_vertices);

//write the node number of each active cell with the layer, row, and column of
//the base grid cell that contains it
void write_node_data_table(Grid * grid, const std::string & filename);

//write the vertices and the cell2d records of each active cell, for the
//VERTICES and CELL2D blocks of MODFLOW 6 DISU and DISV input
void write_vertices_cell2d(Grid * grid, const std::string & prefix);

} // namespace cusg
