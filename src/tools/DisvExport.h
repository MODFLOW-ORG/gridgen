/*! \file DisvExport.h
\brief MODFLOW 6 DISV input for a grid with the same cells in every layer
*/
#pragma once

#include <string>
#include <vector>

#include "Grid.h"

namespace cusg
{

//a cell of a layer of DISV input: its center and the 0-based indices of its
//vertices, clockwise
struct Cell2d
{
    double x;
    double y;
    std::vector<int> vertices;
};

//a horizontal ghost node by 0-based layer and cell of the layer
struct DisvGhostNode
{
    int layer;
    int large_cell;
    int small_cell;
    int contributor1;
    int contributor2;
    double alpha1;
    double alpha2;
};

//MODFLOW 6 DISV input; the arrays of every layer are in the order of the cells
//of a layer, ncpl values for each layer
struct DisvData
{
    int nlay = 0;
    int ncpl = 0;
    std::vector<Point2d> vertices;
    std::vector<Cell2d> cell2d;
    std::vector<double> top;
    std::vector<double> botm;
    std::vector<int> idomain;
    //true for a quadtree grid, whose ghost nodes are written even if there
    //are none
    bool has_ghost_nodes = false;
    std::vector<DisvGhostNode> ghost_nodes;
    //connections from a cell next to an inactive cell to a smaller cell whose
    //ghost node is missing or does not reproduce a linear head field
    int inexact_ghost_nodes = 0;
    //active cells whose top differs from the bottom of the cell above, which
    //DISV takes as their top
    int tops_not_aligned = 0;
};

//the DISV input of a MODFLOW grid, or of a quadtree grid with the same cells in
//every layer; inactive cells are -1 in IDOMAIN if vertical_pass_through is
//true and 0 otherwise. Throws GridgenError if a layer of a quadtree grid is
//refined differently from the first
DisvData disv_data(Grid * grid, bool vertical_pass_through);

//write the DISV input to files named for prefix: dimensions, vertices, cell2d,
//top, botm, idomain, and, if there are ghost nodes, gnc
void write_disv_data(const DisvData & data, const std::string & prefix);

} // namespace cusg
