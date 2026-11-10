#ifndef GHOSTNODE_H
#define GHOSTNODE_H

#include "QuadTree3D.h"
namespace cusg
{
//large_node (node number of the large cell),
//	small_node (node number of the small cell),
//	j1 (neighbor of the large cell the ghost node is interpolated from),
//	j2 (second neighbor the ghost node is interpolated from),
//	alpha1 (contributing factor of j1),
//	alpha2 (contributing factor of j2)
//	(when only one cell contributes, j1 and j2 are that cell and
//	alpha1 and alpha2 are each half its factor)
struct GhostNode
{
    long large_node;
    long small_node;
    long j1;
    long j2;
    double alpha1;
    double alpha2;
};

//compute the ghost nodes of the horizontal connections between cells of
//different size, and of the vertical connections too if vertical is true
void compute_ghost_nodes(QuadTree3D * qtree, vector<GhostNode> & ghost_nodes,
                         bool vertical = false);

//write ghost node
void write_ghost_node(string file_path, vector<GhostNode> & ghost_nodes, int offset);
}

#endif