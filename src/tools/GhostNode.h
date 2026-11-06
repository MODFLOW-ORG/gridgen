#ifndef _GHOSTNODE_H_
#define _GHOSTNODE_H_

#include "QuadTree3D.h"
namespace cusg
{
//large_node (node number of the large cell),
//	small_node (node number of the small cell),
//	j1 (neighbor in the direction of the ghost node),
//	j2 (neighbor in the direction of the ghost node),
//	alpha1 (parameterized distance to the ghost node from j1),
//	alpha2 (parameterized distance to the ghost node from j2)
//	(if j1 and j2 are defined then alpha1 and alpha2 are
//	half of the case of only j1 is defined)
struct GhostNode
{
    long large_node;
    long small_node;
    long j1;
    long j2;
    double alpha1;
    double alpha2;
};

//compute the ghost node for QuadTree3D grid
void computeGNC(QuadTree3D * qtree, vector<GhostNode> & ghost_nodes);

//write ghost node
void writeGhostNode(string filePath, vector<GhostNode> & ghost_nodes, int offset);
}

#endif