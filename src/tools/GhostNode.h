#ifndef _GHOSTNODE_H_
#define _GHOSTNODE_H_

#include "QuadTree3D.h"
namespace cusg
{
	//n (node number of big cell),
	//	m (node number of the small cell),
	//	j1 (neighbor in the direction of the ghost node),
	//	j2 (neighbor in the direction of the ghost node),
	//	alpha_1 (parameterized distance to the ghost node from j1),
	//	alpha_2 (parameterized distance to the ghost node from j2)
	//	(if j1 and j2 are defined then alpha_1 and alpha_2 are
	//	half of the case of only j1 is defined)
	struct GhostNode
	{
		long nNum;
		long mNum;
		long j1;
		long j2;
		double apha1;
		double apha2;
	};


	//compute the ghost node for QuadTree3D grid
	void computeGNC(QuadTree3D* quadGrid, vector<GhostNode>& gns);

	//write ghost node
	void writeGhostNode(string filePath, vector<GhostNode>& gns, int offset);
}

#endif