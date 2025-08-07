
#pragma once

#include "Grid.h"
#include "ModflowGrid2D.h"
#include "ModflowGrid.h"
#include "QuadTree3D.h"
#include <list>

using namespace std;
using namespace cusg;


struct VPTPair
{
	VPTPair()
	{
		bSrc = NULL;
		bDst = NULL;
		area = 0.0;
	}
	VPTPair(Box* b1, Box* b2, double d)
	{
		bSrc = b1; bDst = b2;
		area = d;
	}
	VPTPair(const VPTPair& other)
	{
		bSrc = other.bSrc;
		bDst = other.bDst;
		area = other.area;
	}
	Box* bSrc;
	Box* bDst;
	double area;
};


struct VPassBox
{
	VPassBox()
	{
		bSrc = NULL;
		bDst = NULL;
		area = 0.0;
	}
	VPassBox(const VPassBox& other)
	{
		bSrc= other.bSrc;
		bDst = other.bDst;
		bPssList.clear();
		bPssList.insert(bPssList.end(), other.bPssList.begin(), other.bPssList.end());
		area = other.area;
	}

	Box* bSrc;
	Box* bDst;
	list<Box*> bPssList;
	double area;
};

//void boxVPConnection(Box* b, vector<VPassBox>& vpbl, bool topToDown);

void getVPT(Box* b, vector<VPTPair>& vpts, bool topToDown);
