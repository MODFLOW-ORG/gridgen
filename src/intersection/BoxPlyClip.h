#pragma once
#include "intersection.h"
#include "polygon.h"
#include <list>
#include <vector>

using namespace std;

#define EPSION 0.0001
#define SMALLCOORDINATE 1e-8
#define SMALLAREA 1e-7

namespace cusg
{
class BoxPlyClip;

struct VNode
{
	double x,y;
	struct VNode * next;
	struct VNode * prev;
	bool isInDirect;
	bool isOnline;
	bool isVisited;
	bool isScanPnt;
	VNode()
	{
		init();
	}
	VNode(double tx, double ty)
	{
		x =tx; y =ty;
		init();
	}
	~VNode()
	{
		next = NULL;
		prev = NULL;
	}
	void init()
	{
		next = NULL;
		prev = NULL;
		isInDirect = false;
		isOnline = false;
		isVisited =false;
		isScanPnt = false;
	}
	void resetStatus()
	{
		isInDirect = false;
		isOnline = false;
		isVisited =false;
		isScanPnt = false;
	}
	void setDirection(double fval, bool isX, bool isGreater)
	{
		if(isX)
		{
			if(fabs(x - fval) <SMALLNUMBER)
			{
				//isScanPnt = isOnline = true;
				isInDirect = true;return;
			}
			//bool big = (x >= fval);
			bool big = (x >= fval);
			isInDirect = (big == isGreater);
		}
		else
		{
			if(fabs(y - fval) < SMALLNUMBER)
			{
				//isScanPnt = isOnline = true;
				isInDirect = true; return;
			}
			//bool big = (y >= fval);
			bool big = (y >= fval);
			isInDirect = (big == isGreater);
		}
	}
	bool Online(double fval, bool isX)
	{
		if(isX && (fabs(x - fval) < SMALLNUMBER))
			return true;
		else if( (!isX) && (fabs(y - fval) < SMALLNUMBER) )
			return true;
		return false;
	}

};

class BoxPlyClip
{
public:
	BoxPlyClip(void);
	~BoxPlyClip(void);

	void destroy();
	void destroy(VNode* vn);
	void resetStatus();

	//add a vertex
	void add(VNode* v);
	//insert a vertex
	void insert(VNode* s, VNode* e, VNode* v);
	//next
	VNode* next(VNode* v);

	//clip using scanning line
	void clip(double fval, bool isX, bool isBigger, list<BoxPlyClip*>& cres);
	//clip using a box
	void clip(double xmin, double xmax, double ymin, double ymax, list<BoxPlyClip*>& cres);


	//set intersection based on the scanning line
	int setIntesection(VNode* startV, double fval, bool isX, bool isBigger, vector<VNode*>& intersectV);

	//collect 
	void collectPolygonsInInterval(bool isX, bool isBigger, vector<VNode*>& interVs, list<BoxPlyClip*>& cres, int from, int to, bool isNextO);

	//collect polygon results
	void collectPolygons(double fval, bool isX, bool isBigger, vector<int>& voNums, vector<VNode*>& interVs, list<BoxPlyClip*>& cres);
	//collect one ply
	//BoxPlyClip* collectOnePly(vector<VNode*>& interVs, int sIdx, int eIdx, bool isNextOrder);
	void collectOnePly(vector<VNode*>& interVs, int sIdx, int eIdx, bool isNextOrder, list<BoxPlyClip*>& bpcResults,bool isX, bool isBigger);

	void collectSecondPart(vector<VNode*>& interVs, int sIdx, int tevIdx,VNode* sv, VNode* tev, bool isNextOrder, int vsize, list<BoxPlyClip*>& bpcResults,
		bool isX, bool isBigger, vector<VNode*>& secondHalf);

	void constructClipPly(vector<VNode*>& firstHalf, vector<VNode*>& secondHalf, list<BoxPlyClip*>& bpcResults);

	VNode* collectAlongWay(vector<VNode*>& vs, VNode* tmpv, bool nextO);
	//unprocessed
	//bool unprocessed();
	VNode* copy(const VNode* v);

	//this vector stores the array of holes
	VNode* first;//mark the outboundary
	list<VNode*> plys;//make sure the order, out boundary is the first one

	//new added vertices
	list<VNode*> newAdds;
	//influenced vertices
	list<VNode*> influVs;
};
//convert the ply to BoxPlyClip
BoxPlyClip* convertPlyToClip(const c_ply& tply);

//convert BoxPlyClip to c_ply
//void convertClipToPly(BoxPlyClip* tcp, c_ply& resPly);

void convertClipToPly(VNode* fv, c_ply& resPly);

BoxPlyClip* convertCPolygonToClip(const GIS_polygon& polygon);

BoxPlyClip* convertSimplePolygonToClip(const c_polygon& polygon);

void convertClipToPolygon(BoxPlyClip* bpc, GIS_polygon* cpolygon);

VNode* copyOnePly(VNode* head, bool bSetVisited = false);

double getRotAngle(Vector2d& v1, Vector2d& v2);

bool compareTwoBPC(BoxPlyClip* b1, BoxPlyClip* b2);

//make sure the min and max are in correct order
double intersectTwoBox(Point2d ld, Point2d rd, Point2d ru, Point2d lu, 
	double xmin2, double xmax2, double ymin2, double ymax2);


}