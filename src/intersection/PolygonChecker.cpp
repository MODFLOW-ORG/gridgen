#include "PolygonChecker.h"
#include <list>
#include <algorithm>
#include <vector>
#include <map>
using namespace std;
using namespace cusg;

#define  TESTING 0

//use the ply and area pair to record
struct PlyAreaPair
{
	c_ply* cply;
	double area;
};
 bool compPlyAreaGreater(const PlyAreaPair& pa1, const PlyAreaPair& pa2)
 {
	 return pa1.area > pa2.area;
 }
 c_ply& getOutBound(GIS_polygon& gp)
 {
	 for(GIS_polygon::iterator git = gp.begin(); git != gp.end(); ++git)
	 {
		 if(git->getType() == c_ply::POUT)
			 return *git;
	 }
	 return gp.front();
 }

bool CPlyEnclose(const c_ply& cply, const Point2d& p2d)
{
	double wnum = 0.0;
	//if(cply.getSize() < 3)
	//	return false;
	ply_vertex* tpv = cply.getHead();
	do{
		Vector2d vc1(tpv->getPos()[0] - p2d[0], tpv->getPos()[1] - p2d[1]);
		Vector2d vc2(tpv->getNext()->getPos()[0] - p2d[0], tpv->getNext()->getPos()[1] - p2d[1]);
		double ang = getRotAngle(vc1, vc2);
		if(ang > PI)
			ang = ang - 2 * PI;
		wnum += ang;

		tpv = tpv->getNext();
	}while(tpv != cply.getHead());

	if(fabs(wnum - 2* PI)< EPSION || fabs(wnum + 2*PI) < EPSION)
		return true;

	return false;
}

bool CPolygonEnclose(const c_polygon& cplygon,const Point2d& p2d)
{
	bool inside=true;
	for(c_polygon::const_iterator pit = cplygon.begin(); pit != cplygon.end(); pit++)
	{
		const c_ply& cp = *pit;
		if(cp.getType()==c_ply::POUT)
		{
			if(!CPlyEnclose(cp, p2d))
				return false;
		}
		else
		{
			if(CPlyEnclose(cp,p2d))
				return false;
		}
	}
	return true;
}

//here we use the center ot c_ply to decide the containing relationship
bool checkHole(GIS_polygon& plygon, vector<GIS_polygon*>& plygons)
{
	//if(plygon.size() < 2)
	//	return true;
#if TESTING
	int id = 0;
#endif


	vector<PlyAreaPair> plyareas;
	for(GIS_polygon::iterator git = plygon.begin(); git != plygon.end(); ++git)
	{
		PlyAreaPair pap;
		pap.cply = &(*git);
		pap.area = fabs(git->getArea());
		plyareas.push_back(pap);
	}
	//sort the ply based on the area descendingly
	sort(plyareas.begin(), plyareas.end(), compPlyAreaGreater);

	//set the c_ply having largest area
	bool allInPlygon = true;
	plyareas.front().cply->set(c_ply::POUT, plyareas.front().cply->getHead());
	GIS_polygon* ogp = new GIS_polygon();
	ogp->m_id = plygon.m_id;
	ogp->push_back(*(plyareas.front().cply));

#if TESTING
	id ++;
	ogp->m_id = id;
#endif

	plygons.push_back(ogp);



	vector<PlyAreaPair>::iterator pit = plyareas.begin();
	++pit;
	for(; pit != plyareas.end(); ++pit)
	{
		PlyAreaPair& pap = *pit;
		const Point2d& p2d = pap.cply->getCenter();
		if(CPlyEnclose(*(plyareas.front().cply), p2d))
		{
			pap.cply->set(c_ply::PIN, pap.cply->getHead());
			ogp->push_back(*(pap.cply));
		}
		else
		{
			allInPlygon = false;
			bool b = false;
			const Point2d& pad=pap.cply->getCenter();
			for(vector<GIS_polygon*>::iterator vit =plygons.begin(); vit != plygons.end(); ++vit)
			{
				c_ply& outp = getOutBound(*(*vit));
				if(CPlyEnclose(outp, pad))
				{
					pap.cply->set(c_ply::PIN, pap.cply->getHead());
					(*vit)->push_back(*(pap.cply));
					b = true; break;
				}
			}
			if(!b)
			{
				GIS_polygon* gp = new GIS_polygon();
				gp->m_id = plygon.m_id;
				pap.cply->set(c_ply::POUT, pap.cply->getHead());
				gp->push_back(*(pap.cply));

#if TESTING
				id ++;
				gp->m_id = id;
#endif
				plygons.push_back(gp);


			}

		}

	}
		
		return allInPlygon;
}
//
////return true if only one polygon is in plygon
//bool preProPolygon(GIS_polygon& plygon, vector<GIS_polygon*>& gplygons)
//{
//
//
//	return true;
//}

