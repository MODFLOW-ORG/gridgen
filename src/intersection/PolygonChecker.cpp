#include "PolygonChecker.h"
#include <list>
#include <algorithm>
#include <vector>
#include <map>
using namespace std;
using namespace cusg;

//use the ply and area pair to record
struct PlyAreaPair
{
    c_ply * cply;
    double area;
};
bool compPlyAreaGreater(const PlyAreaPair & pa1, const PlyAreaPair & pa2)
{
    return pa1.area > pa2.area;
}
c_ply & getOutBound(GIS_polygon & gp)
{
    for (auto & git : gp)
    {
        if (git.getType() == c_ply::POUT) return git;
    }
    return gp.front();
}

bool CPlyEnclose(const c_ply & cply, const Point2d & p2d)
{
    //count the edges crossing a ray from the point; summing the angles they turn through
    //loses the small ones
    bool inside = false;
    ply_vertex * tpv = cply.getHead();
    do
    {
        if (CrossesRay(tpv->getPos().get(), tpv->getNext()->getPos().get(), p2d.get()))
            inside = !inside;
        tpv = tpv->getNext();
    } while (tpv != cply.getHead());

    return inside;
}

//here we use the center of c_ply to decide the containing relationship
bool checkHole(GIS_polygon & plygon, vector<GIS_polygon *> & plygons)
{

    vector<PlyAreaPair> plyareas;
    for (auto & ply : plygon)
    {
        PlyAreaPair pap;
        pap.cply = &ply;
        pap.area = fabs(ply.getArea());
        plyareas.push_back(pap);
    }
    //sort the ply based on the area descendingly
    sort(plyareas.begin(), plyareas.end(), compPlyAreaGreater);

    //set the c_ply having largest area
    bool allInPlygon = true;
    plyareas.front().cply->set(c_ply::POUT, plyareas.front().cply->getHead());
    GIS_polygon * ogp = new GIS_polygon();
    ogp->m_id = plygon.m_id;
    ogp->push_back(*(plyareas.front().cply));

    plygons.push_back(ogp);

    vector<PlyAreaPair>::iterator pit = plyareas.begin();
    ++pit;
    for (; pit != plyareas.end(); ++pit)
    {
        PlyAreaPair & pap = *pit;
        const Point2d & p2d = pap.cply->getCenter();
        if (CPlyEnclose(*(plyareas.front().cply), p2d))
        {
            pap.cply->set(c_ply::PIN, pap.cply->getHead());
            ogp->push_back(*(pap.cply));
        }
        else
        {
            allInPlygon = false;
            bool b = false;
            const Point2d & pad = pap.cply->getCenter();
            for (auto & outer : plygons)
            {
                c_ply & outp = getOutBound(*outer);
                if (CPlyEnclose(outp, pad))
                {
                    pap.cply->set(c_ply::PIN, pap.cply->getHead());
                    outer->push_back(*(pap.cply));
                    b = true;
                    break;
                }
            }
            if (!b)
            {
                GIS_polygon * gp = new GIS_polygon();
                gp->m_id = plygon.m_id;
                pap.cply->set(c_ply::POUT, pap.cply->getHead());
                gp->push_back(*(pap.cply));

                plygons.push_back(gp);
            }
        }
    }

    return allInPlygon;
}
//
////return true if only one polygon is in plygon
