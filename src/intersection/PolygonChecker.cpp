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
    c_ply * ply;
    double area;
};
bool compPlyAreaGreater(const PlyAreaPair & pa1, const PlyAreaPair & pa2)
{
    return pa1.area > pa2.area;
}
c_ply & getOutBound(GIS_polygon & polygon)
{
    for (auto & git : polygon)
    {
        if (git.getType() == c_ply::POUT) return git;
    }
    return polygon.front();
}

bool CPlyEnclose(const c_ply & ply, const Point2d & point)
{
    //count the edges crossing a ray from the point; summing the angles they turn through
    //loses the small ones
    bool inside = false;
    ply_vertex * vertex = ply.getHead();
    do
    {
        if (CrossesRay(vertex->getPos().get(), vertex->getNext()->getPos().get(), point.get()))
            inside = !inside;
        vertex = vertex->getNext();
    } while (vertex != ply.getHead());

    return inside;
}

//here we use the center of c_ply to decide the containing relationship
bool checkHole(GIS_polygon & polygon, vector<GIS_polygon *> & polygons)
{

    vector<PlyAreaPair> rings_by_area;
    for (auto & ply : polygon)
    {
        PlyAreaPair pap;
        pap.ply = &ply;
        pap.area = fabs(ply.getArea());
        rings_by_area.push_back(pap);
    }
    //sort the ply based on the area descendingly
    sort(rings_by_area.begin(), rings_by_area.end(), compPlyAreaGreater);

    //set the c_ply having largest area
    bool all_in_first = true;
    rings_by_area.front().ply->set(c_ply::POUT, rings_by_area.front().ply->getHead());
    GIS_polygon * outer_polygon = new GIS_polygon();
    outer_polygon->m_id = polygon.m_id;
    outer_polygon->push_back(*(rings_by_area.front().ply));

    polygons.push_back(outer_polygon);

    vector<PlyAreaPair>::iterator pit = rings_by_area.begin();
    ++pit;
    for (; pit != rings_by_area.end(); ++pit)
    {
        PlyAreaPair & pap = *pit;
        const Point2d & point = pap.ply->getCenter();
        if (CPlyEnclose(*(rings_by_area.front().ply), point))
        {
            pap.ply->set(c_ply::PIN, pap.ply->getHead());
            outer_polygon->push_back(*(pap.ply));
        }
        else
        {
            all_in_first = false;
            bool enclosed = false;
            const Point2d & center = pap.ply->getCenter();
            for (auto & outer : polygons)
            {
                c_ply & outer_ring = getOutBound(*outer);
                if (CPlyEnclose(outer_ring, center))
                {
                    pap.ply->set(c_ply::PIN, pap.ply->getHead());
                    outer->push_back(*(pap.ply));
                    enclosed = true;
                    break;
                }
            }
            if (!enclosed)
            {
                GIS_polygon * new_polygon = new GIS_polygon();
                new_polygon->m_id = polygon.m_id;
                pap.ply->set(c_ply::POUT, pap.ply->getHead());
                new_polygon->push_back(*(pap.ply));

                polygons.push_back(new_polygon);
            }
        }
    }

    return all_in_first;
}
//
////return true if only one polygon is in polygon
