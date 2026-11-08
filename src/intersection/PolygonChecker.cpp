#include "PolygonChecker.h"
#include <list>
#include <algorithm>
#include <vector>
#include <map>
using namespace std;
using namespace cusg;

//use the ring and area pair to record
struct RingArea
{
    Ring * ring;
    double area;
};
bool ring_area_greater(const RingArea & pa1, const RingArea & pa2)
{
    return pa1.area > pa2.area;
}
Ring & getOutBound(GIS_polygon & polygon)
{
    for (auto & git : polygon)
    {
        if (git.getType() == Ring::POUT) return git;
    }
    return polygon.front();
}

bool ring_encloses(const Ring & ring, const Point2d & point)
{
    //count the edges crossing a ray from the point; summing the angles they turn through
    //loses the small ones
    bool inside = false;
    LinkedVertex * vertex = ring.getHead();
    do
    {
        if (CrossesRay(vertex->getPos().get(), vertex->getNext()->getPos().get(), point.get()))
            inside = !inside;
        vertex = vertex->getNext();
    } while (vertex != ring.getHead());

    return inside;
}

//here we use the center of Ring to decide the containing relationship
bool checkHole(GIS_polygon & polygon, vector<GIS_polygon *> & polygons)
{

    vector<RingArea> rings_by_area;
    for (auto & ring : polygon)
    {
        RingArea entry;
        entry.ring = &ring;
        entry.area = fabs(ring.getArea());
        rings_by_area.push_back(entry);
    }
    //sort the ring based on the area descendingly
    sort(rings_by_area.begin(), rings_by_area.end(), ring_area_greater);

    //set the Ring having largest area
    bool all_in_first = true;
    rings_by_area.front().ring->set(Ring::POUT, rings_by_area.front().ring->getHead());
    GIS_polygon * outer_polygon = new GIS_polygon();
    outer_polygon->m_id = polygon.m_id;
    outer_polygon->push_back(*(rings_by_area.front().ring));

    polygons.push_back(outer_polygon);

    vector<RingArea>::iterator entry_it = rings_by_area.begin();
    ++entry_it;
    for (; entry_it != rings_by_area.end(); ++entry_it)
    {
        RingArea & entry = *entry_it;
        const Point2d & point = entry.ring->getCenter();
        if (ring_encloses(*(rings_by_area.front().ring), point))
        {
            entry.ring->set(Ring::PIN, entry.ring->getHead());
            outer_polygon->push_back(*(entry.ring));
        }
        else
        {
            all_in_first = false;
            bool enclosed = false;
            const Point2d & center = entry.ring->getCenter();
            for (auto & outer : polygons)
            {
                Ring & outer_ring = getOutBound(*outer);
                if (ring_encloses(outer_ring, center))
                {
                    entry.ring->set(Ring::PIN, entry.ring->getHead());
                    outer->push_back(*(entry.ring));
                    enclosed = true;
                    break;
                }
            }
            if (!enclosed)
            {
                GIS_polygon * new_polygon = new GIS_polygon();
                new_polygon->m_id = polygon.m_id;
                entry.ring->set(Ring::POUT, entry.ring->getHead());
                new_polygon->push_back(*(entry.ring));

                polygons.push_back(new_polygon);
            }
        }
    }

    return all_in_first;
}
//
////return true if only one polygon is in polygon
