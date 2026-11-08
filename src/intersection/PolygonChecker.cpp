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
Ring & get_out_bound(GisPolygon & polygon)
{
    for (auto & git : polygon)
    {
        if (git.get_type() == Ring::POUT) return git;
    }
    return polygon.front();
}

bool ring_encloses(const Ring & ring, const Point2d & point)
{
    //count the edges crossing a ray from the point; summing the angles they turn through
    //loses the small ones
    bool inside = false;
    LinkedVertex * vertex = ring.get_head();
    do
    {
        if (crosses_ray(vertex->get_pos().get(), vertex->get_next()->get_pos().get(), point.get()))
            inside = !inside;
        vertex = vertex->get_next();
    } while (vertex != ring.get_head());

    return inside;
}

//here we use the center of Ring to decide the containing relationship
bool check_hole(GisPolygon & polygon, vector<GisPolygon *> & polygons)
{

    vector<RingArea> rings_by_area;
    for (auto & ring : polygon)
    {
        RingArea entry;
        entry.ring = &ring;
        entry.area = fabs(ring.get_area());
        rings_by_area.push_back(entry);
    }
    //sort the ring based on the area descendingly
    sort(rings_by_area.begin(), rings_by_area.end(), ring_area_greater);

    //set the Ring having largest area
    bool all_in_first = true;
    rings_by_area.front().ring->set(Ring::POUT, rings_by_area.front().ring->get_head());
    GisPolygon * outer_polygon = new GisPolygon();
    outer_polygon->m_id = polygon.m_id;
    outer_polygon->push_back(*(rings_by_area.front().ring));

    polygons.push_back(outer_polygon);

    vector<RingArea>::iterator entry_it = rings_by_area.begin();
    ++entry_it;
    for (; entry_it != rings_by_area.end(); ++entry_it)
    {
        RingArea & entry = *entry_it;
        const Point2d & point = entry.ring->get_center();
        if (ring_encloses(*(rings_by_area.front().ring), point))
        {
            entry.ring->set(Ring::PIN, entry.ring->get_head());
            outer_polygon->push_back(*(entry.ring));
        }
        else
        {
            all_in_first = false;
            bool enclosed = false;
            const Point2d & center = entry.ring->get_center();
            for (auto & outer : polygons)
            {
                Ring & outer_ring = get_out_bound(*outer);
                if (ring_encloses(outer_ring, center))
                {
                    entry.ring->set(Ring::PIN, entry.ring->get_head());
                    outer->push_back(*(entry.ring));
                    enclosed = true;
                    break;
                }
            }
            if (!enclosed)
            {
                GisPolygon * new_polygon = new GisPolygon();
                new_polygon->m_id = polygon.m_id;
                entry.ring->set(Ring::POUT, entry.ring->get_head());
                new_polygon->push_back(*(entry.ring));

                polygons.push_back(new_polygon);
            }
        }
    }

    return all_in_first;
}
//
////return true if only one polygon is in polygon
