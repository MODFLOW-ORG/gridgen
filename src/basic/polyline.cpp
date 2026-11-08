

#include "polyline.h"
#include "GridgenError.h"

#include <vector>

namespace cusg
{

// clean up the space allocated
void Polyline::destroy()
{
    if (m_head == nullptr) return;
    LinkedVertex * ptr = m_head;
    do
    {
        LinkedVertex * n = ptr->get_next();
        delete ptr;
        ptr = n;
    } while (ptr != nullptr);
    m_head = m_tail = nullptr;

    m_all.clear();
}

// Create a empty polygon
void Polyline::begin_poly()
{
    m_head = m_tail = nullptr;
    m_all.clear();
}

// Add a vertex to the polygonal chain
void Polyline::add_vertex(double x, double y, bool remove_duplicate)
{
    Point2d pt(x, y);

    if (m_tail != nullptr)
    {
        if (m_tail->get_pos() == pt && remove_duplicate) return; //don't add
    }

    LinkedVertex * v = new LinkedVertex(pt);
    if (m_tail != nullptr)
    {
        m_tail->set_next(v);
    }
    m_tail = v;
    if (m_head == nullptr) m_head = m_tail;
    v->set_vid(m_all.size()); //id of the vertex in this ring
    m_all.push_back(v);
}

// finish building the polygon
void Polyline::end_poly()
{
    do_init();
}

// initialize property of the this polychain
// Compute normals and find reflective vertices
void Polyline::do_init()
{
    //compute normals
    LinkedVertex * ptr = m_head;
    do
    {
        ptr->compute_extra_info();
        ptr = ptr->get_next();
    } while (ptr != nullptr);
}

///////////////////////////////////////////////////////////////////////////

void Polyline::build_all()
{
    uint vid = 0;
    m_all.clear();
    LinkedVertex * ptr = m_head;
    do
    {
        ptr->set_vid(vid++);
        m_all.push_back(ptr);
        ptr = ptr->get_next();
    } while (ptr != nullptr); //end while
}

//return the i-th line segment from this polyline
GisLineSeg2d GisPolyline::get_line_segment(int i)
{
    int vsize = get_size();
    if (i >= vsize)
    {
        throw GridgenError("internal error: a line has no segment " + to_string(i + 1));
    }

    if (geodesic_dist.empty()) compute_geodesic_dist();

    //create a GisLineSeg2d from the polyline...
    GisLineSeg2d line((*this)[i]->get_pos(), (*this)[i + 1]->get_pos());
    line.m_curve_id = this->m_id;
    line.m_seg_id = i;
    line.m_start_geo_dist = this->geodesic_dist[i];
    line.m_end_geo_dist = this->geodesic_dist[i + 1];

    return line;
}

void GisPolyline::compute_geodesic_dist()
{
    int vsize = get_size();
    geodesic_dist = vector<double>(vsize, 0);

    LinkedVertex * ptr = m_head->get_next();
    double dist = 0;
    for (int i = 1; i < vsize; i++)
    {
        LinkedVertex * pre = ptr->get_pre();
        dist += (ptr->get_pos() - pre->get_pos()).norm();
        geodesic_dist[i] = dist;
        ptr = ptr->get_next();
    }
}
//
// Compute the center and the box of a list of rings
//

//
//
// IO operators
//
istream & operator>>(istream & is, Polyline & poly)
{
    int vsize;
    is >> vsize;

    poly.begin_poly();

    //read in all the vertices
    int iv;
    vector<pair<double, double>> pts;
    pts.reserve(vsize);
    for (iv = 0; iv < vsize; iv++)
    {
        double x, y;
        is >> x >> y;
        pts.push_back(pair<double, double>(x, y));
    }

    int id;
    for (iv = 0; iv < vsize; iv++)
    {
        is >> id;
        id = id - 1;
        poly.add_vertex(pts[id].first, pts[id].second);
    }

    poly.end_poly();
    return is;
}

void Polyline::rotate(double cx, double cy, double angle)
{
    LinkedVertex * head = this->get_head();
    LinkedVertex * p = head;
    do
    {
        p->rotate(cx, cy, angle);
        p = p->get_next();
    } while (p != nullptr);
}

} //end of namespace cusg
