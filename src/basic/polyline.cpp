

#include "polyline.h"
#include "GridgenError.h"

#include <vector>

namespace cusg
{

// clean up the space allocated
void Polyline::destroy()
{
    if (head == nullptr) return;
    LinkedVertex * ptr = head;
    do
    {
        LinkedVertex * n = ptr->getNext();
        delete ptr;
        ptr = n;
    } while (ptr != nullptr);
    head = tail = nullptr;

    all.clear();
}

// Create a empty polygon
void Polyline::beginPoly()
{
    head = tail = nullptr;
    all.clear();
}

// Add a vertex to the polygonal chain
void Polyline::addVertex(double x, double y, bool remove_duplicate)
{
    Point2d pt(x, y);

    if (tail != nullptr)
    {
        if (tail->getPos() == pt && remove_duplicate) return; //don't add
    }

    LinkedVertex * v = new LinkedVertex(pt);
    if (tail != nullptr)
    {
        tail->setNext(v);
    }
    tail = v;
    if (head == nullptr) head = tail;
    v->setVID(all.size()); //id of the vertex in this ring
    all.push_back(v);
}

// finish building the polygon
void Polyline::endPoly()
{
    doInit();
}

// initialize property of the this polychain
// Compute normals and find reflective vertices
void Polyline::doInit()
{
    //compute normals
    LinkedVertex * ptr = head;
    do
    {
        ptr->computeExtraInfo();
        ptr = ptr->getNext();
    } while (ptr != nullptr);
}

///////////////////////////////////////////////////////////////////////////

void Polyline::build_all()
{
    uint vid = 0;
    all.clear();
    LinkedVertex * ptr = head;
    do
    {
        ptr->setVID(vid++);
        all.push_back(ptr);
        ptr = ptr->getNext();
    } while (ptr != nullptr); //end while
}

//return the i-th line segment from this polyline
GIS_LineSeg2d GIS_polyline::getLineSegment(int i)
{
    int vsize = getSize();
    if (i >= vsize)
    {
        throw GridgenError("internal error: a line has no segment " + to_string(i + 1));
    }

    if (geodesic_dist.empty()) computeGeodesicDist();

    //create a GIS_LineSeg2d from the polyline...
    GIS_LineSeg2d line((*this)[i]->getPos(), (*this)[i + 1]->getPos());
    line.m_curve_id = this->m_id;
    line.m_seg_id = i;
    line.m_start_geo_dist = this->geodesic_dist[i];
    line.m_end_geo_dist = this->geodesic_dist[i + 1];

    return line;
}

void GIS_polyline::computeGeodesicDist()
{
    int vsize = getSize();
    geodesic_dist = vector<double>(vsize, 0);

    LinkedVertex * ptr = head->getNext();
    double dist = 0;
    for (int i = 1; i < vsize; i++)
    {
        LinkedVertex * pre = ptr->getPre();
        dist += (ptr->getPos() - pre->getPos()).norm();
        geodesic_dist[i] = dist;
        ptr = ptr->getNext();
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

    poly.beginPoly();

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
        poly.addVertex(pts[id].first, pts[id].second);
    }

    poly.endPoly();
    return is;
}

void Polyline::rotate(double cx, double cy, double angle)
{
    LinkedVertex * head = this->getHead();
    LinkedVertex * p = head;
    do
    {
        p->rotate(cx, cy, angle);
        p = p->getNext();
    } while (p != nullptr);
}

} //end of namespace cusg
