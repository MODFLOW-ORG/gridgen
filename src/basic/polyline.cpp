

#include "polyline.h"

#include <vector>

namespace cusg
{

// clean up the space allocated
void c_plyline::destroy()
{
    if (head == NULL) return;
    ply_vertex * ptr = head;
    do
    {
        ply_vertex * n = ptr->getNext();
        delete ptr;
        ptr = n;
    } while (ptr != NULL);
    head = tail = NULL;

    all.clear();
}

// Create a empty polygon
void c_plyline::beginPoly()
{
    head = tail = NULL;
    all.clear();
}

// Add a vertex to the polygonal chian
void c_plyline::addVertex(double x, double y, bool remove_duplicate)
{
    Point2d pt(x, y);

    if (tail != NULL)
    {
        if (tail->getPos() == pt && remove_duplicate) return; //don't add
    }

    ply_vertex * v = new ply_vertex(pt);
    if (tail != NULL)
    {
        tail->setNext(v);
    }
    tail = v;
    if (head == NULL) head = tail;
    v->setVID(all.size()); //id of the vertex in this ply
    all.push_back(v);
}

// finish building the polygon
void c_plyline::endPoly()
{
    doInit();
}

// initialize property of the this polychain
// Compute normals and find reflective vertices
void c_plyline::doInit()
{
    //compute normals
    ply_vertex * ptr = head;
    do
    {
        ptr->computeExtraInfo();
        ptr = ptr->getNext();
    } while (ptr != NULL);
}

///////////////////////////////////////////////////////////////////////////

void c_plyline::build_all()
{
    uint vid = 0;
    all.clear();
    ply_vertex * ptr = head;
    do
    {
        ptr->setVID(vid++);
        all.push_back(ptr);
        ptr = ptr->getNext();
    } while (ptr != NULL); //end while
}

//return the i-th line segment from this polyline
GIS_LineSeg2d GIS_plyline::getLineSegment(int i)
{
    int vsize = getSize();
    if (i >= vsize)
    {
        cerr << "! Error: GIS_plyline::getLineSegment index out of bound i=" << i
             << " polyline seg size=" << vsize - 1 << endl;
        assert(false);
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

void GIS_plyline::computeGeodesicDist()
{
    int vsize = getSize();
    geodesic_dist = vector<double>(vsize, 0);

    ply_vertex * ptr = head->getNext();
    double dist = 0;
    for (int i = 1; i < vsize; i++)
    {
        ply_vertex * pre = ptr->getPre();
        dist += (ptr->getPos() - pre->getPos()).norm();
        geodesic_dist[i] = dist;
        ptr = ptr->getNext();
    }
}
//
// Compute the center and the box of a list of plys
//

//
//
// IO operators
//
istream & operator>>(istream & is, c_plyline & poly)
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

void c_plyline::rotate(double cx, double cy, double angle)
{
    ply_vertex * head = this->getHead();
    ply_vertex * p = head;
    do
    {
        p->rotate(cx, cy, angle);
        p = p->getNext();
    } while (p != NULL);
}

} //end of namespace cusg
