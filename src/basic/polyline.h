
#pragma once
#ifndef POLYLINE_H
#define POLYLINE_H

#include "polygon.h"
#include "lineseg.h"

namespace cusg
{

//
// Polyline
//
class Polyline
{

public:
    ///////////////////////////////////////////////////////////////////////////
    Polyline() { m_head = m_tail = nullptr; }

    void destroy();

    ///////////////////////////////////////////////////////////////////////////
    // create Ring
    void begin_poly();
    void add_vertex(double x, double y, bool remove_duplicate = false);
    void end_poly();

    ///////////////////////////////////////////////////////////////////////////
    // Access
    LinkedVertex * get_head() const { return m_head; }
    LinkedVertex * get_tail() const { return m_tail; }

    int get_size()
    {
        if (m_all.empty()) build_all();
        return m_all.size();
    }

    LinkedVertex * operator[](unsigned int id)
    {
        if (m_all.empty()) build_all();
        return m_all[id];
    }

    ///////////////////////////////////////////////////////////////////////////
    // Operator
    //check if give poly line is the same as this
    bool operator==(const Polyline & other) { return other.m_head == m_head; }
    friend istream & operator>>(istream &, Polyline &);

    //rotate function
    void rotate(double cx, double cy, double angle);

protected:
    ///////////////////////////////////////////////////////////////////////////
    void do_init(); /*return # of vertices in this poly*/

    //build elements in vector<LinkedVertex*> all
    void build_all();

    LinkedVertex * m_head; //the head of vertex list
    LinkedVertex * m_tail; //end of the vertex list

    vector<LinkedVertex *> m_all; //all vertices
};

class GisPolyline : public Polyline
{
public:
    GisPolyline() : Polyline() { m_id = 0; }
    GisPolyline(const GisPolyline & other) : Polyline((Polyline &)other) { m_id = other.m_id; }

    //return the i-th line segment from this polyline
    GisLineSeg2d get_line_segment(int i);

    //compute distance for geodesic_dist
    void compute_geodesic_dist();

    //store attributes from GIS
    int m_id;
    vector<double> geodesic_dist; ///geodesic distance of each vertex
};

} //end of namespace cusg

#endif //POLYLINE_H
