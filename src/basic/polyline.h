
#pragma once
#ifndef _POLYLINE_H_
#define _POLYLINE_H_

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
    Polyline() { head = tail = nullptr; }

    void destroy();

    ///////////////////////////////////////////////////////////////////////////
    // create Ring
    void beginPoly();
    void addVertex(double x, double y, bool remove_duplicate = false);
    void endPoly();

    ///////////////////////////////////////////////////////////////////////////
    // Access
    LinkedVertex * getHead() const { return head; }
    LinkedVertex * getTail() const { return tail; }

    int getSize()
    {
        if (all.empty()) build_all();
        return all.size();
    }

    LinkedVertex * operator[](unsigned int id)
    {
        if (all.empty()) build_all();
        return all[id];
    }

    ///////////////////////////////////////////////////////////////////////////
    // Operator
    //check if give poly line is the same as this
    bool operator==(const Polyline & other) { return other.head == head; }
    friend istream & operator>>(istream &, Polyline &);

    //rotate function
    void rotate(double cx, double cy, double angle);

protected:
    ///////////////////////////////////////////////////////////////////////////
    void doInit(); /*return # of vertices in this poly*/

    //build elements in vector<LinkedVertex*> all
    void build_all();

    LinkedVertex * head; //the head of vertex list
    LinkedVertex * tail; //end of the vertex list

    vector<LinkedVertex *> all; //all vertices
};

class GIS_polyline : public Polyline
{
public:
    GIS_polyline() : Polyline() { m_id = 0; }
    GIS_polyline(const GIS_polyline & other) : Polyline((Polyline &)other) { m_id = other.m_id; }

    //return the i-th line segment from this polyline
    GIS_LineSeg2d getLineSegment(int i);

    //compute distance for geodesic_dist
    void computeGeodesicDist();

    //store attributes from GIS
    int m_id;
    vector<double> geodesic_dist; ///geodesic distance of each vertex
};

} //end of namespace cusg

#endif //_POLYLINE_H_
