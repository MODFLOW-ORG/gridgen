/*! \file polygon.h
 
\brief Classes for representing polygons
\author Jyh-Ming Lien
\author <a href="http://masc.cs.gmu.edu">MASC group</a>, George Mason University
\copyright 2010-2012 by Jyh-Ming Lien and George Mason University
*/

#pragma once
#ifndef _CUSG_POLYGON_H_
#define _CUSG_POLYGON_H_

#if defined(_WIN32) || defined(_WIN64)
#pragma warning(disable : 4786)
#endif

#include <Basic.h>
#include <Point.h>
#include <Vector.h>

#include <list>
#include <cassert>
#include <vector>
using namespace std;

#include <limits.h>
#include <float.h>
using uint = unsigned int;

namespace cusg
{
/// Vertex of polygon
class LinkedVertex
{
public:
    ///////////////////////////////////////////////////////////////////////////
    LinkedVertex() { init(); }
    LinkedVertex(const Point2d & p)
    {
        pos = p;
        init();
    }
    virtual ~LinkedVertex();
    void setNext(LinkedVertex * n)
    {
        next = n;
        if (n != nullptr) n->pre = this;
    }
    void setPre(LinkedVertex * n)
    {
        pre = n;
        if (n != nullptr) n->next = this;
    }
    void computeExtraInfo();

    /// reverse the order
    void reverse();

    ///////////////////////////////////////////////////////////////////////////
    void setPos(const Point2d & p) { pos = p; }
    virtual const Point2d & getPos() const { return pos; }

    void translate(const Vector2d & v) { pos = pos + v; }

    virtual LinkedVertex * getNext() const { return next; }
    virtual LinkedVertex * getPre() const { return pre; }

    const Vector2d & getNormal() const { return normal; }
    bool isReflex() const { return reflex; }

    /// get extra information
    uint getVID() const { return vid; }
    void setVID(uint id) { vid = id; }

private:
    void init()
    {
        next = pre = nullptr;
        reflex = false;
        vid = UINT_MAX;
    }

    //basic info
    Point2d pos;         //!< position
    LinkedVertex * next; //!< next vertex in the polygon
    LinkedVertex * pre;  //!< previous vertex in the polygon
    Vector2d normal;     //!< normal, the segment normal from this v to the next.
    bool reflex;
    uint vid;

public:
    void rotate(double cx, double cy, double angle) { pos.rotate(cx, cy, angle); }
};

/// Polygon chain
class Ring
{
public:
    enum POLYTYPE
    {
        UNKNOWN,
        PIN,
        POUT
    };

    ///////////////////////////////////////////////////////////////////////////
    //Ring(const Ring& cp);//add

    ///////////////////////////////////////////////////////////////////////////
    Ring()
    {
        head = tail = nullptr;
        type = UNKNOWN;
        radius = -1;
        area = -FLT_MAX;
    }
    Ring(POLYTYPE t)
    {
        head = tail = nullptr;
        type = t;
        radius = -1;
        area = -FLT_MAX;
    }

    void destroy();

    ///////////////////////////////////////////////////////////////////////////
    /// create Ring
    void beginPoly();
    void addVertex(double x, double y, bool remove_duplicate = false);
    void endPoly(bool remove_duplicate = false);

    void reverse(); //!< reverse vertex order

    ///////////////////////////////////////////////////////////////////////////
    // Access
    LinkedVertex * getHead() const { return head; }
    POLYTYPE getType() const { return type; }
    void set(POLYTYPE t, LinkedVertex * h)
    {
        type = t;
        head = h;
        if (h != nullptr)
        {
            tail = h->getPre();
        }
        else
        {
            tail = nullptr;
        }
    }

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
    // additional functions
    const Point2d & getCenter();

    //
    const Point2d getCenter() const;

    //
    float getRadius() const;

    /// area
    float getArea();

    ///////////////////////////////////////////////////////////////////////////
    // Operator
    /// check if given poly line is the same as this
    bool operator==(const Ring & other) { return other.head == head; }
    friend istream & operator>>(istream &, Ring &);

protected:
    ///////////////////////////////////////////////////////////////////////////
    void doInit(); /*return # of vertices in this poly*/

    /// build elements in vector<LinkedVertex*> all
    void build_all();

private:
    LinkedVertex * head; //!< the head of vertex list
    LinkedVertex * tail; //!< end of the vertex list

    vector<LinkedVertex *> all; //!< all vertices

    //additional info
    Point2d center;
    float radius;
    float area;

    /// In, out or unknown.
    POLYTYPE type;
};

/// a list of \c Ring
class RingList : public list<Ring>
{
    friend istream & operator>>(istream &, RingList &);

public:
    RingList()
    {
        box[0] = box[1] = box[2] = box[3] = 0;
        is_buildboxandcenter_called = false;
    }

    //access
    void buildBoxAndCenter();
    double * getBBox()
    {
        buildBoxAndCenter();
        return box;
    }
    const double * getBBox() const
    {
        assert(is_buildboxandcenter_called);
        return box;
    }
    const Point2d & getCenter()
    {
        buildBoxAndCenter();
        return center;
    }
    const Point2d & getCenter() const
    {
        assert(is_buildboxandcenter_called);
        return center;
    }

protected:
    Point2d center;
    double box[4]; //min_x, max_x, min_y, max_y

private:
    bool is_buildboxandcenter_called;
};

/*! \brief a restricted kind of \c RingList
     
        this defines a simple polygon so that
        the first element must be a \c POUT \c Ring and
        the rest ring lines are a list of holes
    */
class c_polygon : public RingList
{
public:
    c_polygon() { area = 0; }

    /// returns the number of vertices in the polygon
    int getSize()
    {
        if (all.empty()) build_all();
        return all.size();
    }

    /// access the vertices of the polygon as an array
    LinkedVertex * operator[](unsigned int id)
    {
        if (all.empty()) build_all();
        return all[id];
    }

    double getArea();

    /// destroy
    void destroy();

    /// get number of vertices
    uint getSize() const;

    ///rotate this polygon
    void rotate(double cx, double cy, double angle)
    {
        list<Ring>::iterator pit = this->begin();
        for (; pit != this->end(); ++pit)
        {
            Ring & ring = *pit;
            LinkedVertex * phead = ring.getHead();
            LinkedVertex * pcur = phead;
            do
            {
                pcur->rotate(cx, cy, angle);
                pcur = pcur->getNext();
            } while (pcur != phead);
        }
    }

protected:
    void build_all();

    vector<LinkedVertex *> all; //!< all vertices

    float area;
};

class GIS_polygon : public c_polygon
{
public:
    GIS_polygon() : c_polygon() { m_id = 0; }
    int m_id;
};

} //namespace cusg

#endif //_CUSG_POLYGON_H_
