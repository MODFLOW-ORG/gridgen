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
typedef unsigned int uint;

namespace cusg
{
/// a triangle
struct triangle
{
    uint v[3]; //!< id to the vertices
};

/// Vertex of polygon
class ply_vertex
{
public:
    ///////////////////////////////////////////////////////////////////////////
    ply_vertex() { init(); }
    ply_vertex(const Point2d & p)
    {
        pos = p;
        init();
    }
    virtual ~ply_vertex();
    void setNext(ply_vertex * n)
    {
        next = n;
        if (n != NULL) n->pre = this;
    }
    void setPre(ply_vertex * n)
    {
        pre = n;
        if (n != NULL) n->next = this;
    }
    void computeExtraInfo();

    /// reverse the order
    void reverse();

    ///////////////////////////////////////////////////////////////////////////
    void setPos(const Point2d & p) { pos = p; }
    virtual const Point2d & getPos() const { return pos; }

    void translate(const Vector2d & v) { pos = pos + v; }

    virtual ply_vertex * getNext() const { return next; }
    virtual ply_vertex * getPre() const { return pre; }

    const Vector2d & getNormal() const { return normal; }
    bool isReflex() const { return reflex; }

    /// get extra information
    uint getVID() const { return vid; }
    void setVID(uint id) { vid = id; }

private:
    void init()
    {
        next = pre = NULL;
        reflex = false;
        vid = UINT_MAX;
    }

    //basic info
    Point2d pos;       //!< position
    ply_vertex * next; //!< next vertex in the polygon
    ply_vertex * pre;  //!< previous vertex in the polygon
    Vector2d normal;   //!< normal, the segment normal from this v to the next.
    bool reflex;
    uint vid;

public:
    void rotate(double cx, double cy, double angle) { pos.rotate(cx, cy, angle); }
};

/// Polygon chain
class c_ply
{
public:
    enum POLYTYPE
    {
        UNKNOWN,
        PIN,
        POUT
    };

    ///////////////////////////////////////////////////////////////////////////
    //c_ply(const c_ply& cp);//add

    ///////////////////////////////////////////////////////////////////////////
    c_ply()
    {
        head = tail = NULL;
        type = UNKNOWN;
        radius = -1;
        area = -FLT_MAX;
    }
    c_ply(POLYTYPE t)
    {
        head = tail = NULL;
        type = t;
        radius = -1;
        area = -FLT_MAX;
    }

    void destroy();

    ///////////////////////////////////////////////////////////////////////////
    /// create c_ply
    void beginPoly();
    void addVertex(double x, double y, bool remove_duplicate = false);
    void endPoly(bool remove_duplicate = false);

    void reverse(); //!< reverse vertex order

    ///////////////////////////////////////////////////////////////////////////
    // Access
    ply_vertex * getHead() const { return head; }
    POLYTYPE getType() const { return type; }
    void set(POLYTYPE t, ply_vertex * h)
    {
        type = t;
        head = h;
        if (h != NULL)
        {
            tail = h->getPre();
        }
        else
        {
            tail = NULL;
        }
    }

    int getSize()
    {
        if (all.empty()) build_all();
        return all.size();
    }

    ply_vertex * operator[](unsigned int id)
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
    bool operator==(const c_ply & other) { return other.head == head; }
    friend istream & operator>>(istream &, c_ply &);

#if TRIANGULATION

    ///////////////////////////////////////////////////////////////////////////
    // functionality provided if triangulation is available.

    /// find a point that is enclosed by this polychain
    Point2d findEnclosedPt();

    /// triangulate the polygon
    void triangulate(vector<triangle> & tris);

    /*! \brief check if a point is enclosed
            \remarks the behavior is unknown if pt is 
                     on the boundary of the polygon
        */
    bool enclosed(const Point2d & pt);

#endif //TRIANGULATION

protected:
    ///////////////////////////////////////////////////////////////////////////
    void doInit(); /*return # of vertice in this poly*/

    /// build elements in vector<ply_vertex*> all
    void build_all();

private:
    ply_vertex * head; //!< the head of vertex list
    ply_vertex * tail; //!< end of the vertex list

    vector<ply_vertex *> all; //!< all vertices

    //additional info
    Point2d center;
    float radius;
    float area;

    /// In, out or unknown.
    POLYTYPE type;

    //triangulation
    vector<triangle> triangulation; //!< cached triangulation, calculated by triangulate
};

/// a list of \c c_ply
class c_plylist : public list<c_ply>
{
    friend istream & operator>>(istream &, c_plylist &);

public:
    c_plylist()
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

/*! \brief a restricted kind of \c c_plylist
     
        this defines a simple polygon so that
        the first element must be a \c POUT \c c_ply and
        the rest ply lines are a list of holes
    */
class c_polygon : public c_plylist
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
    ply_vertex * operator[](unsigned int id)
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
        list<c_ply>::iterator pit = this->begin();
        for (; pit != this->end(); ++pit)
        {
            c_ply & ply = *pit;
            ply_vertex * phead = ply.getHead();
            ply_vertex * pcur = phead;
            do
            {
                pcur->rotate(cx, cy, angle);
                pcur = pcur->getNext();
            } while (pcur != phead);
        }
    }

#if TRIANGULATION

    ///////////////////////////////////////////////////////////////////////////
    // functionality provided if triangulation is available.

    /// triangulate the polygon
    void triangulate(vector<triangle> & tris);

    /*! \brief check if a point is enclosed
            \remarks the behavior is unknown if pt 
                     is on the boundary of the polygon
        */
    bool enclosed(const Point2d & pt);

    /// find a point inside the polygon
    Point2d findEnclosedPt();

#endif //TRIANGULATION

protected:
    void build_all();

    vector<ply_vertex *> all; //!< all vertices

    //triangulation
    vector<triangle> triangulation; //!< cached triangulation, calculated by triangulate

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
