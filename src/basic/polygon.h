/*! \file polygon.h
 
\brief Classes for representing polygons
\author Jyh-Ming Lien
\author <a href="http://masc.cs.gmu.edu">MASC group</a>, George Mason University
\copyright 2010-2012 by Jyh-Ming Lien and George Mason University
*/

#pragma once
#ifndef CUSG_POLYGON_H
#define CUSG_POLYGON_H

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
        m_pos = p;
        init();
    }
    virtual ~LinkedVertex();
    void set_next(LinkedVertex * n)
    {
        m_next = n;
        if (n != nullptr) n->m_pre = this;
    }
    void set_pre(LinkedVertex * n)
    {
        m_pre = n;
        if (n != nullptr) n->m_next = this;
    }
    void compute_extra_info();

    /// reverse the order
    void reverse();

    ///////////////////////////////////////////////////////////////////////////
    void set_pos(const Point2d & p) { m_pos = p; }
    virtual const Point2d & get_pos() const { return m_pos; }

    void translate(const Vector2d & v) { m_pos = m_pos + v; }

    virtual LinkedVertex * get_next() const { return m_next; }
    virtual LinkedVertex * get_pre() const { return m_pre; }

    const Vector2d & get_normal() const { return m_normal; }
    bool is_reflex() const { return m_reflex; }

    /// get extra information
    uint get_vid() const { return m_vid; }
    void set_vid(uint id) { m_vid = id; }

private:
    void init()
    {
        m_next = m_pre = nullptr;
        m_reflex = false;
        m_vid = UINT_MAX;
    }

    //basic info
    Point2d m_pos;         //!< position
    LinkedVertex * m_next; //!< next vertex in the polygon
    LinkedVertex * m_pre;  //!< previous vertex in the polygon
    Vector2d m_normal;     //!< normal, the segment normal from this v to the next.
    bool m_reflex;
    uint m_vid;

public:
    void rotate(double cx, double cy, double angle) { m_pos.rotate(cx, cy, angle); }
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
        m_head = m_tail = nullptr;
        m_type = UNKNOWN;
        m_radius = -1;
        m_area = -FLT_MAX;
    }
    Ring(POLYTYPE t)
    {
        m_head = m_tail = nullptr;
        m_type = t;
        m_radius = -1;
        m_area = -FLT_MAX;
    }

    void destroy();

    ///////////////////////////////////////////////////////////////////////////
    /// create Ring
    void begin_poly();
    void add_vertex(double x, double y, bool remove_duplicate = false);
    void end_poly(bool remove_duplicate = false);

    void reverse(); //!< reverse vertex order

    ///////////////////////////////////////////////////////////////////////////
    // Access
    LinkedVertex * get_head() const { return m_head; }
    POLYTYPE get_type() const { return m_type; }
    void set(POLYTYPE t, LinkedVertex * h)
    {
        m_type = t;
        m_head = h;
        if (h != nullptr)
        {
            m_tail = h->get_pre();
        }
        else
        {
            m_tail = nullptr;
        }
    }

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
    // additional functions
    const Point2d & get_center();

    //
    const Point2d get_center() const;

    //
    float get_radius() const;

    /// area
    float get_area();

    ///////////////////////////////////////////////////////////////////////////
    // Operator
    /// check if given poly line is the same as this
    bool operator==(const Ring & other) { return other.m_head == m_head; }
    friend istream & operator>>(istream &, Ring &);

protected:
    ///////////////////////////////////////////////////////////////////////////
    void do_init(); /*return # of vertices in this poly*/

    /// build elements in vector<LinkedVertex*> all
    void build_all();

private:
    LinkedVertex * m_head; //!< the head of vertex list
    LinkedVertex * m_tail; //!< end of the vertex list

    vector<LinkedVertex *> m_all; //!< all vertices

    //additional info
    Point2d m_center;
    float m_radius;
    float m_area;

    /// In, out or unknown.
    POLYTYPE m_type;
};

/// a list of \c Ring
class RingList : public list<Ring>
{
    friend istream & operator>>(istream &, RingList &);

public:
    RingList()
    {
        m_box[0] = m_box[1] = m_box[2] = m_box[3] = 0;
        m_build_box_and_center_called = false;
    }

    //access
    void build_box_and_center();
    double * get_bounding_box()
    {
        build_box_and_center();
        return m_box;
    }
    const double * get_bounding_box() const
    {
        assert(m_build_box_and_center_called);
        return m_box;
    }
    const Point2d & get_center()
    {
        build_box_and_center();
        return m_center;
    }
    const Point2d & get_center() const
    {
        assert(m_build_box_and_center_called);
        return m_center;
    }

protected:
    Point2d m_center;
    double m_box[4]; //min_x, max_x, min_y, max_y

private:
    bool m_build_box_and_center_called;
};

/*! \brief a restricted kind of \c RingList
     
        this defines a simple polygon so that
        the first element must be a \c POUT \c Ring and
        the rest ring lines are a list of holes
    */
class Polygon : public RingList
{
public:
    Polygon() { m_area = 0; }

    /// returns the number of vertices in the polygon
    int get_size()
    {
        if (m_all.empty()) build_all();
        return m_all.size();
    }

    /// access the vertices of the polygon as an array
    LinkedVertex * operator[](unsigned int id)
    {
        if (m_all.empty()) build_all();
        return m_all[id];
    }

    double get_area();

    /// destroy
    void destroy();

    /// get number of vertices
    uint get_size() const;

    ///rotate this polygon
    void rotate(double cx, double cy, double angle)
    {
        list<Ring>::iterator pit = this->begin();
        for (; pit != this->end(); ++pit)
        {
            Ring & ring = *pit;
            LinkedVertex * phead = ring.get_head();
            LinkedVertex * pcur = phead;
            do
            {
                pcur->rotate(cx, cy, angle);
                pcur = pcur->get_next();
            } while (pcur != phead);
        }
    }

protected:
    void build_all();

    vector<LinkedVertex *> m_all; //!< all vertices

    float m_area;
};

class GisPolygon : public Polygon
{
public:
    GisPolygon() : Polygon() { m_id = 0; }
    int m_id;
};

} //namespace cusg

#endif //CUSG_POLYGON_H
