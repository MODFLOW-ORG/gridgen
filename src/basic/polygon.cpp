
#include <vector>
#include "polygon.h"
#include "intersection.h"

namespace cusg
{
//    Polygon g_P;
//    Polygon g_Q;
//    Polygon& getP(){ return g_P; }
//    Polygon& getQ(){ return g_Q; }

LinkedVertex::~LinkedVertex()
{
    //doing nothing for now
}

// - compute normal
// - check if the vertex is reflex or not
void LinkedVertex::compute_extra_info()
{
    m_reflex = true;

    //compute normal direction
    if (m_next != nullptr)
    {
        Vector2d v = m_next->m_pos - m_pos;
        if (v[0] == 0)
        {
            if (v[1] > 0)
            {
                m_normal[0] = 1;
                m_normal[1] = 0;
            }
            else
            {
                m_normal[0] = -1;
                m_normal[1] = 0;
            }
        }
        else if (v[0] > 0)
        {
            m_normal[1] = -1;
            m_normal[0] = (v[1] / v[0]);
        }
        else
        { //v[0]<0
            m_normal[1] = 1;
            m_normal[0] = -(v[1] / v[0]);
        }

        m_normal = m_normal.normalize();

        if (m_pre != nullptr)
        {
            //compute if left or right turn
            Vector2d u = m_pos - m_pre->m_pos;
            float z = u[0] * v[1] - u[1] * v[0];
            if (z > 0) m_reflex = false;
        }
    } //end if(next!=NULL)
}

void LinkedVertex::reverse()
{
    swap(m_next, m_pre);
    compute_extra_info();
}

///////////////////////////////////////////////////////////////////////////////

// clean up the space allocated
void Ring::destroy()
{
    if (m_head == nullptr) return;
    LinkedVertex * ptr = m_head;
    do
    {
        LinkedVertex * n = ptr->get_next();
        delete ptr;
        ptr = n;
    } while (ptr != m_head);
    m_head = m_tail = nullptr;

    m_all.clear();
}

// Create a empty polygon
void Ring::begin_poly()
{
    m_head = m_tail = nullptr;
    m_all.clear();
}

// Add a vertex to the polygonal chain
void Ring::add_vertex(double x, double y, bool remove_duplicate)
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
void Ring::end_poly(bool remove_duplicate)
{
    if (m_head != nullptr && m_tail != nullptr)
    {
        if (remove_duplicate)
        {
            if (m_head->get_pos() == m_tail->get_pos())
            { //remove tail..
                delete m_tail;
                m_all.pop_back();
                m_tail = m_all.back();
            }
        } //
    }

    m_tail->set_next(m_head);
    do_init();
}

// initialize property of the this polychain
// Compute normals and find reflective vertices
void Ring::do_init()
{
    //compute area
    get_area();
    if (this->m_area < 0 && m_type == POUT)
    {
        reverse();
    }
    else if (this->m_area > 0 && m_type == PIN)
    {
        reverse();
    }

    //compute normals
    LinkedVertex * ptr = m_head;
    do
    {
        ptr->compute_extra_info();
        ptr = ptr->get_next();
    } while (ptr != m_head);
}

const Point2d & Ring::get_center()
{
    if (m_radius <= 0)
    {
        this->m_center.set(0, 0);
        LinkedVertex * ptr = m_head;
        const Point2d & first = ptr->get_pos();
        uint size = 0;
        do
        {
            size++;
            Vector2d v = ptr->get_pos() - first;
            this->m_center[0] += v[0];
            this->m_center[1] += v[1];
            ptr = ptr->get_next();
        } while (ptr != m_head); //end while
        this->m_center[0] = (m_center[0] / size) + first[0];
        this->m_center[1] = (m_center[1] / size) + first[1];

        m_radius = 0;
    }

    return this->m_center;
}

const Point2d Ring::get_center() const
{
    if (m_radius < 0)
    {
        Point2d center;
        LinkedVertex * ptr = m_head;
        const Point2d & first = ptr->get_pos();
        uint size = 0;
        do
        {
            size++;
            Vector2d v = ptr->get_pos() - first;
            center[0] += v[0];
            center[1] += v[1];
            ptr = ptr->get_next();
        } while (ptr != m_head); //end while
        center[0] = (center[0] / size) + first[0];
        center[1] = (center[1] / size) + first[1];
        return center;
    }
    else
        return this->m_center;
}

//reverse the order of the vertices
void Ring::reverse()
{
    LinkedVertex * ptr = m_head;
    do
    {
        ptr->reverse();
        ptr = ptr->get_next();
    } while (ptr != m_head); //end while

    this->m_area = -this->m_area;
    m_all.clear();
}

float Ring::get_radius() const
{
    if (m_radius <= 0)
    {
        double radius = 0;
        get_center();
        LinkedVertex * ptr = m_head;
        do
        {
            float d = (m_center - ptr->get_pos()).normsqr();
            if (d > radius) radius = d;
            ptr = ptr->get_next();
        } while (ptr != m_head); //end while
        radius = sqrt(radius);

        return radius;
    }
    else
        return m_radius;
}

float Ring::get_area()
{
    if (m_area == -FLT_MAX)
    {
        m_area = 0;
        get_center();
        LinkedVertex * ptr = m_head;
        do
        {
            LinkedVertex * next = ptr->get_next();
            const Point2d & p1 = ptr->get_pos();
            const Point2d & p2 = next->get_pos();
            m_area += area(p1.get(), p2.get(), m_center.get());
            ptr = next;
        } while (ptr != m_head); //end while
    }

    return m_area;
}

void Ring::build_all()
{
    uint vid = 0;
    m_all.clear();
    LinkedVertex * ptr = m_head;
    do
    {
        ptr->set_vid(vid++);
        m_all.push_back(ptr);
        ptr = ptr->get_next();
    } while (ptr != m_head); //end while
}

//
// Compute the center and the box of a list of rings
//

void RingList::build_box_and_center()
{
    if (m_build_box_and_center_called == false)
    {
        m_box[0] = m_box[2] = FLT_MAX;
        m_box[1] = m_box[3] = -FLT_MAX;
        for (auto & i : *this)
        {

            LinkedVertex * ptr = i.get_head();
            do
            {
                const Point2d & p = ptr->get_pos();
                if (p[0] < m_box[0]) m_box[0] = p[0];
                if (p[0] > m_box[1]) m_box[1] = p[0];
                if (p[1] < m_box[2]) m_box[2] = p[1];
                if (p[1] > m_box[3]) m_box[3] = p[1];
                ptr = ptr->get_next();
            } while (ptr != i.get_head()); //end while
        }

        m_center[0] = (m_box[0] + m_box[1]) / 2;
        m_center[1] = (m_box[2] + m_box[3]) / 2;

        m_build_box_and_center_called = true;
    }
}

void Polygon::destroy()
{
    for (auto & i : *this)
    {
        i.destroy();
    }
    clear(); //remove all ring from this list
    m_all.clear();
}

void Polygon::build_all()
{
    uint vid = 0;
    for (auto & i : *this)
    {
        uint vsize = i.get_size();
        for (uint j = 0; j < vsize; j++)
        {
            i[j]->set_vid(vid++); //id of vertex in the polygons
            m_all.push_back(i[j]);
        }
    } //end for i
}

double Polygon::get_area()
{
    if (m_area == 0)
    {
        for (auto & i : *this)
        {
            m_area += i.get_area();
        } //end for i
    }

    return m_area;
}

istream & operator>>(istream & is, Ring & poly)
{
    int vsize;
    string str_type;
    is >> vsize >> str_type;

    if (str_type.find("out") != string::npos)
        poly.m_type = Ring::POUT;
    else
        poly.m_type = Ring::PIN;

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

istream & operator>>(istream & is, RingList & p)
{
    //remove header comments
    do
    {
        char line[1024];
        char c = is.peek();
        if (isspace(c))
            is.get(c); //eat it
        else if (c == '#')
        {
            is.getline(line, 1024);
        }
        else
            break;
    } while (true);

    //start reading
    uint size;
    is >> size;
    uint vid = 0;
    for (uint i = 0; i < size; i++)
    {
        Ring poly(Ring::UNKNOWN);
        is >> poly;
        p.push_back(poly);
        uint vsize = poly.get_size();
        for (uint j = 0; j < vsize; j++)
        {
            poly[j]->set_vid(vid++); //id of vertex in the polygons
        }
    }
    return is;
}

} // end namespace mksum2d
