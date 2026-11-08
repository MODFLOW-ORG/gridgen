
#include <vector>
#include "polygon.h"
#include "intersection.h"

namespace cusg
{
//    c_polygon g_P;
//    c_polygon g_Q;
//    c_polygon& getP(){ return g_P; }
//    c_polygon& getQ(){ return g_Q; }

LinkedVertex::~LinkedVertex()
{
    //doing nothing for now
}

// - compute normal
// - check if the vertex is reflex or not
void LinkedVertex::computeExtraInfo()
{
    reflex = true;

    //compute normal direction
    if (next != nullptr)
    {
        Vector2d v = next->pos - pos;
        if (v[0] == 0)
        {
            if (v[1] > 0)
            {
                normal[0] = 1;
                normal[1] = 0;
            }
            else
            {
                normal[0] = -1;
                normal[1] = 0;
            }
        }
        else if (v[0] > 0)
        {
            normal[1] = -1;
            normal[0] = (v[1] / v[0]);
        }
        else
        { //v[0]<0
            normal[1] = 1;
            normal[0] = -(v[1] / v[0]);
        }

        normal = normal.normalize();

        if (pre != nullptr)
        {
            //compute if left or right turn
            Vector2d u = pos - pre->pos;
            float z = u[0] * v[1] - u[1] * v[0];
            if (z > 0) reflex = false;
        }
    } //end if(next!=NULL)
}

void LinkedVertex::reverse()
{
    swap(next, pre);
    computeExtraInfo();
}

///////////////////////////////////////////////////////////////////////////////

// clean up the space allocated
void Ring::destroy()
{
    if (head == nullptr) return;
    LinkedVertex * ptr = head;
    do
    {
        LinkedVertex * n = ptr->getNext();
        delete ptr;
        ptr = n;
    } while (ptr != head);
    head = tail = nullptr;

    all.clear();
}

// Create a empty polygon
void Ring::beginPoly()
{
    head = tail = nullptr;
    all.clear();
}

// Add a vertex to the polygonal chain
void Ring::addVertex(double x, double y, bool remove_duplicate)
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
void Ring::endPoly(bool remove_duplicate)
{
    if (head != nullptr && tail != nullptr)
    {
        if (remove_duplicate)
        {
            if (head->getPos() == tail->getPos())
            { //remove tail..
                delete tail;
                all.pop_back();
                tail = all.back();
            }
        } //
    }

    tail->setNext(head);
    doInit();
}

// initialize property of the this polychain
// Compute normals and find reflective vertices
void Ring::doInit()
{
    //compute area
    getArea();
    if (this->area < 0 && type == POUT)
    {
        reverse();
    }
    else if (this->area > 0 && type == PIN)
    {
        reverse();
    }

    //compute normals
    LinkedVertex * ptr = head;
    do
    {
        ptr->computeExtraInfo();
        ptr = ptr->getNext();
    } while (ptr != head);
}

const Point2d & Ring::getCenter()
{
    if (radius <= 0)
    {
        this->center.set(0, 0);
        LinkedVertex * ptr = head;
        const Point2d & first = ptr->getPos();
        uint size = 0;
        do
        {
            size++;
            Vector2d v = ptr->getPos() - first;
            this->center[0] += v[0];
            this->center[1] += v[1];
            ptr = ptr->getNext();
        } while (ptr != head); //end while
        this->center[0] = (center[0] / size) + first[0];
        this->center[1] = (center[1] / size) + first[1];

        radius = 0;
    }

    return this->center;
}

const Point2d Ring::getCenter() const
{
    if (radius < 0)
    {
        Point2d center;
        LinkedVertex * ptr = head;
        const Point2d & first = ptr->getPos();
        uint size = 0;
        do
        {
            size++;
            Vector2d v = ptr->getPos() - first;
            center[0] += v[0];
            center[1] += v[1];
            ptr = ptr->getNext();
        } while (ptr != head); //end while
        center[0] = (center[0] / size) + first[0];
        center[1] = (center[1] / size) + first[1];
        return center;
    }
    else
        return this->center;
}

//reverse the order of the vertices
void Ring::reverse()
{
    LinkedVertex * ptr = head;
    do
    {
        ptr->reverse();
        ptr = ptr->getNext();
    } while (ptr != head); //end while

    this->area = -this->area;
    all.clear();
}

float Ring::getRadius() const
{
    if (radius <= 0)
    {
        double radius = 0;
        getCenter();
        LinkedVertex * ptr = head;
        do
        {
            float d = (center - ptr->getPos()).normsqr();
            if (d > radius) radius = d;
            ptr = ptr->getNext();
        } while (ptr != head); //end while
        radius = sqrt(radius);

        return radius;
    }
    else
        return radius;
}

float Ring::getArea()
{
    if (area == -FLT_MAX)
    {
        area = 0;
        getCenter();
        LinkedVertex * ptr = head;
        do
        {
            LinkedVertex * next = ptr->getNext();
            const Point2d & p1 = ptr->getPos();
            const Point2d & p2 = next->getPos();
            area += Area(p1.get(), p2.get(), center.get());
            ptr = next;
        } while (ptr != head); //end while
    }

    return area;
}

void Ring::build_all()
{
    uint vid = 0;
    all.clear();
    LinkedVertex * ptr = head;
    do
    {
        ptr->setVID(vid++);
        all.push_back(ptr);
        ptr = ptr->getNext();
    } while (ptr != head); //end while
}

//
// Compute the center and the box of a list of rings
//

void RingList::buildBoxAndCenter()
{
    if (is_buildboxandcenter_called == false)
    {
        box[0] = box[2] = FLT_MAX;
        box[1] = box[3] = -FLT_MAX;
        for (auto & i : *this)
        {

            LinkedVertex * ptr = i.getHead();
            do
            {
                const Point2d & p = ptr->getPos();
                if (p[0] < box[0]) box[0] = p[0];
                if (p[0] > box[1]) box[1] = p[0];
                if (p[1] < box[2]) box[2] = p[1];
                if (p[1] > box[3]) box[3] = p[1];
                ptr = ptr->getNext();
            } while (ptr != i.getHead()); //end while
        }

        center[0] = (box[0] + box[1]) / 2;
        center[1] = (box[2] + box[3]) / 2;

        is_buildboxandcenter_called = true;
    }
}

void c_polygon::destroy()
{
    for (auto & i : *this)
    {
        i.destroy();
    }
    clear(); //remove all ring from this list
    all.clear();
}

void c_polygon::build_all()
{
    uint vid = 0;
    for (auto & i : *this)
    {
        uint vsize = i.getSize();
        for (uint j = 0; j < vsize; j++)
        {
            i[j]->setVID(vid++); //id of vertex in the polygons
            all.push_back(i[j]);
        }
    } //end for i
}

double c_polygon::getArea()
{
    if (area == 0)
    {
        for (auto & i : *this)
        {
            area += i.getArea();
        } //end for i
    }

    return area;
}

istream & operator>>(istream & is, Ring & poly)
{
    int vsize;
    string str_type;
    is >> vsize >> str_type;

    if (str_type.find("out") != string::npos)
        poly.type = Ring::POUT;
    else
        poly.type = Ring::PIN;

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
        uint vsize = poly.getSize();
        for (uint j = 0; j < vsize; j++)
        {
            poly[j]->setVID(vid++); //id of vertex in the polygons
        }
    }
    return is;
}

} // end namespace mksum2d
