
#include <vector>
#include "polygon.h"
#include "intersection.h"

namespace cusg
{
//    c_polygon g_P;
//    c_polygon g_Q;
//    c_polygon& getP(){ return g_P; }
//    c_polygon& getQ(){ return g_Q; }

//
//ply_vertex * ply_vertex_extra::getDihedralPre(){
//    assert(concavity_bpc);
//    ply_vertex * v=concavity_bpc->getConcavity();
//    for(short i=0;i<4;i++) v=v->getPre();
//    return v;
//}
//ply_vertex * ply_vertex_extra::getDihedralNext(){
//    assert(concavity_bpc);
//    ply_vertex * v=concavity_bpc->getConcavity();
//    for(short i=0;i<4;i++) v=v->getNext();
//    return v;
//}

ply_vertex::~ply_vertex()
{
    //doing nothing for now
}

// - compute normal
// - check if the vertex is reflex or not
void ply_vertex::computeExtraInfo()
{
    reflex = true;

    //compute normal direction
    if (next != NULL)
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

        if (pre != NULL)
        {
            //compute if left or right turn
            Vector2d u = pos - pre->pos;
            float z = u[0] * v[1] - u[1] * v[0];
            if (z > 0) reflex = false;
        }
    } //end if(next!=NULL)
}

void ply_vertex::reverse()
{
    swap(next, pre);
    //normal=-normal;
    computeExtraInfo();
}

/*
    ply_vertex::~ply_vertex()
    {
        //doing nothing for now
    }

    // - compute normal
    // - check if the vertex is reflex or not
    void ply_vertex::computeExtraInfo()
    {
        reflex=true;

         //compute normal direction
         if(next!=NULL){
             Vector2d v=next->pos-pos;
             if( v[0]==0 ){
                 if(v[1]>0){ normal[0]=1; normal[1]=0; }
                 else{ normal[0]=-1; normal[1]=0; }
             }
             else if( v[0]>0 ){
                 normal[1]=-1;
                 normal[0]=(v[1]/v[0]);
             }
             else{//v[0]<0
                 normal[1]=1;
                 normal[0]=-(v[1]/v[0]);
             }

             normal=normal.normalize();

             if(pre!=NULL){
                 //compute if left or right turn
                 Vector2d u=pos-pre->pos;
                 float z=u[0]*v[1]-u[1]*v[0];
                 if(z>0) reflex=false;
             }
         }//end if(next!=NULL)
    }

    void ply_vertex::negate()
    {
        normal=-normal;
        pos[0]=-pos[0];
        pos[1]=-pos[1];
    }

    void ply_vertex::reverse()
    {
        swap(next,pre);
        //normal=-normal;
        computeExtraInfo();
    }

    void ply_vertex::copy(ply_vertex * other)
    {
        pos=other->pos;
        normal=other->normal;
        reflex=other->reflex;
        vid=other->vid;
    }

    void ply_vertex::rotate(double r)
    {
        double cos_r=cos(r);
        double sin_r=sin(r);
        //rotate pos
        double x=pos[0]*cos_r-pos[1]*sin_r;
        double y=pos[0]*sin_r+pos[1]*cos_r;
        pos.set(x,y);

        //rotate normal
        x=normal[0]*cos_r-normal[1]*sin_r;
        y=normal[0]*sin_r+normal[1]*cos_r;
        normal.set(x,y);
    }
    */

///////////////////////////////////////////////////////////////////////////////

// clean up the space allocated
void c_ply::destroy()
{
    if (head == NULL) return;
    ply_vertex * ptr = head;
    do
    {
        ply_vertex * n = ptr->getNext();
        delete ptr;
        ptr = n;
    } while (ptr != head);
    head = tail = NULL;

    all.clear();
    triangulation.clear();
}

// Create a empty polygon
void c_ply::beginPoly()
{
    head = tail = NULL;
    all.clear();
    triangulation.clear();
}

// Add a vertex to the polygonal chian
void c_ply::addVertex(double x, double y, bool remove_duplicate)
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
void c_ply::endPoly(bool remove_duplicate)
{
    if (head != NULL && tail != NULL)
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
void c_ply::doInit()
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
    ply_vertex * ptr = head;
    do
    {
        ptr->computeExtraInfo();
        ptr = ptr->getNext();
    } while (ptr != head);
}

const Point2d & c_ply::getCenter()
{
    if (radius <= 0)
    {
        this->center.set(0, 0);
        ply_vertex * ptr = head;
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

const Point2d c_ply::getCenter() const
{
    if (radius < 0)
    {
        Point2d center;
        ply_vertex * ptr = head;
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
void c_ply::reverse()
{
    ply_vertex * ptr = head;
    do
    {
        ptr->reverse();
        ptr = ptr->getNext();
    } while (ptr != head); //end while

    this->area = -this->area;
    all.clear();
    triangulation.clear();
}

float c_ply::getRadius() const
{
    if (radius <= 0)
    {
        double radius = 0;
        getCenter();
        ply_vertex * ptr = head;
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

float c_ply::getArea()
{
    if (area == -FLT_MAX)
    {
        area = 0;
        getCenter();
        ply_vertex * ptr = head;
        do
        {
            ply_vertex * next = ptr->getNext();
            const Point2d & p1 = ptr->getPos();
            const Point2d & p2 = next->getPos();
            area += Area(p1.get(), p2.get(), center.get());
            ptr = next;
        } while (ptr != head); //end while
    }

    return area;
}

void c_ply::build_all()
{
    uint vid = 0;
    all.clear();
    ply_vertex * ptr = head;
    do
    {
        ptr->setVID(vid++);
        all.push_back(ptr);
        ptr = ptr->getNext();
    } while (ptr != head); //end while
}

//
// Compute the center and the box of a list of plys
//

void c_plylist::buildBoxAndCenter()
{
    if (is_buildboxandcenter_called == false)
    {
        //typedef list<c_ply>::iterator IT;
        box[0] = box[2] = FLT_MAX;
        box[1] = box[3] = -FLT_MAX;
        for (iterator i = begin(); i != end(); i++)
        {

            ply_vertex * ptr = i->getHead();
            do
            {
                const Point2d & p = ptr->getPos();
                if (p[0] < box[0]) box[0] = p[0];
                if (p[0] > box[1]) box[1] = p[0];
                if (p[1] < box[2]) box[2] = p[1];
                if (p[1] > box[3]) box[3] = p[1];
                ptr = ptr->getNext();
            } while (ptr != i->getHead()); //end while
        }

        center[0] = (box[0] + box[1]) / 2;
        center[1] = (box[2] + box[3]) / 2;

        is_buildboxandcenter_called = true;
    }
}

void c_polygon::destroy()
{
    for (iterator i = begin(); i != end(); i++)
    {
        i->destroy();
    }
    clear(); //remove all ply from this list
    all.clear();
    triangulation.clear();
}

void c_polygon::build_all()
{
    uint vid = 0;
    for (iterator i = begin(); i != end(); i++)
    {
        uint vsize = i->getSize();
        for (uint j = 0; j < vsize; j++)
        {
            (*i)[j]->setVID(vid++); //id of vertex in the polygons
            all.push_back((*i)[j]);
        }
    } //end for i
}

double c_polygon::getArea()
{
    if (area == 0)
    {
        for (iterator i = begin(); i != end(); i++)
        {
            area += i->getArea();
        } //end for i
    }

    return area;
}

istream & operator>>(istream & is, c_ply & poly)
{
    int vsize;
    string str_type;
    is >> vsize >> str_type;

    if (str_type.find("out") != string::npos)
        poly.type = c_ply::POUT;
    else
        poly.type = c_ply::PIN;

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
        //double d=x*x+y*y;
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

istream & operator>>(istream & is, c_plylist & p)
{
    //remove header comments
    do
    {
        char tmp[1024];
        char c = is.peek();
        if (isspace(c))
            is.get(c); //eat it
        else if (c == '#')
        {
            is.getline(tmp, 1024);
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
        c_ply poly(c_ply::UNKNOWN);
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
