
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
class c_plyline
{

public:

    ///////////////////////////////////////////////////////////////////////////
    c_plyline(){ head=tail=NULL; }
    c_plyline(const c_ply& ply);

    ///////////////////////////////////////////////////////////////////////////
    void copy(const c_plyline& plyline); //copy from the other ply
    void destroy();

    ///////////////////////////////////////////////////////////////////////////
    // create c_ply
    void beginPoly();
    void addVertex( double x, double y, bool remove_duplicate=false );
    void addVertex( ply_vertex * v );
    void endPoly();

    ///////////////////////////////////////////////////////////////////////////
    //convert plyline to ply
    c_ply toply() const;

    ///////////////////////////////////////////////////////////////////////////
    void negate();
    void translate(const Vector2d& v);
    void rotate(double radius);
    void expand(c_plyline& plyline); //add the given plyline to the end of this plyline

    ///////////////////////////////////////////////////////////////////////////
    // Access
    ply_vertex * getHead() const { return head; }
    ply_vertex * getTail() const { return tail; }

    int getSize() {
        if(all.empty()) build_all();
        return all.size();
    }

    ply_vertex * operator[](unsigned int id){
        if(all.empty()) build_all();
        return all[id];
    }

    ///////////////////////////////////////////////////////////////////////////
    // Operator
    //check if give poly line is the same as this
    bool operator==( const c_plyline& other ){ return other.head==head; }
    friend istream& operator>>( istream&, c_plyline& );
    friend ostream& operator<<( ostream&, c_plyline& );

	//rotate function
	void rotate(double cx, double cy, double angle);

protected:

    ///////////////////////////////////////////////////////////////////////////
    void doInit(); /*return # of vertice in this poly*/

    //build elements in vector<ply_vertex*> all
    void build_all();

    ply_vertex * head; //the head of vertex list
    ply_vertex * tail; //end of the vertex list

    vector<ply_vertex*> all; //all vertices
};

class GIS_plyline : public c_plyline
{
public:

    GIS_plyline() : c_plyline(){ m_id=0; }
    GIS_plyline(const GIS_plyline& other): c_plyline((c_plyline&)other){ m_id=other.m_id; }

    //return the i-th line segment from this polyline
    GIS_LineSeg2d getLineSegment(int i);

    //compute distance for geodesic_dist
    void computeGeodesicDist();

    //store attributes from GIS
    int m_id;
    vector<double> geodesic_dist; ///geodesic distance of each vertex


};

//    //a c_plylist is a list of c_ply
//    class c_plylinelist : public list<c_plyline>
//    {
//        friend ostream& operator<<( ostream&, c_plylinelist& );
//        friend istream& operator>>( istream&, c_plylinelist& );
//
//    public:
//
//        c_plylinelist(){is_buildboxandcenter_called=false;}
//        void negate();
//        void translate(const Vector2d& v);
//        void rotate(double r);
//
//        //access
//        void buildBoxAndCenter();
//        double * getBBox() { assert(is_buildboxandcenter_called); return box; }
//        const Point2d& getCenter() { assert(is_buildboxandcenter_called); return center; }
//
//    protected:
//
//        Point2d center;
//        double box[4];
//
//    private:
//
//        bool is_buildboxandcenter_called;
//    };


} //end of namespace cusg

#endif //_POLYLINE_H_


