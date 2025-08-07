
#pragma once
#if !defined(_MATHTOOL_Point_H_)
#define _MATHTOOL_Point_H_

#include <math.h>
#include "Vector.h"

namespace cusg
{

    //define k-D point
    template<class T, int D=3>
    class Point 
    {
    public:

        ///////////////////////////////////////////////////////////////////////////////////
        //
        //  Constructors
        //
        ///////////////////////////////////////////////////////////////////////////////////

        Point(const T& x=0,const T& y=0,const T& z=0, const T& w=0)
        { 
            set(x,y,z,w); 
        }

        Point(const T V[D])
        {
            set(V);
        }

        Point( const Point & other ) { *this = other; }

        Point( const Vector<T,D>& vec ){ 
            set(vec.get());
        }

        ///////////////////////////////////////////////////////////////////////////////////
        //
        //  Operators
        //
        ///////////////////////////////////////////////////////////////////////////////////
        const Point & operator=( const Point & other ) {
            //memcpy(v,other.v,sizeof(T)*D);
            for(int i=0;i<D;i++) v[i]=other.v[i];
            return *this;
        }

        bool operator==( const Point & other ) const {
            for( int i=0;i<D;i++ )
                if( fabs(other.v[i]-v[i])>SMALLNUMBER ) return false;
            return true;
        }

        bool almost_equ( const Point & other ) const {
            for( int i=0;i<D;i++ )
                if( fabs(other.v[i]-v[i])>1e-10 ) return false;
            return true;
        }

        bool operator!=( const Point & other ) const {
            return !(*this==other);
        }

        Vector<T,D> operator-(const Point & other) const {
            Vector<T,D> vec;
            for(int i=0;i<D;i++)
                vec[i]=v[i]-other.v[i];
            return vec;
        }

        Point operator+(const Vector<T,D> & vect) const {
            Point newp;
            for( int i=0;i<D;i++ )
                newp.v[i]=v[i]+vect[i];
            return newp;
        }

        ///////////////////////////////////////////////////////////////////////////////////
        //
        //  Access
        //
        ///////////////////////////////////////////////////////////////////////////////////
        
        void set(const T& x=0,const T& y=0,const T& z=0, const T& w=0) {
            if( D>0 ) v[0]=x; if( D>1 ) v[1]=y;
            if( D>2 ) v[2]=z; if( D>3 ) v[3]=w;
            for(int i=4;i<D;i++) v[i]=0;
        }

        void set(const T other[D]) { for(int i=0;i<D;i++) v[i]=other[i]; } //memcpy(v,other,sizeof(T)*D); }

        void get(T other[D]) const { for(int i=0;i<D;i++) other[i]=v[i]; } //memcpy(other,v,sizeof(T)*D); }
        const T * get() const { return v; }
        T * get() { return v; }
        
        T& operator[]( int i ){ return v[i]; }
        const T& operator[]( int i ) const{ return v[i]; }

    private:
        T v[D];

	public:
		void rotate(double cx, double cy, double angle)
		{
			if(D < 2)
				return;
			double* v = get();
			double rx = cx + (v[0]  - cx) * cos(angle) - (v[1] - cy) * sin(angle);
			double ry = cy + (v[1] - cy) * cos(angle) + (v[0] - cx) * sin(angle);
			v[0] = rx;
			v[1] = ry;
		}

    };

    template<class T, int D>
    ostream & operator<<(ostream & out, const Point<T,D> & point) {
        for( int d=0;d<D;d++ ) out<<point[d]<<" ";
        return out;
    }

    template<class T, int D>
    istream & operator>>(istream & in, Point<T,D> & point) {
        for( int d=0;d<D;d++ ) 
            in>>point[d];

        return in;
    }
    
    //small helper function to store GIS related data
    template<class T, int D>
    class GIS_Point : public Point<T,D>
    {
    public:
		
        GIS_Point(const T& x=0,const T& y=0,const T& z=0, const T& w=0) : Point<T,D>(x,y,z,w){ m_id=0; }
        GIS_Point(const T V[2]) :  Point<T,D>(V){ m_id=0; }
        GIS_Point( const GIS_Point<T,D> & other ) : Point<T,D>(other){ m_id=other.m_id;}
        GIS_Point( const Vector<T,D>& vec ) :  Point<T,D>(vec){m_id=0;}

        int m_id;

    };

    // Typedef common used vector type
    typedef Point<REAL,2> Point2d;
    typedef Point<REAL,3> Point3d;
    typedef Point<REAL,4> Point4d;

    typedef GIS_Point<REAL,2> GIS_Point2d;
    typedef GIS_Point<REAL,3> GIS_Point3d;
    typedef GIS_Point<REAL,4> GIS_Point4d;

}//end of namespace

#endif // _MATHTOOL_Point_H_
