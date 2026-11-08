
#pragma once
#if !defined(MATHTOOL_POINT_H)
#define MATHTOOL_POINT_H

#include <math.h>
#include "Vector.h"

namespace cusg
{

//define k-D point
template <class T, int D = 3> class Point
{
public:
    ///////////////////////////////////////////////////////////////////////////////////
    //
    //  Constructors
    //
    ///////////////////////////////////////////////////////////////////////////////////

    Point(const T & x = 0, const T & y = 0, const T & z = 0, const T & w = 0) { set(x, y, z, w); }

    Point(const T v[D]) { set(v); }

    Point(const Point & other) { *this = other; }

    Point(const Vector<T, D> & vec) { set(vec.get()); }

    ///////////////////////////////////////////////////////////////////////////////////
    //
    //  Operators
    //
    ///////////////////////////////////////////////////////////////////////////////////
    const Point & operator=(const Point & other)
    {
        for (int i = 0; i < D; i++) m_v[i] = other.m_v[i];
        return *this;
    }

    bool operator==(const Point & other) const
    {
        for (int i = 0; i < D; i++)
            if (fabs(other.m_v[i] - m_v[i]) > SMALLNUMBER) return false;
        return true;
    }

    bool almost_equ(const Point & other) const
    {
        for (int i = 0; i < D; i++)
            if (fabs(other.m_v[i] - m_v[i]) > 1e-10) return false;
        return true;
    }

    bool operator!=(const Point & other) const { return !(*this == other); }

    Vector<T, D> operator-(const Point & other) const
    {
        Vector<T, D> vec;
        for (int i = 0; i < D; i++) vec[i] = m_v[i] - other.m_v[i];
        return vec;
    }

    Point operator+(const Vector<T, D> & vect) const
    {
        Point newp;
        for (int i = 0; i < D; i++) newp.m_v[i] = m_v[i] + vect[i];
        return newp;
    }

    ///////////////////////////////////////////////////////////////////////////////////
    //
    //  Access
    //
    ///////////////////////////////////////////////////////////////////////////////////

    void set(const T & x = 0, const T & y = 0, const T & z = 0, const T & w = 0)
    {
        if (D > 0) m_v[0] = x;
        if (D > 1) m_v[1] = y;
        if (D > 2) m_v[2] = z;
        if (D > 3) m_v[3] = w;
        for (int i = 4; i < D; i++) m_v[i] = 0;
    }

    void set(const T other[D])
    {
        for (int i = 0; i < D; i++) m_v[i] = other[i];
    } //memcpy(v,other,sizeof(T)*D); }

    void get(T other[D]) const
    {
        for (int i = 0; i < D; i++) other[i] = m_v[i];
    } //memcpy(other,v,sizeof(T)*D); }
    const T * get() const { return m_v; }
    T * get() { return m_v; }

    T & operator[](int i) { return m_v[i]; }
    const T & operator[](int i) const { return m_v[i]; }

private:
    T m_v[D];

public:
    void rotate(double cx, double cy, double angle)
    {
        if (D < 2) return;
        double * v = get();
        double rx = cx + (v[0] - cx) * cos(angle) - (v[1] - cy) * sin(angle);
        double ry = cy + (v[1] - cy) * cos(angle) + (v[0] - cx) * sin(angle);
        v[0] = rx;
        v[1] = ry;
    }
};

template <class T, int D> ostream & operator<<(ostream & out, const Point<T, D> & point)
{
    for (int d = 0; d < D; d++) out << point[d] << " ";
    return out;
}

template <class T, int D> istream & operator>>(istream & in, Point<T, D> & point)
{
    for (int d = 0; d < D; d++) in >> point[d];

    return in;
}

//small helper function to store GIS related data
template <class T, int D> class GisPoint : public Point<T, D>
{
public:
    GisPoint(const T & x = 0, const T & y = 0, const T & z = 0, const T & w = 0)
        : Point<T, D>(x, y, z, w)
    {
        m_id = 0;
    }
    GisPoint(const T v[2]) : Point<T, D>(v) { m_id = 0; }
    GisPoint(const GisPoint<T, D> & other) : Point<T, D>(other) { m_id = other.m_id; }
    GisPoint(const Vector<T, D> & vec) : Point<T, D>(vec) { m_id = 0; }

    int m_id;
};

// Typedef common used vector type
using Point2d = Point<REAL, 2>;
using Point3d = Point<REAL, 3>;
using Point4d = Point<REAL, 4>;

using GisPoint2d = GisPoint<REAL, 2>;
using GisPoint3d = GisPoint<REAL, 3>;
using GisPoint4d = GisPoint<REAL, 4>;

} //end of namespace

#endif // MATHTOOL_POINT_H
