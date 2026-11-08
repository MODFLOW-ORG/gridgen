/*! \file Vector.h

\brief Vector Descriptions and Operations
\author Jyh-Ming Lien
\author <a href="http://masc.cs.gmu.edu">MASC Group</a>, George Mason University
\copyright 2010-2012 by Jyh-Ming Lien and George Mason University
*/

#pragma once
#ifndef MATHTOOL_VECTOR
#define MATHTOOL_VECTOR

#include "Basic.h"

namespace cusg
{

/*! \brief Vector Description and Operations
        \tparam T   The type of data stored in the factor (e.g., \c real)
        \tparam D   The dimension of the vector (default: 3)
    */
template <class T, int D = 3> class Vector
{
public:
    Vector(const T v[D]);
    Vector(const Vector & v);
    Vector(const T & x = 0, const T & y = 0, const T & z = 0, const T & w = 0);

    T & operator[](int i) { return m_v[i]; }
    const T & operator[](int i) const { return m_v[i]; }

    T norm() const;           //!< magnitude of vector
    T normsqr() const;        //!< magnitude squared
    Vector normalize() const; //!< normalize

    void set(const Vector & other)
    {
        for (int i = 0; i < D; i++) m_v[i] = other.m_v[i];
    }

    void set(const T other[D])
    {
        for (int i = 0; i < D; i++) m_v[i] = other[i];
    }

    void set(const T & x = 0, const T & y = 0, const T & z = 0, const T & w = 0)
    {
        if (D > 0) m_v[0] = x;
        if (D > 1) m_v[1] = y;
        if (D > 2) m_v[2] = z;
        if (D > 3) m_v[3] = w;
    }

    void get(T other[D]) const
    {
        for (int i = 0; i < D; i++) other[i] = m_v[i];
    }

    const T * get() const { return m_v; }

    T * get() { return m_v; }

    // Vector operator prototypes
    const Vector & operator=(const Vector & v2); //!< assignment
    Vector operator-() const;                    //!< unary negation
    Vector operator+(const Vector & v2) const;   //!< vector add
    Vector operator-(const Vector & v2) const;   //!< vector sub
    Vector operator*(const T & s) const;         //!< scalar multiply
    Vector operator/(const T & s) const;         //!< division by scalar

    Vector operator^(const Vector & v2) const; //!< component *
    T operator*(const Vector & v2) const;      //!< dot product
    Vector operator%(const Vector & v2) const; //!< cross product
    bool operator==(const Vector & two) const; //!< equality

protected:
    T m_v[D];
};

//Implementation
template <class T, int D> Vector<T, D>::Vector(const T v[D])
{
    set(v);
}

template <class T, int D> Vector<T, D>::Vector(const T & x, const T & y, const T & z, const T & w)
{
    set(x, y, z, w);
}

template <class T, int D> Vector<T, D>::Vector(const Vector<T, D> & v)
{
    for (int i = 0; i < D; i++) m_v[i] = v.m_v[i];
}

template <class T, int D> T Vector<T, D>::norm() const // magnitude of vector
{
    return (T)sqrt((float)normsqr());
}

template <class T, int D> T Vector<T, D>::normsqr() const // magnitude squared
{
    T sumsqr = 0;
    for (int i = 0; i < D; i++) sumsqr += SQR(m_v[i]);
    return sumsqr;
}

template <class T, int D> Vector<T, D> Vector<T, D>::normalize() const // normalize
{
    Vector<T, D> newv;
    T magnitude = norm();
    int i;
    for (i = 0; i < D; i++)
    {
        newv.m_v[i] = m_v[i] / magnitude;
        //a zero vector has no direction
        if (fabs(m_v[i]) > magnitude * HUGENUMBER) break;
    }
    for (; i < D; i++)
    {
        newv.m_v[i] = m_v[i] / magnitude;
    }
    return newv;
}

template <class T, int D> const Vector<T, D> & Vector<T, D>::operator=(const Vector<T, D> & v2)
{
    for (int i = 0; i < D; i++) m_v[i] = v2.m_v[i];
    return *this;
}

template <class T, int D> Vector<T, D> Vector<T, D>::operator-() const
{
    Vector<T, D> newv(*this);
    for (int i = 0; i < D; i++) newv.m_v[i] = -newv.m_v[i];
    return newv;
}

template <class T, int D> Vector<T, D> Vector<T, D>::operator+(const Vector & v2) const
{
    Vector<T, D> newv(*this);
    for (int i = 0; i < D; i++) newv.m_v[i] += v2.m_v[i];
    return newv;
}

template <class T, int D> Vector<T, D> Vector<T, D>::operator-(const Vector & v2) const
{
    Vector<T, D> newv(*this);
    for (int i = 0; i < D; i++) newv.m_v[i] -= v2.m_v[i];
    return newv;
}

template <class T, int D> Vector<T, D> Vector<T, D>::operator*(const T & s) const
{
    Vector<T, D> newv(*this);
    for (int i = 0; i < D; i++) newv.m_v[i] *= s;
    return newv;
}

template <class T, int D> Vector<T, D> Vector<T, D>::operator/(const T & s) const
{
    Vector<T, D> newv(*this);
    for (int i = 0; i < D; i++) newv.m_v[i] /= s;
    return newv;
}

template <class T, int D> Vector<T, D> Vector<T, D>::operator^(const Vector<T, D> & v2) const
{
    Vector<T, D> newv(*this);
    for (int i = 0; i < D; i++) newv.m_v[i] *= v2.m_v[i];
    return newv;
}

template <class T, int D> T Vector<T, D>::operator*(const Vector<T, D> & v2) const
{
    T dot = 0;
    for (int i = 0; i < D; i++) dot += (m_v[i] * v2.m_v[i]);
    return dot;
}

template <class T, int D> Vector<T, D> Vector<T, D>::operator%(const Vector<T, D> & v2) const
{
    static_assert(D <= 3, "cannot take cross product of a Vector of more than 3 dimensions");
    if (D == 2) return Vector<T, D>(0, 0);
    //D==3
    Vector<T, D> newv(*this);
    newv.m_v[0] = m_v[1] * v2.m_v[2] - m_v[2] * v2.m_v[1];
    newv.m_v[1] = m_v[2] * v2.m_v[0] - m_v[0] * v2.m_v[2];
    newv.m_v[2] = m_v[0] * v2.m_v[1] - m_v[1] * v2.m_v[0];
    return newv;
}

template <class T, int D> bool Vector<T, D>::operator==(const Vector<T, D> & other) const
{
    for (int i = 0; i < D; i++)
        if (m_v[i] != other.m_v[i]) return false;
    return true;
}

template <class T, int D> Vector<T, D> operator*(T & s, const Vector<T, D> & v)
{
    return v * s;
}

template <class T, int D> ostream & operator<<(ostream & out, const Vector<T, D> & vec)
{
    for (int d = 0; d < D; d++) out << vec[d] << " ";
    return out;
}

template <class T, int D> istream & operator>>(istream & in, Vector<T, D> & vec)
{
    for (int d = 0; d < D; d++) in >> vec[d];

    return in;
}

//Typedef commonly used vector type

/// A 2D vector of reals
using Vector2d = Vector<REAL, 2>;

/// A 3D vector of reals
using Vector3d = Vector<REAL, 3>;

/// A 4D vector of reals
using Vector4d = Vector<REAL, 4>;

} //end of namespace

#endif //MATHTOOL_VECTOR
