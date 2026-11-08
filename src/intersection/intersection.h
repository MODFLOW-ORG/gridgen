
/*! \file intersection.h

\brief Computational geometry functions.

Contains various functions for computational geometry in 2D space.
\cite ORourke:1998:CGC:521378

\author Jyh-Ming Lien
\author <a href="http://masc.cs.gmu.edu/">MASC group</a>, George Mason University
\bug    No known bugs.

*/

#pragma once
#ifndef CUSG_INTERSECTION_H
#define CUSG_INTERSECTION_H
#include <Basic.h>
#include <cassert>

namespace cusg
{
// This is copied/modified from CG in C

/*! \brief Copies the values in \a a to \a p
        \param[out] p    The target of the assignment
        \param[in]  a   The pair of values to assign
     */
template <class Real> void assign(Real p[2], Real a[2])
{
    p[0] = a[0];
    p[1] = a[1];
}

/*! \brief Determines if edge ab crosses the ray from p in the +x direction

       Counting the edges of a ring that cross the ray tells whether p is
       inside the ring: an odd count is inside.
     */
inline bool crosses_ray(const double a[2], const double b[2], const double p[2])
{
    if ((a[1] > p[1]) == (b[1] > p[1])) return false;
    double x = a[0] + (p[1] - a[1]) * (b[0] - a[0]) / (b[1] - a[1]);
    return p[0] < x;
}

/*! \brief Returns the area of a triangle

       Returns the area of a 2D triangle (i.e., three sets of vertices
       in two-dimensional space), based on the coordinates.
       \param[in] a         The (x,y) coordinates of the first vertex
       \param[in] b    The (x,y) coordinates of the second vertex
       \param[in] c    The (x,y) coordinates of the third vertex
       \returns The area of the region enclosed by the
       three vertices (may be negative)
     */
template <class Real> Real area(const Real a[2], const Real b[2], const Real c[2])
{
    return ((b[0] - a[0]) * (c[1] - a[1]) - (c[0] - a[0]) * (b[1] - a[1])) / 2;
}

/*! \brief Determines the sign of a triangle's area

       Computes the area of a 2D triangle using the \c area function,
       then returns a value indicating whether the area is positive,
       negative or zero.
       \param[in] a         The (x,y) coordinates of the first vertex
       \param[in] b    The (x,y) coordinates of the second vertex
       \param[in] c    The (x,y) coordinates of the third vertex
       \returns \b 1 if the area is positive, \b -1 if the area
       is negative, or \b 0 if the area is
       zero (or very close to zero).
       \remarks This function takes floating-point imprecision
       into account (e.g., if the area is -1.0E-11, it is considered
       to be zero).

     */
template <class Real> int area_sign(const Real a[2], const Real b[2], const Real c[2])
{
    Real signed_area = area(a, b, c);

    if (signed_area > SMALLNUMBER)
        return 1;
    else if (signed_area < -SMALLNUMBER)
        return -1;
    else
        return 0;
}

/*! \brief Determines if three vertices are collinear

       Determines the area of the triangle created by three
       vertices in 2D space. If and only if the area is zero,
       the points are considered collinear and \b 1 is returned.
       Otherwise, the points are not considered collinear,
       and \b 0 is returned.
       \param[in] a     The (x,y) coordinates of the first vertex
       \param[in] b     The (x,y) coordinates of the second vertex
       \param[in] c     The (x,y) coordinates of the third vertex
       \returns \b 1 if the points are collinear; \b 0 otherwise.
       \remarks Whether or not the area is zero is determined
                via the \c area_sign function.
     */
template <class Real> int collinear(const Real a[2], const Real b[2], const Real c[2])
{
    return area_sign(a, b, c) == 0;
}

/*! \brief Returns \b true iff point \a c lies on the closed segment \a ab
       \param[in] a     The (x,y) coordinates of the first half of line segment \a ab
       \param[in] b    The (x,y) coordinates of the second half of line segment \a ab
       \param[in] c    The (x,y) coordinates of the vertex to check
       \returns \b True if \a c lies between \a a and \a b, or on them; \b false otherwise.
       \remarks Assumes it is already known that \a a and \b b are collinear.
     */
template <class Real> bool between(const Real a[2], const Real b[2], const Real c[2])
{
    // If ab not vertical, check betweenness on x; else on y.
    //if ( a[0]!=b[0] )
    if (fabs(a[0] - b[0]) > SMALLNUMBER)
        return ((a[0] <= c[0]) && (c[0] <= b[0])) || ((a[0] >= c[0]) && (c[0] >= b[0]));
    else
        return ((a[1] <= c[1]) && (c[1] <= b[1])) || ((a[1] >= c[1]) && (c[1] >= b[1]));
}

/*! \brief Returns \b true iff point \a c lies strictly between \a a and \a b
        \param[in] a     The (x,y) coordinates of the first half of line segment \a ab
        \param[in] b    The (x,y) coordinates of the second half of line segment \a ab
        \param[in] c    The (x,y) coordinates of the vertex to check
        \returns \b True iff \a c lies between \a a and \a b.
                    (If \a c is the same vertex as \a a or \a b, this function
                    returns \b false.)
     */
template <class Real> bool between_strict(const Real a[2], const Real b[2], const Real c[2])
{
    // If ab is not vertical enough, check betweenness on x; else on y.
    if (fabs(a[0] - b[0]) > SMALLNUMBER)
    {
        return ((a[0] < c[0]) && (c[0] < b[0])) || ((a[0] > c[0]) && (c[0] > b[0]));
    }
    else
    {
        return ((a[1] < c[1]) && (c[1] < b[1])) || ((a[1] > c[1]) && (c[1] > b[1]));
    }
}

/*! \brief Checks if two 3D points are the essentially the same

       Checks if two points in 3D space are the same, taking floating-point
       imprecision into account (i.e., if the difference between
       their (x,y,z) coordinates are smaller than 1.0E-10, the
       points are considered the same).
       \returns \b True if \a a is essentially equal to \a b; \b false otherwise.
     */
template <class Real> bool almost_equal3(const Real a[3], const Real b[3])
{
    return (fabs(a[0] - b[0]) < SMALLNUMBER && fabs(a[1] - b[1]) < SMALLNUMBER &&
            fabs(a[2] - b[2]) < SMALLNUMBER);
}

/*! \brief Checks if two 2D points are the essentially the same

       Checks if two points in 2D space are the same, taking floating-point
       imprecision into account (i.e., if the difference between
       their (x,y) coordinates are smaller than 1.0E-10, the
       points are considered the same).
       \returns \b True if \a a is essentially equal to \a b; \b false otherwise.
     */
template <class Real> bool almost_equal(const Real a[2], const Real b[2])
{
    return (fabs(a[0] - b[0]) < SMALLNUMBER && fabs(a[1] - b[1]) < SMALLNUMBER);
}

/*! \brief Checks if two 3D points are exactly the same

       Checks if two points in 3D space are the same. Does \b not
       take floating-point imprecision into account.
       \returns \b True if vertex \a a equals vertex \a b; \b false otherwise.
     */
template <class Real> bool equal3(const Real a[3], const Real b[3])
{
    return (a[0] == b[0]) && (a[1] == b[1]) && (a[2] == b[2]);
}

/*! \brief Checks if two 2D points are exactly the same

       Checks if two points in 2D space are the same. Does \b not
       take floating-point imprecision into account.
       \returns \b True if vertex \a a equals vertex \a b; \b false otherwise.
     */
template <class Real> bool equal(const Real a[2], const Real b[2])
{
    return (a[0] == b[0]) && (a[1] == b[1]);
}

/*!
      \brief Compute union of two colinear  segments ab and cd  place the union in p
       return false if the union is degenerated
     */
inline bool segment_union(const double a[2], const double b[2], const double c[2],
                          const double d[2], const double * p[2])
{
    int id = 0;
    if (almost_equal(a, c))
    {
        p[0] = a;

        //can only be b or d
        if (almost_equal(b, d))
        {
            p[1] = b;
            return true;
        }

        //
        else if (between_strict(a, b, d))
        {
            p[1] = d;
            return true;
        }

        //
        else if (between_strict(c, d, b))
        {
            p[1] = b;
            return true;
        }

        //
        return false;
    }

    //ok, now a!=c

    if (almost_equal(a, d)) //a==d
    {
        p[0] = a; //a==d and a!=c, so can only be ab, ac,

        if (almost_equal(b, c))
        {
            p[1] = b;
            return true;
        }

        //
        else if (between_strict(a, b, c))
        {
            p[1] = c;
            return true;
        }

        //
        else if (between_strict(c, d, b))
        {
            p[1] = b;
            return true;
        }

        //
        return false;
    }

    if (almost_equal(b, c))
    {
        p[id] = b; //b==c. a!=c, a!=d

        //can only have a and d

        //
        if (between_strict(a, b, d))
        {
            p[1] = d;
            return true;
        }

        //
        else if (between_strict(c, d, a))
        {
            p[1] = a;
            return true;
        }

        return false;
    }

    if (almost_equal(b, d))
    {
        p[id] = b; //b==d, b!=c. a!=c, a!=d

        //can only have a and c

        if (between_strict(a, b, c))
        {
            p[1] = c;
            return true;
        }

        //
        else if (between_strict(c, d, a))
        {
            p[1] = a;
            return true;
        }

        return false;
    }

    if (between_strict(a, b, c))
    {
        p[0] = c;
        if (between_strict(a, b, d))
        {
            p[1] = d;
            return true;
        }

        if (between_strict(c, d, a))
        {
            p[1] = a;
            return true;
        }

        if (between_strict(c, d, b))
        {
            p[1] = b;
            return true;
        }

        return false;
    }

    if (between_strict(a, b, d))
    {
        p[0] = d;

        if (between_strict(a, b, c))
        {
            p[1] = c;
            return true;
        }

        if (between_strict(c, d, a))
        {
            p[1] = a;
            return true;
        }

        if (between_strict(c, d, b))
        {
            p[1] = b;
            return true;
        }

        return false;
    }

    if (between_strict(c, d, a))
    {
        p[0] = a;

        if (between_strict(a, b, c))
        {
            p[1] = c;
            return true;
        }

        if (between_strict(a, b, d))
        {
            p[1] = d;
            return true;
        }

        if (between_strict(c, d, b))
        {
            p[1] = b;
            return true;
        }

        return false;
    }

    if (between_strict(c, d, b))
    {
        p[0] = b;
        if (between_strict(a, b, c))
        {
            p[1] = c;
            return true;
        }

        if (between_strict(a, b, d))
        {
            p[1] = d;
            return true;
        }

        if (between_strict(c, d, a))
        {
            p[1] = a;
            return true;
        }

        return false;
    }

    return false;
}

/*! \brief Checks if parallel segments intersect

       \param[in] a    The first half of the first line segment, \a ab
       \param[in] b    The second half of the first line segment, \a ab
       \param[in] c    The first half of the second line segment, \a cd
       \param[in] d    The second half of the second line segment, \a cd
       \param[out] p   If the end points overlap, the coordinates of
                       the vertex where they overlap is placed here
       \returns
       - \c 0 if the segments do not intersect
       - \c e if the segments collinearly intersect
       - \c v if an endpoint of one segment is on the other segment
     */
template <class Real>
char parallel_int(const Real a[2], const Real b[2], const Real c[2], const Real d[2], Real p[2])
{
    if (!collinear(a, b, c)) return '0';

    //check if they overlap..
    if (between(a, b, c)) return 'e';
    if (between(a, b, d)) return 'e';
    if (between(c, d, a)) return 'e';
    if (between(c, d, b)) return 'e';

    //they don't overlap but the end points may..
    //check if the end points overlap
    if (almost_equal(a, c))
    {
        p[0] = a[0];
        p[1] = a[1];
        return 'v';
    }
    if (almost_equal(b, c))
    {
        p[0] = b[0];
        p[1] = b[1];
        return 'v';
    }
    if (almost_equal(a, d))
    {
        p[0] = a[0];
        p[1] = a[1];
        return 'v';
    }
    if (almost_equal(b, d))
    {
        p[0] = b[0];
        p[1] = b[1];
        return 'v';
    }

    return '0';
}

/*! \brief Computes the union of two segments

       Computes the union of two collinear segments \a ab and \a cd.
       \param[in] a    The first half of the first line segment, \a ab
       \param[in] b    The second half of the first line segment, \a ab
       \param[in] c    The first half of the second line segment, \a cd
       \param[in] d    The second half of the second line segment, \a cd
       \param[out] p   The union will be placed in this parameter
       \returns \b False if the union is degenerated. Otherwise returns \b true.
     */
template <class Real>
bool segment_union(const Real a[2], const Real b[2], const Real c[2], const Real d[2],
                   const Real * p[2])
{
    int id = 0;
    if (equal(a, c))
    {
        p[id] = a;
        id++;
    }

    if (equal(a, d))
    {
        p[id] = a;
        id++;
    }
    if (id == 2) return p[0] != p[1];

    if (equal(b, c))
    {
        p[id] = b;
        id++;
    }
    if (id == 2) return p[0] != p[1];

    if (equal(b, d))
    {
        p[id] = b;
        id++;
    }
    if (id == 2) return p[0] != p[1];

    if (between_strict(a, b, c))
    {
        p[id] = c;
        id++;
    }
    if (id == 2) return p[0] != p[1];

    if (between_strict(a, b, d))
    {
        p[id] = d;
        id++;
    }
    if (id == 2) return p[0] != p[1];

    if (between_strict(c, d, a))
    {
        p[id] = a;
        id++;
    }
    if (id == 2) return p[0] != p[1];

    if (between_strict(c, d, b))
    {
        p[id] = b;
        id++;
    }
    if (id == 2) return p[0] != p[1];

    return false;
}

/*! \brief Determines if segments ab and cd intersect. */
template <class Real>
bool seg_seg_int(const Real a[2], const Real b[2], const Real c[2], const Real d[2])
{
    //check X and Y coordinates
    for (int i = 0; i < 2; i++)
    {
        if (a[i] < b[i])
        {
            if (c[i] < d[i])
            {
                if (a[i] > d[i]) return false;
                if (c[i] > b[i]) return false;
            }
            else
            { //c[i]>=d[i]
                if (a[i] > c[i]) return false;
                if (d[i] > b[i]) return false;
            }
        }
        else
        { //a[i]>=b[i]
            if (c[i] < d[i])
            {
                if (b[i] > d[i]) return false;
                if (c[i] > a[i]) return false;
            }
            else
            { //c[i]>=d[i]
                if (b[i] > c[i]) return false;
                if (d[i] > a[i]) return false;
            }
        }
    } //end for i

    //OK potential intersection
    int abc = area_sign(a, b, c);
    int abd = area_sign(a, b, d);

    if (abc == 0 && abd == 0)
    { //collinear
        //check if they overlap..
        if (between(a, b, c)) return true;
        if (between(a, b, d)) return true;
        if (between(c, d, a)) return true;
        if (between(c, d, b)) return true;
        return false;
    }
    else if (abc == 0)
    {
        if (between(a, b, c)) return true;
        return false;
    }
    else if (abd == 0)
    {
        if (between(a, b, d)) return true;
        return false;
    }
    //
    else
    { // if(abc!=0 && abd!=0)
        if (abc > 0 && abd > 0) return false;
        if (abc < 0 && abd < 0) return false;
    }

    int cda = area_sign(c, d, a);
    int cdb = area_sign(c, d, b);

    assert(cda != 0 || cdb != 0);

    if (cda == 0)
    {
        if (between(c, d, a)) return true;
        return false;
    }
    else if (cdb == 0)
    {
        if (between(c, d, b)) return true;
        return false;
    }
    else
    {
        if (cda > 0 && cdb > 0) return false;
        if (cda < 0 && cdb < 0) return false;
    }

    return true;
}

/*! \brief Finds the point of intersection between two segments

       Finds the point of intersection \a p between two closed
       segments \a ab and \a cd
       \param[in] a    The first half of the first line segment, \a ab
       \param[in] b    The second half of the first line segment, \a ab
       \param[in] c    The first half of the second line segment, \a cd
       \param[in] d    The second half of the second line segment, \a cd
       \param[out] p   The point of intersection will be placed in this parameter
       \returns
           - \c e:  The segments collinearly overlap, sharing a point.
           - \c v:  An endpoint (vertex) of one segment is on the other segment,
                    but \c e doesn't hold.
           - \c 1:  The segments intersect properly (i.e., they share a point and
                    neither \c v nor \c e holds).
           - \c 0:  The segments do not intersect (i.e., they share no points).
       \remarks Two collinear segments that share just one point, an endpoint
                of each, returns \c e rather than \c v as one might expect.
     */
template <class Real>
char seg_seg_int(const Real a[2], const Real b[2], const Real c[2], const Real d[2], Real p[2])
{
    Real s, t;                // The two parameters of the parametric eqns.
    Real num_s, num_t, denom; // Numerator and denominator of equations.
    char code = '?';          // Return char characterizing intersection.

    const Real small_number = SMALLNUMBER;

    //
    if (a[0] == c[0] && a[1] == c[1])
    {
        p[0] = a[0];
        p[1] = a[1];
        return 'v';
    }
    if (b[0] == c[0] && b[1] == c[1])
    {
        p[0] = b[0];
        p[1] = b[1];
        return 'v';
    }
    if (a[0] == d[0] && a[1] == d[1])
    {
        p[0] = a[0];
        p[1] = a[1];
        return 'v';
    }
    if (b[0] == d[0] && b[1] == d[1])
    {
        p[0] = b[0];
        p[1] = b[1];
        return 'v';
    }
    //

    denom =
        a[0] * (d[1] - c[1]) + b[0] * (c[1] - d[1]) + d[0] * (b[1] - a[1]) + c[0] * (a[1] - b[1]);

    // If denom is zero, then segments are parallel: handle separately.
    if (denom == 0) return parallel_int(a, b, c, d, p);

    if (fabs(denom) < small_number)
    {
        //ok, we are using double, may need more precision
        if (sizeof(a[0]) <= sizeof(double) && sizeof(REAL) != sizeof(double))
        {
            REAL precise_a[2];
            precise_a[0] = a[0];
            precise_a[1] = a[1];
            REAL precise_b[2];
            precise_b[0] = b[0];
            precise_b[1] = b[1];
            REAL precise_c[2];
            precise_c[0] = c[0];
            precise_c[1] = c[1];
            REAL precise_d[2];
            precise_d[0] = d[0];
            precise_d[1] = d[1];
            REAL precise_p[2];
            char r = seg_seg_int<REAL>(precise_a, precise_b, precise_c, precise_d, precise_p);
            p[0] = precise_p[0];
            p[1] = precise_p[1];
            return r;
        }
    }

    //compute s
    num_s = a[0] * (d[1] - c[1]) + c[0] * (a[1] - d[1]) + d[0] * (c[1] - a[1]);

    s = num_s / denom;

    if (fabs(s) < small_number)
        s = 0;
    else if (fabs(1 - s) < small_number)
        s = 1;

    if (s < 0 || s > 1) return '0';

    //compute t
    num_t = -(a[0] * (c[1] - b[1]) + b[0] * (a[1] - c[1]) + c[0] * (b[1] - a[1]));

    t = num_t / denom;

    if (fabs(t) < small_number)
        t = 0;
    else if (fabs(1 - t) < small_number)
        t = 1;

    if (t < 0 || t > 1) return '0';

    //decide the code
    if ((0.0 < s) && (s < 1) && (0.0 < t) && (t < 1))
        code = '1';
    else
        code = 'v';

    if (code != 'v')
    {
        p[0] = (a[0] + s * (b[0] - a[0]));
        p[1] = (a[1] + s * (b[1] - a[1]));
    }
    else
    {
        if (s == 0)
        {
            p[0] = a[0];
            p[1] = a[1];
        }
        else if (s == 1)
        {
            p[0] = b[0];
            p[1] = b[1];
        }
        else if (t == 0)
        {
            p[0] = c[0];
            p[1] = c[1];
        }
        else if (t == 1)
        {
            p[0] = d[0];
            p[1] = d[1];
        }
        else
            assert(false);
    }

    return code;
}
} //namespace cusg

#endif //CUSG_INTERSECTION_H
