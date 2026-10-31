
#pragma once
#ifndef _MATHTOOL_UTILITY
#define _MATHTOOL_UTILITY

#include <math.h>
#include <string.h>
#include <stdlib.h>
#include <string>
#include <iostream> // define C++ stream I/O routines
#include <iomanip>
using namespace std;

#if defined(_WIN32) || defined(_WIN64)
#pragma warning(disable : 4244 4800 4267)
#endif

namespace cusg
{

//data type
#define REAL double

/* range of real numbers */
#define SMALLNUMBER 1.0e-10
#define HUGENUMBER 1.0e10

/* Miscellaneous Scalar Math */
//#define abs(x)      (((x) < 0) ? (-(x)) : (x))
#define sqr(x) ((x) * (x))

inline int round(REAL x, REAL p)
{
    return (long int)(((long int)((x)*pow(10.0, p) + ((x) < 0 ? -0.5 : 0.5))) / pow(10.0, p));
}

inline int round(REAL v)
{
    // Visual C++ won't accept a call to floor with a long int.
    // Calling with long double instead.
#if defined(_MSC_VER)
    int integer = (int)floor((long double)(v));
#else
    int integer = (int)floor((long int)(v));
#endif
    REAL fraction = v - integer;

    if (v > 0)
        return (fraction >= 0.5) ? integer + 1 : integer;
    else
        return (fraction >= -0.5) ? integer : integer + 1;
}

#define sign(x) ((x) >= 0 ? 1 : -1)
//#define swap(x1, x2)  {int tmp=x1; x1=x2; x2=tmp}
#define applysign(x, y) ((y) >= 0 ? abs(x) : -abs(x))

/* Angle Conversions & Constants */

#ifndef PI
#define PI ((REAL)3.14159265358979323846)
#endif

#ifndef PI2
#define PI2 (2 * PI)
#endif

#define RAD2DEG (180 / PI)
#define DEG2RAD (PI / 180)

#define DegToRad(x) ((x) * DEG2RAD)
#define RadToDeg(x) ((x) * RAD2DEG)

#if defined(_WIN32) || defined(_WIN64)

////////////////////////////////////////////////////////////////////////////////////////
// Following functions define M_PI and drand48, which are not starndard c library and
// definitions. In addition, rint used to round off float points to int is also here.
/////////////////////////////////////////////////////////////////////////////////////////

extern "C"
{
    //Implementation of these functions are located in util.cpp
    REAL drand48();
    REAL erand48(register unsigned short * xsubi);
    long irand48(register unsigned short m);
    long krand48(register unsigned short * xsubi, unsigned short m);
    long lrand48();
    long mrand48();
    static void next();
    void srand48(long seedval);
    unsigned short * seed48(unsigned short seed16v[3]);
    void lcong48(unsigned short param[7]);
    long nrand48(register unsigned short * xsubi);
    long jrand48(register unsigned short * xsubi);

    /**Round to closest integer.
          *The rint() function rounds x to an integer value according
          *to the prevalent rounding mode.  The default rounding mode
          *is to round to the nearest integer.
          *@return The  rint() function returns the integer value as a float-
          *ing-point number.
          */
    REAL rint(REAL x);

} //end extern "C"

#endif //_WIN32

} //end of namespace

#endif //_MATHTOOL_UTILITY
