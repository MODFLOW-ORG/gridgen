/*! \file Util.cpp

\brief Utility function implementations

\author <a href="http://masc.cs.gmu.edu/">MASC group</a>, George Mason University

*/

#include <cstdlib>
#include <cstring>
#include <stdexcept>
#include <cassert>
#include "Util.h"
using namespace std;

namespace cusg
{
double * clone_array(double * a, int n)
{
    double * a2 = nullptr;

    if (n > 0)
    {
        if (a != nullptr)
        {
            a2 = new double[n];
            assert(a2);
            memcpy(a2, a, sizeof(double) * n);
        }
        else
        {
            throw logic_error("internal error: an array of values to copy has no data");
        }
    }

    return a2;
}

}
