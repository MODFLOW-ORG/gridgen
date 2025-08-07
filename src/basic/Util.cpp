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
double* clone_array(double *a, int n) 
{
    double *a2=NULL;

    if (n>0) 
    {
        if (a!=NULL) 
        {
            a2 = new double[n];
            assert(a2);
            memcpy(a2, a, sizeof(double) * n);
        } 
        else 
        {
            throw logic_error("!Error: clone_array: n > 0, but 'a' points to nothing");
        }
    } 

    return a2;
}

int* clone_array(int *a, int n) 
{
    int *a2=NULL;

    if (n>0) 
    {
        if (a!=NULL) 
        {
            a2 = new int[n];
            assert(a2);
            memcpy(a2, a, sizeof(int) * n);
        } 
        else 
        {
            throw logic_error("!Error: clone_array: n > 0, but 'a' points to nothing");
        }
    } 

    return a2;
}

double *doubleZeros(int n)
{
    double *a;

    if (n>0) { 
        a = new double[n];
        assert(a);
        memset(a, 0, sizeof(double) * n);
    } else {
        a = NULL;
    }

    return a;
}

int *intZeros(int n)
{
    int *a;

    if (n>0) { 
        a = new int[n];
        assert(a);
        memset(a, 0, sizeof(int) * n);
    } else {
        a = NULL;
    }

    return a;
}
}
