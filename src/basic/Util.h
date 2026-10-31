/*! \file Util.h

\brief Utility function definitions

\author <a href="http://masc.cs.gmu.edu/">MASC group</a>, George Mason University

 */

//$Id: $

#pragma once
#include <cstdlib>
#include <vector>

namespace cusg
{
/**
    \brief Clones a C-style array of doubles.

    Creates a new array, filled with the contents of an existing array.
    A "deep" copy is performed, meaning that each element from the source
    array is copied to the new array. 

    @param[in]  a   The array to clone
    @param[in]  n   The number of elements in parameter 'a' 
    @return         A pointer to the new array. If the input array is 
                    empty, a null reference is returned.
*/
double * clone_array(double * a, int n);

}

// Define types used by SWIG
typedef std::vector<int> VECTORI;
