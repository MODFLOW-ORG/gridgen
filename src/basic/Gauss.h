
#pragma once

/*
        Gaussian Random Number Generator

    Originally created by
    Donald H. House         July 1, 1982
    conversion to C -- Nov. 30, 1989

  This function takes as parameters real valued mean and standard-deviation,
  and an integer valued seed.  It returns a real number which may be
  interpreted as a sample of a normally distributed (Gaussian) random
  variable with the specified mean and standard deviation.
  After the first call to gauss, the seed parameter is ignored.

  The computational technique used is to pass a uniformly distributed random
  number through the inverse of the Normal Distribution function.

  Your program must #include <math.h>
*/
#include <math.h>
#include "Basic.h"

namespace cusg
{
	REAL gauss(REAL mean, REAL std);

	inline REAL gauss(REAL mean, REAL std, int seed)
	{
		static bool first_time=true;
		if( first_time ){ srand48(seed); first_time=false; }
		return gauss(mean,std);
	}
}
