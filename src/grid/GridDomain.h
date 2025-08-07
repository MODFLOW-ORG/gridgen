/*! \file GridDomain.h

\brief a class for computing the active/inactive domain

\author <a href="http://masc.cs.gmu.edu">MASC group</a>, George Mason University
\author <a href="http://profile.usgs.gov/langevin/">Christian Langevin</a>, USGS
\bug    No known bugs.
*/

#pragma once
#include "GridIntersection.h"
#include "def_struct.h"

namespace cusg
{

class GridDomain
{
public:
	void build(Grid * grid, int layer, active_domain_raw_data * data);
};


}//end namespace cusg
