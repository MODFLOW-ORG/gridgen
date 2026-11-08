/*! \file QuadTreeBuilder.h

\brief QuadTree Builder

\author <a href="masc.cs.gmu.edu/">MASC group</a>, George Mason University,
        <a href="profile.usgs.gov/langevin/">Christian Langevin</a>, USGS
\bug    No known bugs.

 */

//   $Id: $

#pragma once

#include "QuadTree3D.h"
#include "def_parser.h"

namespace cusg
{

class QuadtreeBuilder
{
public:
    ///build quadtree from QuadtreeBuilderBlock
    void build(QuadTree3D * qtree, QuadtreeBuilderBlock & rawdata);

protected:
    ///analyze the thickness of the top and bottom elevations
    void analyze_elevation(QuadTree3D * qtree, double * top_elev, double * bot_elev);
};

} //namespace cusg
