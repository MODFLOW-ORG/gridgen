
/*! \file Grid.h

\brief A base class for all types of Grids

\author <a href="masc.cs.gmu.edu/">MASC group</a>, George Mason University,
        <a href="profile.usgs.gov/langevin/">Christian Langevin</a>, USGS
\bug    No known bugs.

*/

#pragma once

#if defined(_WIN32) || defined(_WIN64)
#pragma warning(disable : 4355)
#endif

#include "Box.h"
#include "CsrData.h"
#include "NodeGroup.h"

namespace cusg
{

/*! \brief A base class for all types of Grids

    This is the grid superclass from which the other grid classes inherit
    from.  It contains methods that work with the ModflowGrid and
    QuadTree classes.
    
*/
class Grid
{
public:
    Grid() : node_group(this)
    {
        is3d = false;
        sort_connections = true;
        x_centers = y_centers = top = bot = nullptr;
        name = "anonymous-grid";
    }

    virtual ~Grid()
    {
        if (x_centers != nullptr) delete[] x_centers;
        if (y_centers != nullptr) delete[] y_centers;
        if (bot != nullptr) delete[] bot;
        if (top != nullptr) delete[] top;
    }

    /*! \brief Obtain CSRData 
	
	    Return a UsgCsrData object, which contains all of the necessary
        information for creating an unstructured grid model.
     */
    virtual CSRData get_usg_csr_data() = 0;

    ///
    /// Return a sorted list of the connections of cell, where a connection is
    /// a pair (to_box_index, fldir).  fldir is -1, +1, -2, +2 for -x, +x, -y, +y.
    ///
    virtual list<pair<int, int>> get_node_connections(int cell) = 0;

    //
    virtual vector<double> get_cell_areas() = 0;

    //the box of cell, numbered as in the output
    virtual Box * get_cell_box(int cell) = 0;

    ///this is supposed to return number of nodes in CSR/COO data
    ///so for grids such as quadtree/octrees, this returns the number
    ///of leaves
    virtual int nodes() { return node_group.node_size(); }

    ///number the nodes (i.e. leaves) in the grid
    virtual void number_nodes() = 0;

    //
    virtual void deactivate_nodes(vector<Box *> & boxes)
    {
        for (auto & box : boxes)
        {
            box->active = false;
        }
    }

    //
    // data
    //

    ///name of the grid
    string name;

    ///true if there are more than 1 layer
    bool is3d;

    ///should the connection be sorted (used in get_usg_csr_data functions)
    bool sort_connections;

    ///all boxes in the grid
    NodeGroup node_group;

    ///top layer (this is the first of bottom layer...)
    double * top;

    ///all bottom layers (size of nlay+1)
    double * bot;

    ///centers of all boxes
    double * x_centers;
    double * y_centers;
};

} //namespace cusg
