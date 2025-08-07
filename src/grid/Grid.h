
/*! \file Grid.h

\brief A base class for all types of Grids

\author <a href="masc.cs.gmu.edu/">MASC group</a>, George Mason University,
        <a href="profile.usgs.gov/langevin/">Christian Langevin</a>, USGS
\bug    No known bugs.

*/


#pragma once

#if defined(_WIN32) || defined(_WIN64)
#pragma warning( disable : 4355 )
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
	
	Grid():nodegroup(this)
	{
	    is3d=false;
	    sortconnections=true;
	    X=Y=top=bot=NULL;
	    name="anonymous-grid";
	}
	
	virtual ~Grid()
	{
	    if (X!=NULL)   delete [] X;
	    if (Y!=NULL)   delete [] Y;
	    if (bot!=NULL) delete [] bot;
	    if (top!=NULL) delete [] top;
	}

	/*! \brief Obtain CSRData 
	
	    Return a UsgCsrData object, which contains all of the necessary
        information for creating an unstructured grid model.
     */
	virtual CSRData get_usg_csr_data() =0;
	
	///
	/// Return a sorted list of connections, where a connection is
    /// a pair (tonodeid, fldir).  fldir is -1, +1, -2, +2 for -x, +x, -y, +y.  
    ///
	virtual list< pair<int,int> > get_node_connections(int nodeid)=0;
	
	//
	virtual double * get_cell_areas()=0;
	
	//get the box of this grid by node id
	virtual Box * get_nodeobj(int nodeid)=0;

	///this is supposed to return number of nodes in CSR/COO data
	///so for grids such as quadtree/octrees, this returns the number
	///of leaves
	virtual int nodes(){ return nodegroup.node_size(); }

	///number the nodes (i.e. leaves) in the grid
	virtual void number_nodes()=0;

	//
	virtual void deactive_nodes(vector<Box*>& boxes)
	{
	    for(vector<Box*>::iterator i = boxes.begin(); i!=boxes.end(); i++)
	    {
	        (*i)->active=false;
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
	bool sortconnections;
	
    ///all boxes in the grid
    NodeGroup nodegroup;
    
    ///top layer (this is the first of bottom layer...)
    double * top;

    ///all bottom layers (size of nlay+1)
    double * bot;

    ///centers of all boxes
	double * X;
	double * Y;

};


} //namespace cusg

