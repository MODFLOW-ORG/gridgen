
/*!

\brief ...

\author <a href="masc.cs.gmu.edu/">MASC group</a>, George Mason University,
        <a href="profile.usgs.gov/langevin/">Christian Langevin</a>, USGS
\bug    No known bugs.

*/


#pragma once

#include "Grid.h"
#include <Point.h>

namespace cusg
{

/*! \brief 2D Uniform Grid used in Modflow

    ModflowGrid2D is a simple class for representing a two-dimensional
    grid object.
    
*/
class ModflowGrid2D : public Grid
{
public:
	
	/*! \brief Instantiate a ModflowGrid2D object.
        
		\param nrow The number of rows
        \param *ncol* The number of columns.
		\param *delr*: The column spacings along a row.
		\param *delc*: The row spacings along a column.
            
        Optional arguments that can be specified using "...":

        \param *xoffset*: The offset of the grid in the x direction.

        \param *yoffset*: The offset of the grid in the y direction.

        \param *rotation*: The grid rotation angle in degrees relative to the lower  
                      left corner.  DOES NOT WORK YET
                      
        The local coordinate system has a (0, 0) x,y location that corresponds 
        to the lower left corner of the grid.
        
        \note This is a variadic function.
    */
	ModflowGrid2D(int nrow, int ncol, double delr, double delc, double xoffset=0, double yoffset=0, double rotation=0);
	ModflowGrid2D(int nrow, int ncol, double * delr, double * delc, double xoffset=0, double yoffset=0, double rotation=0);
    ModflowGrid2D();
	
	/// # of nodes in all layers (int *)
	int * get_nodelay();
	
	///x coords for the center of all cells
    double * get_local_x_array();
    
    
    /// #  y coord for the center of the cell (and for all cells)
    double * get_local_y_array();
    
	/// x coord of all edges
    double * get_local_Xe_array();
	
	/// y coord of all edges
    double * get_local_Ye_array();
    
    
	/// Return a numpy array of size (nodes) that contains the area for all cells.
	double * get_cell_areas();


        
	/*!
        \brief Return the four vertices for the specified nodeid.

        3-------2
        |       |
        |       |
        0-------1
    */

	vector<Point2d> get_vertices(int nodeid);

	/*! 
        \brief Return a tuple containing two points (the lower left and the upper 
        right) in global grid coordinates.
    */    	
	pair<Point2d,Point2d> get_extent();
	
	///Return the nodeid for the modflow grid.
	int get_nodeid(int i, int j);
    
    ///Return the indices for the specified nodeid. (row major (row_id, col_id))
    Index2d get_indices(int nodeid);

    ///find the index of the node containing pt (row major (row_id, col_id))
    Index2d getIndex(const Point2d& pt);

	int getXIndexByBinarySearch(const Point2d& pt, double* Xe, int left, int right, bool ascending);
	int getYIndexByBinarySearch(const Point2d& pt, double*  Ye, int left, int right, bool ascending);
        
	///get box from id
    ///!!! NOTE This function creates a new BOX
    Box * create_nodeobj(int nodeid);
    

    /// !!! NOTE This function creates a new BOX
    Box * create_nodeobj(int i, int j);
    
    /// return a box containtain the given pt
    Box * find_nodeobj(const Point2d& pt);

	/// return a box containtain the given id
    Box * find_nodeobj(int i, int j);

    ///get box connections from id
    ///Return a sorted list of connections, where a connection is
    ///a list [tonodeid, fldir].  fldir is -1, +1, -2, +2 for -x, +x, -y, +y.  
    list< pair<int,int> > get_node_connections(int nodeid);
    
    
    //	/*!
    //	    Method for intersecting the model grid with a point, line, rectangle,
    //        or polygon.  A rectangle intersection is faster than a polygon
    //        intersection.
    //     */
    //    void intersection()
    //    {
    //		//nothing yet
    //    }

    ///number the nodes (i.e. leaves) in the grid
    virtual void number_nodes();
    
    //access the value of the values in bottom
    double botm(int layer, int row, int col);

	//data
	int nlay;
	int nrow;
	int ncol;
	int nodes;
	
	// DELR is the column widths along a row.  It is of dimension NCOL
	// DELC is the row widths along a column.  It is of dimension NROW.  

	double * delr;
	double * delc;
	double xoffset;
	double yoffset;
	double rotation;
	
	double width;
	double height;

	//x and y coordinate of edges
	double * Xe; //x coord of all edges
	double * Ye; //y coord of all edges

protected:

    void init2D();

    ///create nodes and add them to nodegroup
    void create_nodes();

    ///add connection between boxes
    ///fldir is -1, +1, -2, +2, -3, +3 for -x, +x, -y, +y, -z, +z.
    void add_connections();
    
};


} //namespace cusg

