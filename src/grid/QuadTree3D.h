/*! \file QuadTree3D.h

\brief QuadTree 3D grid

\author <a href="masc.cs.gmu.edu/">MASC group</a>, George Mason University,
        <a href="profile.usgs.gov/langevin/">Christian Langevin</a>, USGS
\bug    No known bugs.

 */

//   $Id: $

#pragma once
#include <iostream>
#include "PriorityQueue.h"
#include "ModflowGrid.h"
#include "Util.h"
#include "Point.h"
#include "lineseg.h"
#include "polygon.h"
#include "polyline.h"
#include "GhostNode.h"
#include <vector>
using namespace std;

namespace cusg
{

	//forward declaration in cusg namespace

	struct GhostNode;

	class ascii_grid; //defined in ascii_grid.h

/*! \brief QuadTree 3D grid

    A QuadTree 3D grid is a subclass of Grid.  
    A QuadTree grid contains layers of Quadtrees.

*/

class QuadTree3D : public Grid
{
public:

    /*!
        Initialize the quadtree grid with a parent grid.  The parent grid should
        be a ModflowGrid object.
	*/
	
    QuadTree3D(ModflowGrid * mfgrid, int * chunksize=NULL);
	
	virtual ~QuadTree3D(void);

	//compute the intersection
	//void intersection(int layer=0);

	//
	void smooth_refinement(int level_tol_horizontal, int level_tol_vertical);

    //
	void refine( Point2d& pt,     int max_level, int at_layer=0);
    void refine( LineSeg2d line,  int max_level, int at_layer=0);
    void refine( c_plyline& arc,  int max_level, int at_layer=0);
    void refine( c_ply& ply,      int max_level, int at_layer=0);
	void refine( c_polygon& pgon, int max_level, int at_layer=0);

    //void refine( max_level, at_layer, Point2d);

	/// Return a numpy array of size (nodes) that contains the area for all cells.
	double * get_cell_areas();

        
	/*!
        \brief Return the four vertices for the specified nodeid.
    */    
    /*
	vector<Point2d> get_vertices(int nodeid)
	{
		
	}
	*/

	/*! 
        \brief Return a tuple containing two points (the lower left and the upper 
        right) in global grid coordinates.
    */    	
	pair<Point2d,Point2d> get_extent();
	
	/*
	///Return the nodeid for the modflow grid.
	int get_nodeid(int i, int j)
	{
        //return i * ncol + j;
    }
    */

	/*! 
        Return an array of size nodes that contains the nodeid for all
        active leaf nodes.  Thus the array contains:
        
        nodeid_array[nodenumber] = nodeid
        
        Note that this array is rebuilt, but only if the
        nodegroup.rebuild_nodeid_array flag is True.  This flag is set to
        true anytime quadtree nodes are added or removed.
    */  
    int * getnodeid_array();
        
    ///get nodelay_array
    int * get_nodelay();
	/*!
	    \brief Create a nodenumber_array of size nodegroup with the following:
        
        node number is the id of the node in CSR or Coo array

        nodenumber_array[nodeid] = nodenumber
     */	
	int * createnodenumber_array();
	
	/*!
        Create an array of size modflowgrid.nlay that has the number of
        nodes in each layer.
     */
    void create_nodelay_array();
	
	/*! 
	    Return an array of len(nodegroup) that contains the nodenumber for all
        active leaf nodes.  Thus the array contains:
        
        nodenumber_array[nodeid] = nodenumber
        
        Note that this array is rebuilt, but only if the
        nodegroup.rebuild_nodeid_array flag is True.  This flag is set to
        true anytime quadtree nodes are added or removed.
    */
    int * getnodenumber_array();
            
	///get box from id
    Box * get_nodeobj(int nodenumber);
    
    ///get box from id
    void getnode(void * node, int nodenumber);
    
    ///get number of nodes
    int get_nodes();
    
    int get_nodenumber(int nodeid);
    
    void set_one_based_node_numbering(bool flag){ one_based_node_numbering=flag; }
    
	bool get_one_based_node_numbering(){
		return one_based_node_numbering;
	}

    ///get box connections from id
    ///Return a sorted list of connections, where a connection is
    ///a list [tonodeid, fldir].  fldir is -1, +1, -2, +2 for -x, +x, -y, +y, -z, +z.
	///notice that this function takes "nodenumber" not "nodeid"
	///remember again that the return is for node id.

	list<pair<int, int> > get_node_connections(int nodenumber);
	list<pair<int, int> > get_node_connections(Box* box);

    list< pair<int,int> > get_node_connections_vp(int nodenumber, vector<float>* area );

    //return the
    list< pair<int,int> > get_node_connections_vp(Box * box , vector<float>* area);
	
	/// direction indicator of size nja (+-1 for x, +-2 for y +-3 for z)
	int boundayNodes(vector<int>& nodeids, int dir);


    /*! \brief Obtain CSRData

        Return a UsgCsrData object, which contains all of the necessary
        information for creating an unstructured grid model.

        If top[nodes] and bot[nodes] are provided they are used in the
        calculation of fahl.

     */
    CSRData get_usg_csr_data();

	//obtain CSRData with vertical pass
	CSRData get_usg_csr_data_vp();

    ///get the underline ModflowGrid
    ModflowGrid * getModflowGrid(){ return m_mfgrid; }
    

    ///CREATE top coordinates
	///Do NOT call this function if you just wanted to get the top array
    double * get_top();

    /// CREATE bottom coordinates
	/// Do NOT call this function if you just wanted to get the bot array
    double * get_bot();

	///align the elevation (z) between the shared surface (bottom of top cell and top of bottom cells)
	bool align_elev();

	///recompute z value and dz value of each cell based on the top and bottom
	bool recompute_z_of_cells();

    ///number the nodes (i.e. leaves) in the grid
    void number_nodes();

	///save the Ghost Node to file
	void saveGhostNode(string filePath);

	///refresh the tree to get the new node id and number arrays
	void refresh();

public:

	///output >>
    friend ostream & operator<<(ostream& out, QuadTree3D &grid);

protected:

    void number_nodes(Box * box, unsigned int& number);
	
	void buildConnectionBetweenLayers();
	
	void _add_root_nodes();
    
    ///fldir is -1, +1, -2, +2, -3, +3 for -x, +x, -y, +y, -z, +z.
    void _add_2d_rootnode_connections();
	
    void _add_rootnode_connections();
    
    ///sort the given connection by the node number
	void sortconnectionbynodenumber(list< pair<int,int> > & connectionlist);

	//sort the area array based on connection id
	void sortareabyconid(const list<pair<int, int> >& connectionlist, vector<float>& area);
	
	bool smooth_refinement(Box * box, vector< pair<int,Box*> >& box_pq, int level_tol_horizontal, int level_tol_vertical);

	//linear interploation the value of the center of the box
	double bilinear_interpolation(Box * box, double * grid);
	
	//linear interploation the value of the center of the box
	double bilinear_interpolation(Box * box, const string & arcinfo_filename, bool useAreaWeighted = true);

	//linear interploation the value of a given point
	double bilinear_interpolation(const Point2d& pt, double * zgrid);
	
	//linear interploation the value of a given point 
	double bilinear_interpolation(const Point2d& pt, const string & arcinfo_filename);

	///align the elevation (z) between the shared surface (bottom of top cell and top of bottom cells) of the given box
	void align_elev(Box * top_box);

public:
	
	//linear interploation the value of the center of the box
	//Instead of having this as static, we should move it to modflow grid...
	static double bilinear_interpolation(const Point2d& pt, ascii_grid * agrid);
	
	//linear interploation the value of a given point in a given layer (either top or bottom)
	double bilinear_interpolation(const Point2d& pt, int layer, bool top);
	
protected:

	///modflow grid
	ModflowGrid * m_mfgrid; 
	int * m_chunksize;
	
	int * nodenumber_array;
	int * nodeid_array;
	int * nodelay_array; //? used, there is a same member in nodegroup...
	bool rebuild_nodenumber_array;
	bool rebuild_nodeid_array;

	bool one_based_node_numbering;

	vector<GhostNode> m_gns;//store the ghost nodes

public:

	enum LAYER_Z_VALUE_OPEARTOR { Z_VALUE_REPLICATE_OPEARTOR, Z_VALUE_LINEAR_INTERPOLATION_OPEARTOR};
	void setTopOperators( int i, LAYER_Z_VALUE_OPEARTOR op ){ m_top_layer_z_operator[i]=op; }
	void setBotOperators( int i, LAYER_Z_VALUE_OPEARTOR op ){ m_bot_layer_z_operator[i]=op; }
	LAYER_Z_VALUE_OPEARTOR getTopOperators( int i ){ return m_top_layer_z_operator[i]; }
	LAYER_Z_VALUE_OPEARTOR getBotOperators( int i ){ return m_bot_layer_z_operator[i]; }

	void setTopSource( int i, const string& source ){ m_top_layer_z_source[i]=source; }
	void setBotSource( int i, const string& source ){ m_bot_layer_z_source[i]=source; }
	string getTopSource( int i )
	{ 
		if( m_top_layer_z_source.find(i)==m_top_layer_z_source.end() ) return string(""); //not defined
		return m_top_layer_z_source[i]; 
	}
	string getBotSource( int i )
	{ 
		if( m_bot_layer_z_source.find(i)==m_bot_layer_z_source.end() ) return string(""); //not defined
		return m_bot_layer_z_source[i]; 
	}

	void setTopAreaWeighted( int i, bool area_weighted ){ m_top_area_weighted_interpolate[i]=area_weighted; }
	void setBotAreaWeighted( int i, bool area_weighted ){ m_bot_area_weighted_interpolate[i]=area_weighted; }

protected:

	vector<LAYER_Z_VALUE_OPEARTOR> m_top_layer_z_operator; //default is replicate (Z_VALUE_REPLICATE_OPEARTOR)
	vector<LAYER_Z_VALUE_OPEARTOR> m_bot_layer_z_operator;

	//source name for each layer
	map<int,string> m_top_layer_z_source;
	map<int,string> m_bot_layer_z_source;
	map<int, bool> m_top_area_weighted_interpolate;
	map<int, bool> m_bot_area_weighted_interpolate;

	//map filename to ascii grid
	map<string, ascii_grid *> m_str2ag;

	public:
		bool vertical_pass_through;

		
};

ostream & operator<<(ostream& out, QuadTree3D &grid);

}//end of namespace

