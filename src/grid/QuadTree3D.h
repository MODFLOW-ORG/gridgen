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
#include <map>
#include <memory>
#include <vector>
using namespace std;

namespace cusg
{

//forward declaration in cusg namespace

struct GhostNode;

class AsciiGrid;        //defined in ascii_grid.h
class ArcinfoAsciiGrid; //defined in ascii_grid.h

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

    QuadTree3D(ModflowGrid * mfgrid);

    ~QuadTree3D(void) override;

    //compute the intersection
    //void intersection(int layer=0);

    //
    void smooth_refinement(int level_tol_horizontal, int level_tol_vertical);

    //
    //refine the cells a feature intersects, and return whether it intersects any
    bool refine(Point2d & pt, int max_level, int at_layer = 0);
    bool refine(Polyline & arc, int max_level, int at_layer = 0);
    bool refine(Polygon & pgon, int max_level, int at_layer = 0);

    /// Return a numpy array of size (nodes) that contains the area for all cells.
    vector<double> get_cell_areas() override;

    /*!
        \brief Return the four vertices for the specified box_index.
    */

    /*! 
        Return an array of size nodes that contains the box_index for all
        active leaf nodes.  Thus the array contains:
        
        box_index_array[cell] = box_index
        
        Note that this array is rebuilt, but only if the
        node_group.rebuild_box_index_array flag is True.  This flag is set to
        true anytime quadtree nodes are added or removed.
    */
    int * get_box_index_array();

    ///the number of cells in each layer
    int * get_layer_cell_counts();
    /*!
	    \brief Create a cell_array of size node_group with the following:
        
        a cell is the number of a node in the CSR or COO arrays

        cell_array[box_index] = cell
     */
    int * create_cell_array();

    /*! 
	    Return an array of len(node_group) that contains the cell for all
        active leaf nodes.  Thus the array contains:
        
        cell_array[box_index] = cell
        
        Note that this array is rebuilt, but only if the
        node_group.rebuild_box_index_array flag is True.  This flag is set to
        true anytime quadtree nodes are added or removed.
    */
    int * get_cell_array();

    ///the box of cell
    Box * get_cell_box(int cell) override;

    int get_cell(int box_index);

    void set_one_based_node_numbering(bool flag) { m_one_based_node_numbering = flag; }

    bool get_one_based_node_numbering() { return m_one_based_node_numbering; }

    ///get box connections from id
    ///Return a sorted list of connections, where a connection is
    ///a list [to_box_index, fldir].  fldir is -1, +1, -2, +2 for -x, +x, -y, +y, -z, +z.
    ///the connections are to box indices

    list<pair<int, int>> get_node_connections(int cell) override;
    list<pair<int, int>> get_node_connections(Box * box);

    list<pair<int, int>> get_node_connections_with_pass_through(int cell, vector<double> * area);

    //return the
    list<pair<int, int>> get_node_connections_with_pass_through(Box * box, vector<double> * area);

    /*! \brief Obtain CSRData

        Return a UsgCsrData object, which contains all of the necessary
        information for creating an unstructured grid model.

        If top[nodes] and bot[nodes] are provided they are used in the
        calculation of fahl.

     */
    CSRData get_usg_csr_data() override;

    ///get the underline ModflowGrid
    ModflowGrid * get_modflow_grid() { return m_mfgrid; }

    ///CREATE top coordinates
    ///Do NOT call this function if you just wanted to get the top array
    double * get_top();

    /// CREATE bottom coordinates
    /// Do NOT call this function if you just wanted to get the bot array
    double * get_bot();

    ///align the elevation (z) between the shared surface (bottom of top cell and top of bottom cells)
    bool align_shared_elevations();

    ///recompute z value and dz value of each cell based on the top and bottom
    bool set_cell_elevations();

    ///number the nodes (i.e. leaves) in the grid
    void number_nodes() override;

    ///save the ghost nodes to file, including those of vertical connections if vertical is true
    void save_ghost_node(string file_path, bool vertical = false);

    ///refresh the tree to get the new box index and cell arrays
    void refresh();

public:
    ///output >>
    friend ostream & operator<<(ostream & out, QuadTree3D & grid);

protected:
    void number_nodes(Box * box, unsigned int & number);

    void add_root_nodes();

    ///sort the given connections by cell number
    void sort_connections_by_cell(list<pair<int, int>> & connection_list);

    //sort the area array based on connection id
    void sort_areas_by_cell(const list<pair<int, int>> & connection_list, vector<double> & area);

    bool smooth_refinement(Box * box, vector<pair<int, Box *>> & box_pq, int level_tol_horizontal,
                           int level_tol_vertical);

    //linear interpolation the value of the center of the box
    double bilinear_interpolation(Box * box, double * grid);

    //linear interpolation the value of the center of the box
    double bilinear_interpolation(Box * box, const string & arcinfo_filename,
                                  bool use_area_weighted = true);

    ///align the elevation (z) between the shared surface (bottom of top cell and top of bottom cells) of the given box
    void align_shared_elevations(Box * top_box);

public:
    //linear interpolation the value of the center of the box
    //Instead of having this as static, we should move it to modflow grid...
    static double bilinear_interpolation(const Point2d & pt, AsciiGrid * agrid);

protected:
    ///modflow grid
    ModflowGrid * m_mfgrid;

    vector<int> m_cell_array;
    int * m_box_index_array;
    int * m_layer_cell_counts; //? used, there is a same member in node_group...
    bool m_rebuild_cell_array;
    bool m_rebuild_box_index_array;

    bool m_one_based_node_numbering;

public:
    enum LayerZValueOperator
    {
        Z_VALUE_REPLICATE_OPERATOR,
        Z_VALUE_LINEAR_INTERPOLATION_OPERATOR
    };
    void set_top_operator(int i, LayerZValueOperator op) { m_top_layer_z_operator[i] = op; }
    void set_bottom_operator(int i, LayerZValueOperator op) { m_bot_layer_z_operator[i] = op; }
    LayerZValueOperator get_top_operator(int i) { return m_top_layer_z_operator[i]; }
    LayerZValueOperator get_bottom_operator(int i) { return m_bot_layer_z_operator[i]; }

    void set_top_source(int i, const string & source) { m_top_layer_z_source[i] = source; }
    void set_bot_source(int i, const string & source) { m_bot_layer_z_source[i] = source; }
    string get_top_source(int i)
    {
        if (m_top_layer_z_source.find(i) == m_top_layer_z_source.end())
            return string(""); //not defined
        return m_top_layer_z_source[i];
    }
    string get_bot_source(int i)
    {
        if (m_bot_layer_z_source.find(i) == m_bot_layer_z_source.end())
            return string(""); //not defined
        return m_bot_layer_z_source[i];
    }

    void set_top_area_weighted(int i, bool area_weighted)
    {
        m_top_area_weighted_interpolate[i] = area_weighted;
    }
    void set_bot_area_weighted(int i, bool area_weighted)
    {
        m_bot_area_weighted_interpolate[i] = area_weighted;
    }

protected:
    vector<LayerZValueOperator>
        m_top_layer_z_operator; //default is replicate (Z_VALUE_REPLICATE_OPERATOR)
    vector<LayerZValueOperator> m_bot_layer_z_operator;

    //source name for each layer
    map<int, string> m_top_layer_z_source;
    map<int, string> m_bot_layer_z_source;
    map<int, bool> m_top_area_weighted_interpolate;
    map<int, bool> m_bot_area_weighted_interpolate;

    //map filename to ascii grid
    map<string, unique_ptr<ArcinfoAsciiGrid>> m_str2ag;

public:
    bool vertical_pass_through;
};

ostream & operator<<(ostream & out, QuadTree3D & grid);

} //end of namespace
