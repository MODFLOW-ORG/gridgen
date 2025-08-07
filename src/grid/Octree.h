/*! \file Octree.h

\brief A tree of `Cube` nodes

\author <a href="http://masc.cs.gmu.edu/">MASC group</a>, 
        George Mason University,
\author <a href="http://profile.usgs.gov/langevin/">Christian Langevin</a>,
        USGS
\bug    No known bugs.

*/

//   $Id: $

#pragma once
#include <list>
#include "Point.h"
#include "Cube.h"
#include "CubeQueue.h"
#include "Grid.h"
#include "ModflowGrid.h"

namespace cusg
{

/*! \class Octree
    \brief A tree of `Cube` nodes

    A Octree is a tree data structure used to store `Cube` nodes.
    The root of a Octree is a single `Cube`. Each layer divides
    the parent cube into four quadrants: NW, NE, SE and SW. 
    The diagram below shows a Octree of depth 2:
    \verbatim
        -root------                       root          depth 1
        | NW | NE |           |--> NW | NE | SE | SW    depth 2
        -----------
        | SW | SE |
        -----------
    \endverbatim
*/
class Octree : public Grid
{
public:

	CubeQueue* PQ;   //!< a priority queue used to store the tree data structure
	Cube* pRoot;     //!< the root of the Octree
	Vector3d epsilons; //!< minimum width/height/depth of a cube before it will split

    /*!
        Initialize the octree with a cube.
        \param root the root of the new Octree
        \param e    epsilons, the minimum width/height/depth of a cube node in the tree
	*/
	Octree(Cube* root, Vector3d e);
	
    /*!
        Initialize the octree with a cube.
        \param root the root of the new Octree
        \param e    epsilon, the minimum width, height and depth
                    of a cube node in the tree
	*/
    Octree(Cube* root, REAL e);

	virtual ~Octree()
	{

	}
    
    void _add_root_nodes();
    void _add_rootnode_connections();

	///
	/// balance the quadtree so that the depth difference of leaves
	/// is at most 1
	///
	/// this is called "smooth" in the python code
	///
	void balancePhase();

	/*!
	    \brief Create a nodenumber_array of size nodegroup with the following:
        
        node number is the id of the node in CSR or Coo array

        nodenumber_array[nodeid] = nodenumber
     */	
	int *createnodenumber_array();

    /*! \brief Finds the leaf node containing the specified point
     *
     * Starting at parameter `root`, this function searches for the leaf node
     * in which the point \a (x,y,z) lies.
     * \param[in] root  Starting node for the search 
     * \param[in] q     (x,y,z)-coordinates of the point being queried
     * \return A pointer to the node containing the point, or a null reference
     *         if the point cannot be found.
     */
	Cube* getCube(Cube* root, Point3d &q);

    /*! Finds the leaf node containing the specified point
     *
     * Starting at the root of the Octree, this function searches
     * for the leaf node in which the point \a (x,y,z) lies.
     * \param[in] q     (x,y,z)-coordinates of the point being queried
     * \return A pointer to the node containing the point, or a null reference
     *         if the point cannot be found.
     */
	Cube* getCube(Point3d &q) {
		return getCube(pRoot, q);
	}

    /*! \brief Gets box connections, based on the node number.
        \returns A sorted list of connections. A connection is the pair
                `(tonodeid, fldir)`. 
                `fldir` is -1, +1, -2, +2 for -x, +x, -y, +y, -z, +z.
    */
    list< pair<int,int> > get_node_connections(int nodenumber);

	/*! 
        Return an array of size nodes that contains the nodeid for all
        active leaf nodes.  Thus the array contains:
        
        nodeid_array[nodenumber] = nodeid
        
        Note that this array is rebuilt, but only if the
        nodegroup.rebuild_nodeid_array flag is True.  This flag is set to
        true anytime quadtree nodes are added or removed.
    */  
    int * getnodeid_array();

	/*! 
	    Return an array of len(nodegroup) that contains the nodenumber for all
        active leaf nodes.  Thus the array contains:
        
        nodenumber_array[nodeid] = nodenumber
        
        Note that this array is rebuilt, but only if the
        nodegroup.rebuild_nodeid_array flag is True.  This flag is set to
        true anytime quadtree nodes are added or removed.
    */
    int * getnodenumber_array();

    /*! \brief Attempts to expand the Octree by further subdividing it.
        \returns \b True if the tree was expanded; \b false if the tree
                 could not be expanded. For example, if the leaf nodes are
                 already smaller than `epsilon`, the tree will not be expanded.
     */
	bool expand ();
	
	void intersection(int layer=0)
	{
	
	}

	////
	/// iteratively expanding the cubes in the tree
	///
	void subdividePhase();

    /// Returns the height of the octree
    int treeHeight();

    ///number the nodes (i.e. leaves) in the grid
    void number_nodes(){}

protected:
    void initOctree();

    /// The underlying Modflow grid
    ModflowGrid * m_mfgrid; 

	int * nodenumber_array;
	int * nodeid_array;
	int * nodelay_array;
	bool rebuild_nodenumber_array;
	bool rebuild_nodeid_array;

private:

	///insert Node to PQ
    void insertNode(Cube *c);
};

}//end namespace cusg
