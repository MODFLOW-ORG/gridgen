/*! \file Cube.h

\brief a Cube in Octree

\author <a href="masc.cs.gmu.edu/">MASC group</a>, George Mason University,
        <a href="profile.usgs.gov/langevin/">Christian Langevin</a>, USGS
\bug    No known bugs.

*/

//$Id: $

#pragma once

#include <cassert>
#include <cmath>
#include <cstdlib>
#include <iostream>
#include <list>
#include <set>
#include <string>
#include <vector>
#include <cfloat>
#include "Point.h"
#include "Vector.h"

using namespace std;


namespace cusg
{

class Cube; //defined later in this file

/*!\class CubeIter
 * \brief Iterator for cube nodes
 *
 * An iterator that walks through `Cube` nodes in a `QuadTree`,
 * starting at a specific node and walking in a specific
 * direction (north, east, south or west).
 *
 */
 
class CubeIter
{
private:
	const Cube* b;
	int direction;
	int prev;
	int next;
	int cross;
	Cube* neighbor;
    int layer;
public:
	
	CubeIter(const Cube* bb, int direc, int layer);

	Cube* First();

	Cube* Next();

	Cube* End();
}; //class CubeIter

/*!\class Cube
 * \brief A node in an octree
 *
 * Each node in a \c Octree is represented by a \c Cube.
 * Cubes have a defined width, height and Cartesian coordinates.
 * Cubes are capable of adaptively splitting at runtime.
 */
class Cube
{
private:

	/*! \brief Check if the x, y or z coord of the center of a given cube (nextCube) is outside the base cube (base)
	 * Precondition: none
	 * Postcondition: none
	 * \param[in] base    The base cube
	 * \param[in] nextCube The cube to be checked
	 * \return true if over limit */

	static bool isOverLimit(const Cube* base, const Cube* nextCube)
	{
        /* My program crashes if the line below is uncommented...? */
		if ((nextCube->x > base->x - base->dx / 2 && nextCube->x < base->x + base->dx / 2)
			|| (nextCube->y > base->y - base->dy / 2 && nextCube->y < base->y + base->dy / 2)
			/*|| (nextCube->z > base->z - base->dz / 2 && nextCube->z < base->z + base->dz / 2)*/)
        {
			return false;
		}
		return true;
	}

public:
	friend class CubeIter;

    /// number of nodes in `pChildren` items
    static const int numChildren = 4;

    /// number of arrays in `pChildren`
    static const int numLayers = 2;

	/// depth of the cube in the tree
	int depth; 
	
	/*! \brief x coordinate of the center of the cube
        \remarks You may ask, "Why not store the start coordinate, 
                 rather than the center coordinate?" 
                 Answer: Using the center coordinates makes it easier 
                 to split the cube into four quadrants.
    */
	double x;
	
	/// y coordinate of the center of the cube
	double y; 
	
	/// z coordinate of the center of the cube
	double z; 
	
	/// width of the cube
	double dx;
	
	/// height of the cube
	double dy;
	
	/// depth of the cube
	double dz;
	
    /// Allows a priority to be set when the cube is placed 
    /// in a priority queue
    int priority;

	/// a flag indicating if this cube is a leaf in the tree
	bool isLeaf; 
	
	/*!

        how pChildren points to the neighbors...
        
	        |----|
	        | n0 |
	        |----+----|----|
	        |         | n1 |
	   |----|  cell   |----|
	   | n3 |         |
	   |----|---------|
	             | n2 |
	             |____|

	 */
	
	/*! Pointers to children, but when no children (i.e., leaf),
	 *	the pointers are used as neighbor pointers
	 * where
     *  inner index:
	 *	 0 = NW, 1 = NE, 2 = SE, 3 = SW
     *      (if pChildren are kids)
	 *   0 = N, 1 = E, 2 = S, 3 = W
     *      (for neighbors)
     *  outer index:
     *   0 = front, 1 = back
	 */
	Cube* pChildren[numChildren][numLayers]; 
	
    /// Releases all memory allocated for nodes in the tree.
    void deleteNodes();

    // Releases memory allocated for nodes under the specified cube.
    void deleteNodes(Cube *);

	/// get all neighbors in direction "dir"
	/// 0 = N, 1 = E, 2 = S, 3 = W
    /// (for neighbors)
	void getNeighbors(list<Cube*>& neighbors, int dir, int layer);

    /// If this is a leaf node, this gives the node to the north of this cube
	void getNorthNeighbors(list<Cube*>& neighbors);

    /// If this is a leaf node, this gives the node to the east of this cube
	void getEastNeighbors(list<Cube*>& neighbors);

    /// If this is a leaf node, this gives the node to the south of this cube
	void getSouthNeighbors(list<Cube*>& neighbors);

    /// If this is a leaf node, this gives the node to the west of this cube
	void getWestNeighbors(list<Cube*>& neighbors);

    /// Pointer to the cube's parent in quadtree
	Cube* pParent; 

	//when this dirty flag is true, this cube should be deleted
	bool dirty;


    /*!  \brief Constructor of the Cube
     *
     * Creates a new cube with the specified dimensions at
     * the specified coordinates. The cube is created as a
     * leaf node (i.e., it will have no children).
     *
     * \param xx  x coordinate of the center of the cube
	 * \param yy  y coordinate of the center of the cube
     * \param zz  z coordinate of the center of the cube
	 * \param w   width of the cube
	 * \param h   height of the cube
     * \param d   depth of the cube
     */
	Cube(double xx, double yy, double zz, double w, double h, double d):
	    depth(0), x(xx), y(yy), z(zz), dx(w), dy(h), dz(d),
	    isLeaf(true), pParent(NULL), dirty(false)
	{
		for (int i = 0; i < numChildren; ++i)
		{
            for (int j = 0; j < 2; ++j)
            {
			    pChildren[i][j] = NULL;
            }
		}
		
	}

    /*!  \brief Constructor of the Cube
     *
     * Creates a new cube with the specified dimensions at
     * the specified coordinates. The cube is created as a
     * leaf node (i.e., it will have no children).
     *
     * \param center The Cartesian coordinates of the center of the cube
	 * \param dims   A vector specifying the width and height of the cube
     */
	Cube(const Point3d& center, Vector3d& dims):
	    depth(0), isLeaf(true), pParent(NULL), dirty(false)
	{
		x=center[0];
		y=center[1];
		z=center[2];
		
		dx=dims[0];
		dy=dims[1];
		dz=dims[2];
		
		for (int i = 0; i < numChildren; ++i)
        {
            for (int j = 0; j < numLayers; ++j)
		    {
			    pChildren[i][j] = NULL;
		    }
	    }	
	}
	
    virtual ~Cube();

	/*!
	 * determine the status of the cube (ex: UNKNOWN, IN, OUT, ON)
	 */

	void updateStatus()
	{
	
	}

    
     /// This function is not implemented yet.
     void remove()
     {
        assert(false);
     }
        
    /*!  \brief Splits the cube into sub-cubes
     *
     * Splits this cube off into sub-cubes (one for each quadrant
     * of this cube). The new cubees are created as the child nodes of 
     * this cube; the cube itself remains unchanged, except that it is no
     * longer considered a leaf node.
     * \param epsilons      If the width, height or depth of this cube is not
     *                      at least as the corresponding values in this vector, 
     *                      the cube will not be split, and \b false 
     *                      is returned.
     *
     * \return \b true, unless:
     *     - the width, height or depth of the cube is smaller than 
     *       the corresponding value in `epsilons`
     *     - this cube is an internal node, so it has already been split
     *     - insufficient memory
     */
	virtual bool split(Vector3d epsilons);

    /// \brief reverse operation of split
    /// This function is not implemented yet.
	bool coarsen();
	
    /// Check if a given point `q` is contained in this cube
	bool in(Point3d q)
	{
	    if(q[0]<x-dx/2 || q[0]>x+dx/2) return false;
	    if(q[1]<y-dy/2 || q[1]>y+dy/2) return false;
        if(q[2]<z-dz/2 || q[2]>z+dz/2) return false;
	    return true;
	}

    /// find a leaf cube containing a given point `q`
	Cube * find(Point3d q)
	{
	    if(isLeaf && in(q)) return this;

	    for(int i=0;i<numChildren;i++){
            for (int j=0;j<numLayers;j++) {
                if(pChildren[i][j]->in(q))
                    return pChildren[i][j]->find(q);
            }
	    }

	    return NULL;
	}

    /*!\brief recursively obtain all leaves of this cube
     
        Recursive method for finding all of the leaf nodes of this object.
        The leaves list should be empty when first called, and will then
        be recursively constructed by traversing the tree.  The leaves
        will contain the pointer for all the leaves of this cube.

     */
     
	void getLeaves(list<Cube*>& leaves)
	{
	    if(isLeaf){
	        leaves.push_back(this);
	        return;
	    }

        for(int i=0;i<numChildren;i++) {
            for (int j=0;j<numLayers;j++) {
                pChildren[i][j]->getLeaves(leaves);
            }
        }    
	}
	
    /*! \brief Connects this cube to another cube.
     *
     * Connects this cube to the cube pointed to by `nei`. The
     * current cube must be a leaf node.
     * \param nei   The new neighbor to connect to this node
     * \param dir   The direction where the neighbor should be placed 
     *              w.r.t. this cube:
     *              -1, +1, -2, +2, -3, +3 for -x, +x, -y, +y, -z, +z
     */
	void addConnection(Cube * nei, int dir)
	{
		assert(isLeaf);
		switch(dir)
		{
			case -1: pChildren[3][0]=nei; break;
			case  1: pChildren[1][0]=nei; break;
			case -2: pChildren[2][0]=nei; break;
			case  2: pChildren[0][0]=nei; break;
			case -3: assert(false); break;
			case  3: assert(false); break;
			default:
				assert(false); break;
		}
	}

    /*! \brief Compare cube pointers by depth
    
        \returns \b True if the depth of `a` is less than the depth of `b`
    */                 
    static bool comparePtrDepths(const Cube *a, const Cube *b) {
        if (!a || !b) return false;
        return a->depth < b->depth;
    }

    /// Returns the coordinates and dimensions of this cube.
    string str();

protected:
    /// Creates new child nodes within this cube
    virtual Cube **spawn(int back);
};//end of class Cube


} //end namespace cusg
