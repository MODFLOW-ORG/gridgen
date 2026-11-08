/*! \file Box.h

\brief a Box in QuadTree

\author <a href="masc.cs.gmu.edu/">MASC group</a>, George Mason University,
        <a href="profile.usgs.gov/langevin/">Christian Langevin</a>, USGS
*/

//$Id: $

#pragma once

#include <time.h>
#include <stdlib.h>
#include <iostream>
#include <assert.h>
#include <math.h>
#include <vector>
#include <list>
#include <set>
#include <float.h>
//#include "GridIntersection.h"
//#include "ModflowGrid.h"

using namespace std;

#define INEPSILON 1e-10

#include <Point.h>
#include <Vector.h>
namespace cusg
{

using Index2d = Point<int, 2>; //!< 2d index
using Index3d = Point<int, 3>; //!< 3d index

class ModflowGrid;
class Box; //defined later in this file

/*!\class BoxIter
 * \brief Iterator for box nodes
 *
 * An iterator that walks through `Box` nodes in a `QuadTree`,
 * starting at a specific node and walking in a specific
 * direction (north, east, south or west).
 *
 */

class BoxIter
{
public:
    BoxIter(const Box * starting_node, int direction);

    Box * first();

    Box * next();

    Box * end();

private:
    const Box * m_start_box;

    int m_direction;
    int m_prev;
    int m_next;
    int m_cross;

    Box * m_neighbor;
}; //class BoxIter

/*!\class Box
 * \brief A node in QuadTree
 *
 * Each node in a \c QuadTree is represented by a \c Box.
 * Boxes have a defined width, height and Cartesian coordinates.
 * Boxes are capable of adaptively splitting at runtime.
 * Once a box is split, it is divided into four separate boxes -- a box
 * for each quadrant (NW, NE, SW, SE). 
 */
class Box
{
public:
    /*! \brief Check if the x coord or y coord of the center of a given box (next_box) is outside the base box (base)
	 * Precondition: none
	 * Postcondition: none
	 * \param[in] base    The base box
	 * \param[in] next_box The box to be checked
	 * \return true if over limit */

    static bool is_over_limit(const Box * base, const Box * next_box);

    /*!  \brief Constructor of the Box
     *
     * Creates a new box with the specified dimensions at
     * the specified coordinates. The box is created as a
     * leaf node (i.e., it will have no children).
     *
     * \param xx  x coordinate of the center of the box
     * \param yy  y coordinate of the center of the box
     * \param w   width of the box
     * \param h   height of the box
     */
    Box(double xx, double yy, double w, double h)
        : depth(0), x(xx), y(yy), z(0), dx(w), dy(h), dz(0), is_leaf(true), p_parent(nullptr),
          dirty(false), number(0), layer(0), locator(0), id(0),
          active(true) //, childmfgrid(NULL), bRegular(true)
    {
        for (auto & i : p_children)
        {
            i = nullptr;
        }

        p_up = p_down = nullptr;
        priority = 0;
    }

    /*!  \brief Constructor of the Box
     *
     * Creates a new box with the specified dimensions at
     * the specified coordinates. The box is created as a
     * leaf node (i.e., it will have no children).
     *
     * \param center The Cartesian coordinates of the center of the box
     * \param dims   A vector specifying the width and height of the box
     */
    Box(const Point3d & center, Vector3d & dims)
        : depth(0), is_leaf(true), p_parent(nullptr), dirty(false), number(0), layer(0), locator(0),
          id(0), active(true) //, childmfgrid(NULL), bRegular(true)
    {
        x = center[0];
        y = center[1];
        z = center[2];

        dx = dims[0];
        dy = dims[1];
        dz = dims[2];

        for (auto & i : p_children)
        {
            i = nullptr;
        }

        p_up = p_down = nullptr;
        priority = 0;
    }

    /// get all neighbors in direction "dir"
    /// 0 = N, 1=E, 2=S, W=3 (for neighbors)
    void get_neighbors(list<Box *> & neighbors, int dir);

    /*!
	 * determine the status of the box (ex: UNKNOWN, IN, OUT, ON)
	 */

    void update_status() {}

    /// This function is not implemented yet. see node.py for detail.
    void remove()
    {
        /* remove the box from a group

            def remove(self):
                '''
                Remove the node.  It cannot be recovered once it is removed.
                '''
                self.node_group.remove_node(self)
                return

                We may implement this some place else...
         */
    }

    /*!  \brief Splits the box into sub-boxes
     *
     * Splits this box off into sub-boxes (one for each quadrant
     * of this box). The new boxes are created as the child nodes of 
     * this box; the box itself remains unchanged, except that it is no
     * longer considered a leaf node.
     * \param epsilon      If the width or height of this box is not
     *                     at least as big as this value, the box will
     *                     not be split, and \b false is returned.
     * \return \b true, unless:
     *     - the width or height of the box is smaller than `epsilon`
     *     - this box is an internal node, so it has already been split
     *     - insufficient memory
     */
    bool split(double epsilon);

    /// Check if a given point (qx,qy) is contained in this box
    bool in(double qx, double qy)
    {
        if (dx > 0 && (qx < x - dx / 2 - INEPSILON || qx > x + dx / 2 + INEPSILON)) return false;
        if (dx < 0 && (qx < x + dx / 2 - INEPSILON || qx > x - dx / 2 + INEPSILON)) return false;
        if (dy > 0 && (qy < y - dy / 2 - INEPSILON || qy > y + dy / 2 + INEPSILON)) return false;
        if (dy < 0 && (qy < y + dy / 2 - INEPSILON || qy > y - dy / 2 + INEPSILON)) return false;
        return true;
    }

    //Given the point is in box, test whether the point is on the boundaries of the box
    bool on_bound(double qx, double qy)
    {
        if (fabs(qx - (x - dx / 2)) < INEPSILON) return true;
        if (fabs(qx - (x + dx / 2)) < INEPSILON) return true;
        if (fabs(qy - (y + dy / 2)) < INEPSILON) return true;
        if (fabs(qy - (y - dy / 2)) < INEPSILON) return true;
        return false;
        //if((fabs(qx - dx/2) < INEPSILON) || (fabs(qx)
    }

    /// find a leaf box containing a given point (qx,qy)
    Box * find(double qx, double qy);

    /*!\brief recursively obtain all leaves of this box
     
        Recursive method for finding all of the leaf nodes of this object.
        The leaves list should be empty when first called, and will then
        be recursively constructed by traversing the tree.  The leaves
        will contain the pointer for all the leaves of this box.

     */

    void get_leaves(list<Box *> & leaves)
    {
        if (is_leaf)
        {
            leaves.push_back(this);
            return;
        }

        for (auto & i : p_children) i->get_leaves(leaves);
    }

    /*! \brief Connects this box to another box.
     *
     * Connects this box to the box pointed to by `nei`. The
     * current box must be a leaf node.
     * \param nei   The new neighbor to connect to this node
     * \param dir   The direction where the neighbor should be placed 
     *              w.r.t. this box:
     *              -1, +1, -2, +2, -3, +3 for -x, +x, -y, +y, -z, +z
     */

    void add_connection(Box * nei, int dir)
    {
        //assert(is_leaf);//comment by Guilin
        switch (dir)
        {
            case -1: p_children[3] = nei; break;
            case 1: p_children[1] = nei; break;
            case -2: p_children[2] = nei; break;
            case 2: p_children[0] = nei; break;
            case -3: p_down = nei; break;
            case 3: p_up = nei; break;
            default: assert(false); break;
        }
    }

    /*!
	 * Refine p_up (this is called when p_up is at the level lower than the box but is not a leaf)
	 * Don't confuse this with split (sometime it's called refine too...)
	 */
    bool refine_p_up();

    /*!
     * Refine p_down (this is called when p_up is at the level lower than the box but is not a leaf)
     * Don't confuse this with split (sometime it's called refine too...)
     */
    bool refine_p_down();

    ///this implements refine_p_up and refine_p_down
    bool refine_p_up_down(Box ** ppbox);

    //-------------------------------------------------------------------------
    //
    // DATA
    //
    //-------------------------------------------------------------------------

    friend class BoxIter;

    /// depth of the box in the tree
    int depth;

    /*! \brief x coordinate of the center of the box
        \remarks You may ask, "Why not store the start coordinate,
                 rather than the center coordinate?"
                 Answer: Using the center coordinates makes it easier
                 to split the box into four quadrants.
    */
    double x;

    /// y coordinate of the center of the box
    double y;

    /// z coordinate of the center of the box
    double z;

    // potentially z value for each box ()
    // we can represent the z surface using a regular grid.

    /// width of the box
    double dx;

    /// height of the box
    double dy;

    /// depth of the box
    double dz;

    ///?
    int priority;

    /// a flag indicating if this box is a leaf in the tree
    bool is_leaf;

    /// which layer this box is in
    int layer;

    /// is this an active box?
    bool active;

    /*! \brief location of this box

        locator based on following numbering scheme:
        \verbatim
        -------------
        |  0  |  1  |
        |-----|-----|
        |  3  |  2  |
        -------------
        \endverbatim
        The locator for a root box is 0.
    */
    int locator;

    /// id of the box
    int id;

    /// number of the box
    long number; //this value is valid only if the box is a leaf

    /*!

        how p_children points to the neighbors...

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
     *  the pointers are used as neighbor pointers
     * where
     *  0 = NW, 1 = NE, 2 = SE, 3 = SW (if p_children are kids)
     *  0 = N, 1=E, 2=S, W=3 (for neighbors)
     */
    Box * p_children[4];

    /*! Pointer to an adjacent box in a higher layer
     * p_up must be at the same depth or  at a lower depth of this box in the quadtree grid
     *
     * p_up does not need to be a leaf if p_up has the same depth as this box
     * p_up must be a leaf if p_up is at a lower depth than this box
     */
    Box * p_up;

    /*! Pointer to an adjacent box in a lower layer
     * p_down must be at the same depth or at a lower depth of this box in the quadtree grid
     *
     * p_down does not need to be a leaf if p_down has the same depth as this box
     * p_down must be a leaf if p_down is at a lower depth than this box
     */
    Box * p_down;

    /// Pointer to the box's parent in quadtree
    Box * p_parent;

    //when this dirty flag is true, this box should be deleted
    bool dirty;

    //the child modflow grid,used in LGR_Grid
    //ModflowGrid * childmfgrid;

    ///the coordinates of the center after rotation
    void rotate_ctr(double cx, double cy, double rotation, double & rotated_x, double & rotated_y);
};
//end of class Box

} //end namespace cusg
