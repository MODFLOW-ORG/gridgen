/*! \file Box.cpp

\brief `Box` and `BoxIter` class implementations

\author <a href="masc.cs.gmu.edu/">MASC group</a>, George Mason University
\author <a href="profile.usgs.gov/langevin/">Christian Langevin</a>, USGS
\bug    No known bugs.

*/

//$Id: $

#include "Box.h"
#include "ModflowGrid.h"
#include "GridgenError.h"
#include <assert.h>

namespace cusg
{

// ----------------------------------------------------------------------------
//
//
// Box Iterator Implementation
//
//
//
//
// ----------------------------------------------------------------------------

/*!

       |----|  |----|
       |end |..| 1st|
       |----+--+----|
       |            |
       |    cell    |
       |            |
       |            |
       |------------|

       when direction=0 (north)
       then, 1st (returned by first() is the cell, on the far right
       and end is the cell on the far left and is also the cell stored in p_children[0]

*/

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

/*! \brief Creates a new box neighbor iterator.
 *
 * Box neighbor iterator constructor. Runs through QuadTree
 * nodes in the specified direction.
 * \param[in] starting_node  The starting node
 * \param[in] direction      The direction to travel:
 *                           - 0: North
 *                           - 1: East
 *                           - 2: South
 *                           - 3: West
 */
BoxIter::BoxIter(const Box * starting_node, int direction)
    : m_start_box(starting_node), m_direction(direction), m_neighbor(nullptr)
{
    m_prev = (direction + 3) % 4;
    m_next = (direction + 1) % 4;
    m_cross = (direction + 2) % 4;
}

/*! \brief Returns the first box in the node walk.
 *  \returns A pointer to the first node in the direction to be
 *           traveled (as specified when the iterator was instantiated).
 *           If there are no nodes adjacent to the given box
 *           in that direction, a null reference is returned.
 */

Box * BoxIter::first()
{
    // Get the first box in the node walk.
    //
    // If 'b' is a leaf node, then p_children points to the neighbors
    // of the node, with 0 = N, 1 = E, 2 = S, 3 = W. So looking up
    // p_children[direction] will give the adjacent box in the direction
    // to walk, if such a box exists.
    //
    // If 'b' is an internal node, then p_children, as the name implies,
    // holds pointers to the child nodes. For p_children: 0 = NW, 1 = NE,
    // 2 = SE, 3 = SW.
    Box * n = m_start_box->p_children[m_direction];
    if (!n)
    {
        // There's no nodes when we walk in this direction,
        // so return a null pointer.
        return nullptr;
    }

    // If 'b' is a leaf node, its depth will be the same as the
    // node from p_children, since p_children is holding neighbor nodes.
    // In that case, we're done; n points to the neighbor node in
    // the direction to walk, and we'll return n.
    if (n->depth > m_start_box->depth)
    {
        // If 'n' is outside the start box 'b', something's wrong.
        assert(!Box::is_over_limit(m_start_box, n));

        //get the next-most box whose x or y coord is still within the box b

        while (true)
        {
            // If we've made it here, n shouldn't be a null pointer.
            assert(n != nullptr);

            //when the n->p_children[next] is out of bound of b
            if (n->p_children[m_next] == nullptr ||
                Box::is_over_limit(m_start_box, n->p_children[m_next]))
            {
                break;
            }

            //when the n->p_children[next] is bigger...
            if (n->p_children[m_next]->depth < m_start_box->depth)
            {
                break;
            }

            n = n->p_children[m_next];
        }

        //get the box that touches b
        while (true)
        {
            assert(n != nullptr);
            if (n->p_children[m_cross] == m_start_box)
            {
                break;
            }
            n = n->p_children[m_cross];
        }
    }

    m_neighbor = n;
    return m_neighbor;
}

/*! \brief    Gets the end of the node walk 
 *  \returns  The box node at the end of the walk, or a null reference
 *            if the walk contains no nodes.
 */
Box * BoxIter::end()
{
    if (m_start_box->p_children[m_direction] != nullptr)
    {
        return m_start_box->p_children[m_direction]->p_children[m_prev];
    }

    return nullptr;
}

/*! \brief   Advances to the next neighboring box node
 *  \returns A pointer to the next neighbor, or a null reference
 *           if there are no more nodes in the node walk.
 */
Box * BoxIter::next()
{
    if (!m_neighbor)
    {
        return nullptr;
    }

    m_neighbor = m_neighbor->p_children[m_prev];
    return m_neighbor;
}

// ----------------------------------------------------------------------------
//
//
// Box Implementation
//
//
//
//
// ----------------------------------------------------------------------------

bool Box::is_over_limit(const Box * base, const Box * next_box)
{
    if ((next_box->x > base->x - base->dx / 2 && next_box->x < base->x + base->dx / 2) ||
        (next_box->y > base->y - base->dy / 2 && next_box->y < base->y + base->dy / 2))
    //|| (next_box->z > base->z - base->dz / 2 && next_box->z < base->z + base->dz / 2) )
    {
        return false;
    }
    return true;
}
//
// dir: 0:north, 1 east, 2 south, 3 west
//
void Box::get_neighbors(list<Box *> & neighbors, int dir)
{
    //go through all neighboring boxes in direction i
    BoxIter iter(this, dir);
    Box * neighbor = iter.first();

    if (neighbor == nullptr) return; //no neighbors...

    while (neighbor != iter.end())
    {
        neighbors.push_back(neighbor);
        neighbor = iter.next();
    }
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
 */
bool Box::split(double epsilon)
{

    if (this->dy < epsilon || this->dx < epsilon)
    {
        return false;
    }

    if (!this->is_leaf)
    {
        return false;
    }

    //record the time of this split event, will be used to set priority of children
    //++Box::counter;

    Box * children[4];
    children[0] = new Box(x - dx / 4, y + dy / 4, dx / 2, dy / 2);
    children[1] = new Box(x + dx / 4, y + dy / 4, dx / 2, dy / 2);
    children[2] = new Box(x + dx / 4, y - dy / 4, dx / 2, dy / 2);
    children[3] = new Box(x - dx / 4, y - dy / 4, dx / 2, dy / 2);

    //initialize the new nodes
    for (int i = 0; i < 4; ++i)
    {
        children[i]->depth = this->depth + 1;
        children[i]->dz = this->dz;
        children[i]->z = this->z;
        children[i]->layer = layer;
        children[i]->locator = i;
    }

    //
    // for internal nodes, children [i] is the i-th child
    // for leaves, children [i] is the pointer to first node in i-th adj list
    //
    for (int i = 0; i < 4; ++i)
    {
        //find three other directions
        int prev = (i + 3) % 4;
        int next = (i + 1) % 4;
        int cross = (i + 2) % 4;

        //update easy cases
        children[i]->p_children[i] = p_children[i];
        children[i]->p_children[next] = children[next];
        children[i]->p_children[cross] = children[prev];

        //init box neighbor iterator for direction i
        BoxIter iter(this, i);
        Box * neighbor = iter.first();

        if (!neighbor)
        {
            continue;
        }

        // if neighbor are no smaller
        if (neighbor->depth <= this->depth)
        {
            //after split child 'next' should also point to
            //neighbor in direction i
            children[next]->p_children[i] = neighbor;

            //if neighbor's cross direction point to this, it should
            //instead point to child 'next' after split
            if (neighbor->p_children[cross] == this)
            {
                neighbor->p_children[cross] = children[next];
            }
            continue;
        }

        Box * prev_neighbor = neighbor;

        //indicate if we go across the boundary between child 'i'
        //and 'next' the first time
        bool first_time_cross_between_children = true;

        //if neighbor smaller
        while (neighbor != iter.end())
        {
            //within the strip of child next, neighbor's cross direction
            //should point to next
            if (!is_over_limit(children[next], neighbor))
            {
                neighbor->p_children[cross] = children[next];
            }

            //within the strip of child i, neighbor's cross
            //direction should point to i
            else if (!is_over_limit(children[i], neighbor))
            {
                neighbor->p_children[cross] = children[i];

                //first time cross between child i and next,
                //should update next's i direction pointer
                if (first_time_cross_between_children)
                {
                    first_time_cross_between_children = false;
                    children[next]->p_children[i] = prev_neighbor;
                }
            }
            else
            {
                assert(0);
            }

            prev_neighbor = neighbor;
            neighbor = iter.next();

        } //end while
    }

    //setup/update connections to p_up
    if (this->p_up != nullptr)
    {
        this->refine_p_up();

        if (p_up->is_leaf)
        {
            for (auto & i : children)
            {
                i->p_up = p_up;
                assert(i->p_up->in(i->x, i->y));
            }
        }
        else
        { //p_up is at the same depth as this current box
            for (int i = 0; i < 4; ++i)
            {
                children[i]->p_up = p_up->p_children[i];
                p_up->p_children[i]->p_down = children[i];
                //note: to have complete updates we will have to dig into all the nodes in the
                //subtree of p_up and update their p_down. This can be expensive. We resort
                //to delay this until p_up is needed
            }
        }
    }

    //setup/update connections to p_down
    if (this->p_down != nullptr)
    {
        this->refine_p_down();

        if (p_down->is_leaf)
        {
            for (auto & i : children) i->p_down = p_down;
        }
        else
        {
            //p_down is at the same depth as this current box
            for (int i = 0; i < 4; ++i)
            {
                children[i]->p_down = p_down->p_children[i];
                p_down->p_children[i]->p_up = children[i];
                //note: to have complete updates we will have to dig into all the nodes in the
                //subtree of p_down and update their p_down. This can be expensive. We resort
                //to delay this until p_down is needed
            }
        }
    }

    //setup the children/parent relation
    for (int i = 0; i < 4; ++i)
    {
        this->p_children[i] = children[i];
        this->p_children[i]->p_parent = this;
    }

    this->is_leaf = false;

    return true;
}

/*!
 * Refine p_up (this is called when p_up is at the level lower than the box but is not a leaf)
 */
bool Box::refine_p_up()
{
    return refine_p_up_down(&p_up);
}

/*!
 * Refine p_down (this is called when p_up is at the level lower than the box but is not a leaf)
 */
bool Box::refine_p_down()
{
    return refine_p_up_down(&p_down);
}

bool Box::refine_p_up_down(Box ** ppbox)
{
    if ((*ppbox) == nullptr) return true; //nothing to refine...

    if ((*ppbox)->in(x, y) == false)
    {
        throw GridgenError("internal error: the cell below or above cell " + to_string(id) +
                           " does not contain its center");
    }

    if ((*ppbox)->depth > this->depth)
    {
        throw GridgenError("internal error: the cell below or above cell " + to_string(id) +
                           " is smaller than it");
    }
    if ((*ppbox)->depth == this->depth) return true; //nothing to do
    //(*ppbox)->depth<this->depth
    if ((*ppbox)->is_leaf) return true; //nothing to do

    for (auto & i : (*ppbox)->p_children)
    {
        if (i->in(x, y))
        {
            (*ppbox) = i;
            break;
        }
    }

    return refine_p_up_down(ppbox);
}

/// find a leaf box containing a given point (qx,qy)
Box * Box::find(double qx, double qy)
{
    if (is_leaf && in(qx, qy)) return this;

    assert(!is_leaf);

    for (auto & i : p_children)
    {
        if (i->in(qx, qy)) return i->find(qx, qy);
    }

    throw GridgenError("internal error: no cell contains the point (" + to_string(qx) + ", " +
                       to_string(qy) + ")");
}

void Box::rotate_ctr(double cx, double cy, double rotation, double & rotated_x, double & rotated_y)
{
    rotated_x = cx + (x - cx) * cos(rotation) - (y - cy) * sin(rotation);
    rotated_y = cy + (y - cy) * cos(rotation) + (x - cx) * sin(rotation);
}

} //end namespace cusg
