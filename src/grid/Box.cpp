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
       then, 1st (returned by First() is the cell, on the far right
       and end is the cell on the far left and is also the cell stored in pChildren[0]

*/

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
    : b(starting_node), direction(direction), neighbor(nullptr)
{
    prev = (direction + 3) % 4;
    next = (direction + 1) % 4;
    cross = (direction + 2) % 4;
}

/*! \brief Returns the first box in the node walk.
 *  \returns A pointer to the first node in the direction to be
 *           traveled (as specified when the iterator was instantiated).
 *           If there are no nodes adjacent to the given box
 *           in that direction, a null reference is returned.
 */

Box * BoxIter::First()
{
    // Get the first box in the node walk.
    //
    // If 'b' is a leaf node, then pChildren points to the neighbors
    // of the node, with 0 = N, 1 = E, 2 = S, 3 = W. So looking up
    // pChildren[direction] will give the adjacent box in the direction
    // to walk, if such a box exists.
    //
    // If 'b' is an internal node, then pChildren, as the name implies,
    // holds pointers to the child nodes. For pChildren: 0 = NW, 1 = NE,
    // 2 = SE, 3 = SW.
    Box * n = b->pChildren[direction];
    if (!n)
    {
        // There's no nodes when we walk in this direction,
        // so return a null pointer.
        return nullptr;
    }

    // If 'b' is a leaf node, its depth will be the same as the
    // node from pChildren, since pChildren is holding neighbor nodes.
    // In that case, we're done; n points to the neighbor node in
    // the direction to walk, and we'll return n.
    if (n->depth > b->depth)
    {
        // If 'n' is outside the start box 'b', something's wrong.
        assert(!Box::isOverLimit(b, n));

        //get the next-most box whose x or y coord is still within the box b

        while (true)
        {
            // If we've made it here, n shouldn't be a null pointer.
            assert(n != nullptr);

            //when the n->pChildren[next] is out of bound of b
            if (n->pChildren[next] == nullptr || Box::isOverLimit(b, n->pChildren[next]))
            {
                break;
            }

            //when the n->pChildren[next] is bigger...
            if (n->pChildren[next]->depth < b->depth)
            {
                break;
            }

            n = n->pChildren[next];
        }

        //get the box that touches b
        while (true)
        {
            assert(n != nullptr);
            if (n->pChildren[cross] == b)
            {
                break;
            }
            n = n->pChildren[cross];
        }
    }

    neighbor = n;
    return neighbor;
}

/*! \brief    Gets the end of the node walk 
 *  \returns  The box node at the end of the walk, or a null reference
 *            if the walk contains no nodes.
 */
Box * BoxIter::End()
{
    if (b->pChildren[direction] != nullptr)
    {
        return b->pChildren[direction]->pChildren[prev];
    }

    return nullptr;
}

/*! \brief   Advances to the next neighboring box node
 *  \returns A pointer to the next neighbor, or a null reference
 *           if there are no more nodes in the node walk.
 */
Box * BoxIter::Next()
{
    if (!neighbor)
    {
        return nullptr;
    }

    neighbor = neighbor->pChildren[prev];
    return neighbor;
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

bool Box::isOverLimit(const Box * base, const Box * nextBox)
{
    if ((nextBox->x > base->x - base->dx / 2 && nextBox->x < base->x + base->dx / 2) ||
        (nextBox->y > base->y - base->dy / 2 && nextBox->y < base->y + base->dy / 2))
    //|| (nextBox->z > base->z - base->dz / 2 && nextBox->z < base->z + base->dz / 2) )
    {
        return false;
    }
    return true;
}
//
// dir: 0:north, 1 east, 2 south, 3 west
//
void Box::getNeighbors(list<Box *> & neighbors, int dir)
{
    //go through all neighboring boxes in direction i
    BoxIter iter(this, dir);
    Box * neighbor = iter.First();

    if (neighbor == nullptr) return; //no neighbors...

    while (neighbor != iter.End())
    {
        neighbors.push_back(neighbor);
        neighbor = iter.Next();
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

    if (!this->isLeaf)
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
        children[i]->pChildren[i] = pChildren[i];
        children[i]->pChildren[next] = children[next];
        children[i]->pChildren[cross] = children[prev];

        //init box neighbor iterator for direction i
        BoxIter iter(this, i);
        Box * neighbor = iter.First();

        if (!neighbor)
        {
            continue;
        }

        // if neighbor are no smaller
        if (neighbor->depth <= this->depth)
        {
            //after split child 'next' should also point to
            //neighbor in direction i
            children[next]->pChildren[i] = neighbor;

            //if neighbor's cross direction point to this, it should
            //instead point to child 'next' after split
            if (neighbor->pChildren[cross] == this)
            {
                neighbor->pChildren[cross] = children[next];
            }
            continue;
        }

        Box * prevNeighbor = neighbor;

        //indicate if we go across the boundary between child 'i'
        //and 'next' the first time
        bool firstTimeCrossBetweenChildren = true;

        //if neighbor smaller
        while (neighbor != iter.End())
        {
            //within the strip of child next, neighbor's cross direction
            //should point to next
            if (!isOverLimit(children[next], neighbor))
            {
                neighbor->pChildren[cross] = children[next];
            }

            //within the strip of child i, neighbor's cross
            //direction should point to i
            else if (!isOverLimit(children[i], neighbor))
            {
                neighbor->pChildren[cross] = children[i];

                //first time cross between child i and next,
                //should update next's i direction pointer
                if (firstTimeCrossBetweenChildren)
                {
                    firstTimeCrossBetweenChildren = false;
                    children[next]->pChildren[i] = prevNeighbor;
                }
            }
            else
            {
                assert(0);
            }

            prevNeighbor = neighbor;
            neighbor = iter.Next();

        } //end while
    }

    //setup/update connections to pUp
    if (this->pUp != nullptr)
    {
        this->refine_pUp();

        if (pUp->isLeaf)
        {
            for (auto & i : children)
            {
                i->pUp = pUp;
                assert(i->pUp->in(i->x, i->y));
            }
        }
        else
        { //pUp is at the same depth as this current box
            for (int i = 0; i < 4; ++i)
            {
                children[i]->pUp = pUp->pChildren[i];
                pUp->pChildren[i]->pDown = children[i];
                //note: to have complete updates we will have to dig into all the nodes in the
                //subtree of pUp and update their pDown. This can be expensive. We resort
                //to delay this until pUp is needed
            }
        }
    }

    //setup/update connections to pDown
    if (this->pDown != nullptr)
    {
        this->refine_pDown();

        if (pDown->isLeaf)
        {
            for (auto & i : children) i->pDown = pDown;
        }
        else
        {
            //pDown is at the same depth as this current box
            for (int i = 0; i < 4; ++i)
            {
                children[i]->pDown = pDown->pChildren[i];
                pDown->pChildren[i]->pUp = children[i];
                //note: to have complete updates we will have to dig into all the nodes in the
                //subtree of pDown and update their pDown. This can be expensive. We resort
                //to delay this until pDown is needed
            }
        }
    }

    //setup the children/parent relation
    for (int i = 0; i < 4; ++i)
    {
        this->pChildren[i] = children[i];
        this->pChildren[i]->pParent = this;
    }

    this->isLeaf = false;

    return true;
}

/*!
 * Refine pUp (this is called when pUp is at the level lower than the box but is not a leaf)
 */
bool Box::refine_pUp()
{
    return refine_pUpDown(&pUp);
}

/*!
 * Refine pDown (this is called when pUp is at the level lower than the box but is not a leaf)
 */
bool Box::refine_pDown()
{
    return refine_pUpDown(&pDown);
}

bool Box::refine_pUpDown(Box ** ppbox)
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
    if ((*ppbox)->isLeaf) return true; //nothing to do

    for (auto & i : (*ppbox)->pChildren)
    {
        if (i->in(x, y))
        {
            (*ppbox) = i;
            break;
        }
    }

    return refine_pUpDown(ppbox);
}

/// find a leaf box containing a given point (qx,qy)
Box * Box::find(double qx, double qy)
{
    if (isLeaf && in(qx, qy)) return this;

    assert(!isLeaf);

    for (auto & i : pChildren)
    {
        if (i->in(qx, qy)) return i->find(qx, qy);
    }

    throw GridgenError("internal error: no cell contains the point (" + to_string(qx) + ", " +
                       to_string(qy) + ")");
}

void Box::rotateCtr(double cx, double cy, double rotation, double & rotated_x, double & rotated_y)
{
    rotated_x = cx + (x - cx) * cos(rotation) - (y - cy) * sin(rotation);
    rotated_y = cy + (y - cy) * cos(rotation) + (x - cx) * sin(rotation);
}

} //end namespace cusg
