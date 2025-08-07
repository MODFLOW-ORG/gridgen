/*! \file Cube.cpp

\brief `Cube` and `CubeIter` class implementations

\author <a href="masc.cs.gmu.edu/">MASC group</a>, George Mason University
\author <a href="profile.usgs.gov/langevin/">Christian Langevin</a>, USGS
\bug    No known bugs.

*/

//$Id: $

#include <cassert>
#include <sstream>
#include <stack>
#include "Cube.h"
using namespace std;

namespace cusg
{


// ----------------------------------------------------------------------------
//
//
// Cube Iterator Implementation
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

       when direc=0 (north)
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


/*! \brief Creates a new cube neighbor iterator.
 *
 * Cube neighbor iterator constructor. Runs through QuadTree
 * nodes in the specified direction.
 * \param[in] bb     The starting node
 * \param[in] direc  The direction to travel:
 *                   - 0: Front, North
 *                   - 1: Front, East
 *                   - 2: Front, South
 *                   - 3: Front, West
 * \param[in] layer_in The layer to visit 
 *                     (0 for the front layer, 1 for the back)
 *
 */
CubeIter::CubeIter(const Cube* bb, int direc, int layer_in)
 : b(bb), direction(direc), neighbor(0), layer(layer_in)
{
	prev = (direc + 3) % 4;
	next = (direc + 1) % 4;
	cross = (direc + 2) % 4;
}

//double Cube::r0 = 0;
//int Cube::counter = 0;

/*! \brief Returns the first cube in the node walk.
 *  \returns A pointer to the first node in the direction to be
 *           traveled (as specified when the iterator was instantiated).
 *           If there are no nodes adjacent to the given cube
 *           in that direction, a null reference is returned.
 */           

Cube* CubeIter::First()
{
    // Get the first cube in the node walk.
    //
    // If 'b' is a leaf node, then pChildren points to the neighbors
    // of the node, with 0 = N, 1 = E, 2 = S, 3 = W. So looking up
    // pChildren[direction] will give the adjacent cube in the direction
    // to walk, if such a cube exists.
    //
    // If 'b' is an internal node, then pChildren, as the name implies,
    // holds pointers to the child nodes. For pChildren: 0 = NW, 1 = NE,
    // 2 = SE, 3 = SW. 
	Cube* n;
    n = b->pChildren[direction][layer];

	if (!n)
	{
        // There's no nodes when we walk in this direction,
        // so return a null pointer.
		return 0;
	}

    // If 'b' is a leaf node, its depth will be the same as the
    // node from pChildren, since pChildren is holding neighbor nodes.
    // In that case, we're done; n points to the neighbor node in
    // the direction to walk, and we'll return n.
	if (n->depth > b->depth)
	{
        // If 'n' is outside the start cube 'b', something's wrong.
        assert(!Cube::isOverLimit(b, n));

        //get the next-most cube whose x or y coord is still within the cube b

        // TODO: Rest of function, which deals with internal nodes,
        //       needs more comments. Unfortunately, I am not sure
        //       how internal nodes are being handled.
        while (true)
        {
            // If we've made it here, n shouldn't be a null pointer.
            assert(n != 0);

            //when the n->pChildren[next] is out of bound of b
            if (n->pChildren[next][layer]==NULL || Cube::isOverLimit(b, n->pChildren[next][layer]))
            {
                break;
            }

            //when the n->pChildren[next] is bigger...
            if(n->pChildren[next][layer]->depth < b->depth)
            {
                break;
            }

            n = n->pChildren[next][layer];
        }

        //get the cube that touches b
        while (true)
        {
            assert(n != 0);
            if (n->pChildren[cross][layer] == b)
            {
                break;
            }
            n = n->pChildren[cross][layer];
        }
	}

	neighbor = n;
	return neighbor;
}

/*! \brief    Gets the end of the node walk 
 *  \returns  The cube node at the end of the walk, or a null reference
 *            if the walk contains no nodes.
 */            
Cube* CubeIter::End()
{
	if (b->pChildren[direction][layer]!=NULL)
	{
        return b->pChildren[direction][layer]->pChildren[prev][layer];
	}

	return 0;
}

/*! \brief   Advances to the next neighboring cube node
 *  \returns A pointer to the next neighbor, or a null reference
 *           if there are no more nodes in the node walk.
 */           
Cube* CubeIter::Next()
{
	if (!neighbor)
	{
		return 0;
	}

    neighbor = neighbor->pChildren[prev][layer];
	return neighbor;
}

// ----------------------------------------------------------------------------
//
//
// Cube Implementation
//
//
//
//
// ----------------------------------------------------------------------------

Cube::~Cube() {
    deleteNodes();
}

void Cube::deleteNodes() {
    deleteNodes(this);
}

void Cube::deleteNodes(Cube *c) {
    if (!c || c->isLeaf) {
        return;
    }

    for (short i = 0; i < numChildren; i++) {
        for (short j = 0; j < numLayers; j++) {
            if (c->pChildren[i][j]) {
                deleteNodes(c->pChildren[i][j]);
                delete c->pChildren[i][j];
                c->pChildren[i][j] = NULL;
            }
        }
    }
}

//
// dir: 0:north, 1 east, 2 south, 3 west
// layer: 0:front, 1 back
//
void Cube::getNeighbors(list<Cube*>& neighbors, int dir, int layer)
{
    //go through all neighboring cubes in direction i
    CubeIter iter(this, dir, layer);
    Cube* neighbor = iter.First();

    if (neighbor==NULL) return; //no neighbors...

    while(neighbor != iter.End())
    {
        neighbors.push_back(neighbor);
        neighbor = iter.Next();
    }
}

void Cube::getNorthNeighbors(list<Cube*>& neighbors) { getNeighbors(neighbors,0,0); }
void Cube::getEastNeighbors(list<Cube*>& neighbors)   { getNeighbors(neighbors,1,0); }
void Cube::getSouthNeighbors(list<Cube*>& neighbors)  { getNeighbors(neighbors,2,0); }
void Cube::getWestNeighbors(list<Cube*>& neighbors)   { getNeighbors(neighbors,3,0); }

Cube **Cube::spawn(int back) {
    Cube **children = new Cube*[numChildren];
    double newz = back ? z + dz / 4 : z - dz / 4;

    children[0] = new Cube(x - dx / 4, y + dy / 4, newz, 
        dx / 2, dy / 2, dz / 2);
    children[1] = new Cube(x + dx / 4, y + dy / 4, newz,
        dx / 2, dy / 2, dz / 2);
    children[2] = new Cube(x + dx / 4, y - dy / 4, newz,
        dx / 2, dy / 2, dz / 2);
    children[3] = new Cube(x - dx / 4, y - dy / 4, newz,
        dx / 2, dy / 2, dz / 2);

    return children;

}

/*!  \brief Splits the cube into sub-cubes
 *
 * Splits this cube off into sub-cubes (one for each quadrant
 * of this cube). The new cubes are created as the child nodes of
 * this cube; the cube itself remains unchanged, except that it is no
 * longer considered a leaf node.
 * \param epsilons      If the width, height and depth of this cube are
 *                      all smaller than, or equal to, the corresponding values 
 *                      here, the cube not be split, and \b false is returned.
 * \return \b true, unless:
 *     - the width or height of the cube is smaller than `epsilon`
 *     - this cube is an internal node, so it has already been split
 *     - insufficient memory
 */
bool Cube::split(Vector3d epsilons)
{
    if ((epsilons[0] <= 0 || this->dx < epsilons[0]) && 
        (epsilons[1] <= 0 || this->dy < epsilons[1]) &&
        (epsilons[2] <= 0 || this->dz < epsilons[2]))
    {
        return false;
    }

    if (!this->isLeaf)
    {
        return false;
    }

    //
    // for internal nodes, children [i] is the i-th child
    // for leaves, children [i] is the pointer to first node in i-th adj list
    //
    for (int j = 0; j < numLayers; j++)
    {
        Cube **children = spawn(j);

        for (int i = 0; i < numChildren; ++i)
        {
            if(children[i]==NULL)
            {
                cerr<<"! Error: Cube::split: Not enough memory"<<endl;
                return false;
            }

            children[i]->depth = this->depth + 1;
        }

        for (int i = 0; i < numChildren; ++i)
        {
            //find three other directions
            int prev = (i + 3) % 4;
            int next = (i + 1) % 4;
            int cross = (i + 2) % 4;

            //update easy cases
            children[i]->pChildren[i][j] = pChildren[i][j];
            children[i]->pChildren[next][j] = children[next];
            children[i]->pChildren[cross][j] = children[prev];

            //init cube neighbor iterator for direction i
            CubeIter iter(this, i, j);
            Cube* neighbor = iter.First();

            if (!neighbor)
            {
                continue;
            }

            // if neighbor are no smaller
            if (neighbor->depth <= this->depth)
            {
                //after split child 'next' should also point to
                //neighbor in direction i
                children[next]->pChildren[i][j] = neighbor;

                //if neighbor's cross direction point to this, it should
                //instead point to child 'next' after split
                if (neighbor->pChildren[cross][j] == this)
                {
                    neighbor->pChildren[cross][j] = children[next];
                }
                continue;
            }

            Cube* prevNeighbor = neighbor;

            //indicate if we go across the boundary between child 'i'
            //and 'next' the first time
            bool firstTimeCrossBetweenChildren = true;

            //if neighbor smaller
            while(neighbor != iter.End())
            {
                //within the strip of child next, neighbor's cross direction
                //should point to next
                if (!isOverLimit(children[next], neighbor))
                {
                    neighbor->pChildren[cross][j] = children[next];
                }

                //within the strip of child i, neighbor's cross
                //direction should point to i
                else if (!isOverLimit(children[i], neighbor))
                {
                    neighbor->pChildren[cross][j] = children[i];

                    //first time cross between child i and next,
                    //should update next's i direction pointer
                    if (firstTimeCrossBetweenChildren)
                    {
                        firstTimeCrossBetweenChildren = false;
                        children[next]->pChildren[i][j] = prevNeighbor;
                    }
                }
                else
                {
                    break;
                }

                prevNeighbor = neighbor;
                neighbor = iter.Next();
            }//end while
        }

        //setup the children/parent relation
        for (int i = 0; i < numChildren; ++i)
        {
            this->pChildren[i][j] = children[i];
            this->pChildren[i][j]->pParent = this;
        }
    }        

    this->isLeaf = false;

    return true;
}


/// \brief reverse operation of split
/// This function is not implemented yet. see node.py for detail.
bool Cube::coarsen()
{
    //log a debug message
    //cout<< "Coarsening nodeid: "<<id<<endl;

    //get the nodegroup
    //ng = self.nodegroup

    //ensure node can be coarsened
    if(isLeaf){
        cerr<<"! Error: node is a leaf node so cannot be coarsened"<<endl;
        assert(false);
    }


    //make sure it has four branches that are leaf nodes
    //and then collect neighbors
    for(short x=0;x<numLayers;x++)
    {
        for(short i=0;i<numChildren;i++)
        {
            //find other directions
            int next = (i + 1) % 4;
            int cross = (i + 2) % 4;

            if( pChildren[i][x]->isLeaf==false ){
                cerr<<"! Error: node: Not all children are leaves. Cannot be coarsened."<<endl;
                assert(false);
            }

            //init cube neighbor iterator for direction i
            //for each direction i, there will be two children cubes
            for(short j=0;j<2;j++) //for each of these two cubes
            {
                Cube * kid=(j==0)?pChildren[i][x]:pChildren[next][x];

                //go through all neighboring cubes in direction i
                CubeIter iter(kid, i, x);
                Cube* neighbor = iter.First();

                if (neighbor==NULL) continue; //no neighbors...


                do{

                    //change pointer from kid to this
                    if(neighbor->pChildren[cross][x]==kid)
                        neighbor->pChildren[cross][x]=this;
                    neighbor = iter.Next();
                }
                while(neighbor != iter.End());
            }//end for j
        }

        //update the pChildren link
        for(short i=0;i<numChildren;i++)
        {
            pChildren[i][x]->dirty=true; //mark for deletion
            pChildren[i][x]=pChildren[i][x]->pChildren[i][x]; //now, this pChildren[i] becomes the neighbor at i-direction
        }//end for i
    }

    //reset this node to a leaf
    isLeaf=true;

    return true;
}

string Cube::str() {
    ostringstream os;

    os << dx << 'x' << dy << 'x' << dz << " cube with center (";
    os << x << ", " << y << ", " << z << ")";

    return os.str();
}
}//end namespace cusg
