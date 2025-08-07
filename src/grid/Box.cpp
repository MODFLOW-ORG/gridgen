/*! \file Box.cpp

\brief `Box` and `BoxIter` class implementations

\author <a href="masc.cs.gmu.edu/">MASC group</a>, George Mason University
\author <a href="profile.usgs.gov/langevin/">Christian Langevin</a>, USGS
\bug    No known bugs.

*/

//$Id: $

#include "Box.h"
#include "ModflowGrid.h"
#include <assert.h>

namespace cusg
{
	int generateID()
	{
		static int id=1;
		return id++;
	}

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


/*! \brief Creates a new box neighbor iterator.
 *
 * Box neighbor iterator constructor. Runs through QuadTree
 * nodes in the specified direction.
 * \param[in] bb     The starting node
 * \param[in] direc  The direction to travel:
 *                   - 0: North
 *                   - 1: East
 *                   - 2: South
 *                   - 3: West
 */
BoxIter::BoxIter(const Box* bb, int direc):b(bb), direction(direc), neighbor(0)
{
	prev = (direc + 3) % 4;
	next = (direc + 1) % 4;
	cross = (direc + 2) % 4;
}

//double Box::r0 = 0;
//int Box::counter = 0;

/*! \brief Returns the first box in the node walk.
 *  \returns A pointer to the first node in the direction to be
 *           traveled (as specified when the iterator was instantiated).
 *           If there are no nodes adjacent to the given box
 *           in that direction, a null reference is returned.
 */           

Box* BoxIter::First()
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
	Box* n = b->pChildren[direction];
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
        // If 'n' is outside the start box 'b', something's wrong.
        assert(!Box::isOverLimit(b, n));

        //get the next-most box whose x or y coord is still within the box b

        while (true)
        {
            // If we've made it here, n shouldn't be a null pointer.
            assert(n != 0);

            //when the n->pChildren[next] is out of bound of b
            if ( n->pChildren[next]==NULL || Box::isOverLimit(b, n->pChildren[next]))
            {
                break;
            }

            //when the n->pChildren[next] is bigger...
            if(n->pChildren[next]->depth < b->depth)
            {
                break;
            }

            n = n->pChildren[next];
        }

        //get the box that touches b
        while (true)
        {
            assert(n != 0);
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
Box* BoxIter::End()
{
	if (b->pChildren[direction]!=NULL)
	{
		return b->pChildren[direction]->pChildren[prev];
	}

	return 0;
}

/*! \brief   Advances to the next neighboring box node
 *  \returns A pointer to the next neighbor, or a null reference
 *           if there are no more nodes in the node walk.
 */           
Box* BoxIter::Next()
{
	if (!neighbor)
	{
		return 0;
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

bool Box::isOverLimit(const Box* base, const Box* nextBox)
{
    if ((nextBox->x > base->x - base->dx / 2 && nextBox->x < base->x + base->dx / 2)
        || (nextBox->y > base->y - base->dy / 2 && nextBox->y < base->y + base->dy / 2) )
        //|| (nextBox->z > base->z - base->dz / 2 && nextBox->z < base->z + base->dz / 2) )
    {
        return false;
    }
    return true;
}
//
// dir: 0:north, 1 east, 2 south, 3 west
//
void Box::getNeighbors(list<Box*>& neighbors, int dir)
{
    //go through all neighboring boxes in direction i
    BoxIter iter(this, dir);
    Box* neighbor = iter.First();

    if (neighbor==NULL) return; //no neighbors...

    while(neighbor != iter.End())
    {
        neighbors.push_back(neighbor);
        neighbor = iter.Next();
    }
}

void Box::getNorthNeighbors(list<Box*>& neighbors) { getNeighbors(neighbors,0); }
void Box::getEastNeighbors(list<Box*>& neighbors)   { getNeighbors(neighbors,1); }
void Box::getSouthNeighbors(list<Box*>& neighbors)  { getNeighbors(neighbors,2); }
void Box::getWestNeighbors(list<Box*>& neighbors)   { getNeighbors(neighbors,3); }

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
bool Box::split(double epsilon)
{

    //cout<< "refining nodeid: "<<id<<endl;


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

    Box* children[4];
    children[0] = new Box(x - dx / 4, y + dy / 4, dx / 2, dy / 2);
    children[1] = new Box(x + dx / 4, y + dy / 4, dx / 2, dy / 2);
    children[2] = new Box(x + dx / 4, y - dy / 4, dx / 2, dy / 2);
    children[3] = new Box(x - dx / 4, y - dy / 4, dx / 2, dy / 2);

    //initialize the new nodes
    for (int i = 0; i < 4; ++i)
    {
        if(children[i]==NULL)
        {
            cerr<<"! Error: Box::split: Not enough memory"<<endl;
            return false;
        }

        children[i]->depth = this->depth + 1;
        children[i]->dz=this->dz;
        children[i]->z=this->z;
        children[i]->layer=layer;
        children[i]->locator=i;
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
        Box* neighbor = iter.First();

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

        Box* prevNeighbor = neighbor;

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

        }//end while
    }

    //setup/update connections to pUp
    if(this->pUp!=NULL)
    {
		this->refine_pUp();

        if(pUp->isLeaf)
        {
			for (int i = 0; i < 4; ++i){ 
				children[i]->pUp=pUp; 
				assert(children[i]->pUp->in(children[i]->x,children[i]->y)); 
			}
        }
        else
        {   //pUp is at the same depth as this current box
            for (int i = 0; i < 4; ++i)
			{
                children[i]->pUp=pUp->pChildren[i];
                pUp->pChildren[i]->pDown=children[i];
				//note: to have complete updates we will have to dig into all the nodes in the
                //subtree of pUp and update their pDown. This can be expensive. We resort
                //to delay this until pUp is needed
            }
        }
    }

    //setup/update connections to pDown
    if(this->pDown!=NULL)
    {
		this->refine_pDown();

        if(pDown->isLeaf)
        {
            for (int i = 0; i < 4; ++i) children[i]->pDown=pDown;
        }
        else
        {
			//pDown is at the same depth as this current box
            for (int i = 0; i < 4; ++i){
                children[i]->pDown=pDown->pChildren[i];
                pDown->pChildren[i]->pUp=children[i];
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

/// \brief reverse operation of split
bool Box::coarsen()
{
    //log a debug message
    //cout<< "Coarsening nodeid: "<<id<<endl;

    //ensure node can be coarsened
    if(isLeaf){
        cerr<<"! Error: node["<<id<<"] is a leaf node so cannot be coarsened"<<endl;
        assert(false);
    }


    //make sure it has four branches that are leaf nodes
    //and then collect neighbors
    
	//Box * children[4]; //backup

    for(short i=0;i<4;i++)
    {
        //find other directions
        int next = (i + 1) % 4;
        int cross = (i + 2) % 4;

        if( pChildren[i]->isLeaf==false ){
            cerr<<"! Error: node["<<id<<"]: Not all children are leaves. Cannot be coarsened."<<endl;
            assert(false);
        }

        //init box neighbor iterator for direction i
        //for each direction i, there will be two children boxes
        for(short j=0;j<2;j++) //for each of these two boxes
        {
            Box * kid=(j==0)?pChildren[i]:pChildren[next];

            //go through all neighboring boxes in direction i
            BoxIter iter(kid, i);
            Box* neighbor = iter.First();

            if (neighbor==NULL) continue; //no neighbors...


            do{

                //change pointer from kid to this
                if(neighbor->pChildren[cross]==kid)
                    neighbor->pChildren[cross]=this;
                neighbor = iter.Next();
            }
            while(neighbor != iter.End());
        }//end for j
    }

    //update the pChildren link
    for(short i=0;i<4;i++)
    {
        pChildren[i]->dirty=true; //mark for deletion
        pChildren[i]=pChildren[i]->pChildren[i]; //now, this pChildren[i] becomes the neighbor at i-direction
    }//end for i


    //reset this node to a leaf
    isLeaf=true;

    //update up and down pointers
    for(short i=0;i<4;i++)
    {
        if(pChildren[i]->pUp!=NULL) if(pChildren[i]->pUp->pDown==pChildren[i]) pChildren[i]->pUp->pDown=this;
        if(pChildren[i]->pDown!=NULL) if(pChildren[i]->pDown->pUp==pChildren[i]) pChildren[i]->pDown->pUp=this;
    }//end for i

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
    if((*ppbox)==NULL) return true; //nothing to refine...

	if((*ppbox)->in(x,y)==false)
	{
        cout<<"! Error: (*ppbox) refine_pUpDown error: box["<<id<<"] with depth="<<depth
            <<" has (*ppbox)=box["<<(*ppbox)->id<<"] with depth="<<(*ppbox)->depth<<endl;
        assert(false);
	}

    if((*ppbox)->depth>this->depth){
        cout<<"! Error: (*ppbox) depth error: box["<<id<<"] with depth="<<depth
            <<" has (*ppbox)=box["<<(*ppbox)->id<<"] with depth="<<(*ppbox)->depth<<endl;
        assert(false);
    }
    if((*ppbox)->depth==this->depth) return true; //nothing to do
    //(*ppbox)->depth<this->depth
    if((*ppbox)->isLeaf) return true; //nothing to do

    for(int i=0;i<4;i++)
    {
        if( (*ppbox)->pChildren[i]->in(x,y) )
        {
            (*ppbox)=(*ppbox)->pChildren[i];
            break;
        }
    }

    return refine_pUpDown(ppbox);
}

/*!
 * compute the corners (vertices) of this box from its center and dimensions
 * UL: upper-left vertex, UR: upper-right vertex
 * LR: lower-right vertex, LL: lower-left vertex
 */
void Box::getCorners(Point2d& UL, Point2d& UR, Point2d& LR, Point2d& LL ) const
{
    double w2=dx/2;
    double h2=dy/2;
    double x_l=x-w2;
    double x_r=x+w2;
    double y_l=y-h2;
    double y_u=y+h2;

    UL.set(x_l, y_u);
    LL.set(x_l, y_l);
    UR.set(x_r, y_u);
    LR.set(x_r, y_l);
}




/// find a leaf box containing a given point (qx,qy)
Box * Box::find(double qx, double qy)
{
	if(isLeaf && in(qx,qy) ) return this;

	//if(childmfgrid != NULL)
	//{
	//	Point2d pt(qx, qy);
	//	return childmfgrid->find_nodeobj(pt);
	//}


	assert(!isLeaf);

	for(int i=0;i<4;i++){
		if(pChildren[i]->in(qx,qy))
			return pChildren[i]->find(qx,qy);
	}

	cerr<<"! Error: Box::find: Locate a box error!"<<endl;
	return NULL;
}

//rotate
void Box::rotate(double cx, double cy, double rotation, double& ldX, double& ldY, double& luX, double& luY, double& rdX, double& rdY, double& ruX, double& ruY)
{
	ldX = cx + (x - dx / 2  - cx) * cos(rotation) - (y - dy / 2 - cy) * sin(rotation);
	ldY = cy + (y - dy / 2 - cy) * cos(rotation) + (x - dx / 2 - cx) * sin(rotation);

	luX = cx + (x - dx / 2  - cx) * cos(rotation) - (y + dy / 2 - cy) * sin(rotation);
	luY = cy + (y + dy / 2 - cy) * cos(rotation) + (x - dx / 2 - cx) * sin(rotation);

	rdX = cx + (x + dx / 2  - cx) * cos(rotation) - (y - dy / 2 - cy) * sin(rotation);
	rdY = cy + (y - dy / 2 - cy) * cos(rotation) + (x + dx / 2 - cx) * sin(rotation);

	ruX = cx + (x + dx / 2  - cx) * cos(rotation) - (y + dy / 2 - cy) * sin(rotation);
	ruY = cy + (y + dy / 2 - cy) * cos(rotation) + (x + dx / 2 - cx) * sin(rotation);
	
}
void Box::rotateCtr(double cx, double cy, double rotation,double& resx, double& resy)
{
	resx= cx +  (x - cx) * cos(rotation) - (y - cy)* sin(rotation);
	resy = cy + (y - cy)*cos(rotation) + (x - cx)*sin(rotation);
}

ostream& operator<<(ostream& out, const Box& box)
{
    Point2d UL, UR, LR, LL;
    box.getCorners(UL, UR, LR, LL);
    out<<"Box ["<<box.id<<"] = ("<<LL<<") * ("<<UR<<") number="<<box.number<<" depth="<<box.depth;
    return out;
}

}//end namespace cusg
