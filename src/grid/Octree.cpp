/*! \file Octree.cpp
    \author Eric Popelka

    \brief `Octree` class implementation

    \bug No known bugs.
*/

#include <algorithm>
#include <iostream>
#include <stdexcept>
#include "Octree.h"
#include "Util.h"

using namespace std;

namespace cusg {
	Octree::Octree(Cube* root, Vector3d e): pRoot(root), epsilons(e)
	{
        initOctree();
	}

    Octree::Octree(Cube *root, REAL e): pRoot(root)
    {
        epsilons[0] = epsilons[1] = epsilons[2] = e;
        initOctree();
    }

    // TODO: How do we extend this to 3D?
	void Octree::_add_root_nodes()
	{
        nodegroup = NodeGroup(this);

        // Copy the X and Y coordinates for the center of all cells
        // from the parent grid.
        this->X = clone_array(m_mfgrid->X, m_mfgrid->ncol);
        this->Y = clone_array(m_mfgrid->Y, m_mfgrid->nrow);
        // TODO: would we set this->Z here?

        for( int k=0; k< m_mfgrid->nlay; k++)
        {
            for( int i=0; i< m_mfgrid->nrow; i++) 
            {
                double y = Y[i];
                
                for( int j=0; j<m_mfgrid->ncol; j++)
                {
                    double x = X[j];
                    
                    Point3d position;
                    Vector3d dxdydz;
                    
                    //cout<<"HA is3d="<<is3d<<endl;

                    if(is3d)
                    {
                        double z = 0.5 * (m_mfgrid->botm(k,i,j)+m_mfgrid->botm(k+1,i,j));
                        double dz = (m_mfgrid->botm(k,i,j) - m_mfgrid->botm(k+1,i,j));

                        position = Point3d(x, y, z);
                        dxdydz = Vector3d(m_mfgrid->delr[j], m_mfgrid->delc[i], dz);
                    }
                    else
                    {
                        position.set(x, y, 0);
                        dxdydz.set(m_mfgrid->delr[j], m_mfgrid->delc[i], 0);
                    }
                    
                    Box * nodeobj = new Box(position, dxdydz);
        			assert(nodeobj);
        			
        			//TODO: add nodeobj to something?
					//nodegroup.add(nodeobj);
                }//end i
            }//end j
        }//end k
        
        //Create the connections for the root nodes.
        _add_rootnode_connections();
    }

    void Octree::_add_rootnode_connections()
    {
        for( int k=0;k<m_mfgrid->nlay;k++)
        {
            for( int i=0;i<m_mfgrid->nrow;i++)
            {
                for( int j=0;j<m_mfgrid->ncol;j++)
                {
                    int nodeid = m_mfgrid->get_nodeid(k, i, j);
                    
                    //cout<<"Adding connections for nodeid: "<<nodeid<<endl;
                    /* TODO: nodegroup is all Box nodes. We have Cubes... 
                    Cube * nodeobj = nodegroup[nodeid];

                    if( i < m_mfgrid->nrow - 1 )
                    {
                        //make connection to front
                        int tonode = m_mfgrid->get_nodeid(k, i + 1, j);
                        int fldir = -2;
                        nodeobj->addConnection(nodegroup[tonode], fldir);
                    }
                    
                    if( i > 0 )
                    {
                        //make connection to back
                        int tonode = m_mfgrid->get_nodeid(k, i - 1, j);
                        int fldir = 2;
                        nodeobj->addConnection(nodegroup[tonode], fldir);
                    }
                    
                    if( j < m_mfgrid->ncol - 1 )
                    {
                        //make connection to right
                        int tonode = m_mfgrid->get_nodeid(k, i, j + 1);
                        int fldir = 1;
                        nodeobj->addConnection(nodegroup[tonode], fldir);
                    }
                    
                    if( j > 0 )
                    {
                        //make connection to left
                        int tonode = m_mfgrid->get_nodeid(k, i, j - 1);
                        int fldir = -1;
                        nodeobj->addConnection(nodegroup[tonode], fldir);
                    }
                    
                    if( k < m_mfgrid->nlay - 1 )
                    {
                        //make connection down
                        int tonode = m_mfgrid->get_nodeid(k + 1, i, j);
                        int fldir = -3;
                        nodeobj->addConnection(nodegroup[tonode], fldir);
                    }

                    if( k > 0 )
                    {
                        //make connection up
                        int tonode = m_mfgrid->get_nodeid(k - 1, i, j);
                        int fldir = 3;
                        nodeobj->addConnection(nodegroup[tonode], fldir);
                    }
                    */
                }//end j
            }//end i
        }//end k
	}

    void Octree::balancePhase() {
	    list<Cube*> leaves;
	    pRoot->getLeaves(leaves);

	    //we use a heap here to retrieve the deepest leaves
	    vector< pair<int,Cube*> > cube_pq;
	    for(list<Cube*>::iterator i=leaves.begin();i!=leaves.end();i++)
	    {
	        Cube * cube=*i;
	        pair<int,Cube*> tmp(cube->depth,cube);
	        cube_pq.push_back(tmp);
	        push_heap(cube_pq.begin(),cube_pq.end());
	    }


	    //loop until the heap is empty
	    while(cube_pq.empty()==false)
	    {

	        pair<int,Cube*> tmp=cube_pq.front();
	        Cube * cube=tmp.second;
	        pop_heap(cube_pq.begin(),cube_pq.end());
	        cube_pq.pop_back();

	        //check the neighbors
            for(int j=0;j<Cube::numLayers;j++) {
                for(int i=0;i<Cube::numChildren;i++){

                    Cube * nei=cube->pChildren[i][j]; //for leaves, neighbors are stored in pChildren

                    if(nei==NULL) continue;

                    //visiting all neighbors
                    if(nei->depth < cube->depth-1)
                    {
                        bool results=nei->split(epsilons); //ask neighbor to split

                        if(results) //the neighbor did split
                        {
                            for(int k=0;k<Cube::numChildren;k++){ //enqueue neighbors' kid

                                Cube * nei_kid=nei->pChildren[k][j];

                                nei_kid->updateStatus();
                                //if(nei_kid->status==Cube::IN)
                                {
                                    PQ->push(nei_kid);
                                }

                                pair<int,Cube*> tmp(nei_kid->depth,nei_kid);
                                cube_pq.push_back(tmp);
                                push_heap(cube_pq.begin(),cube_pq.end());
                            }//end for k
                        }
                        else{
                            cout<<"! Warning: Octree::balancePhase: split failed"<<endl;
                        }
                    }
                }//end for i
            }// end for j
	    }//end while
    } // end balancePhase

	int *Octree::createnodenumber_array()
	{
		int nsize = nodegroup.size();
		int * nodenumber_array=new int[nsize];
		assert(nodenumber_array);
		memset(nodenumber_array,-1,sizeof(int)*nsize);
		
		int leaf_size=nodes();

		for( int nodenumber=0;  nodenumber < leaf_size; nodenumber++)
		{
			int nodeid=nodeid_array[nodenumber];
			nodenumber_array[nodeid] = nodenumber;
		}
		
		return nodenumber_array;
	}

    bool Octree::expand()
    {
        while(!PQ->empty()) {
            Cube* b = PQ->extract();
			//b might not be a leaf since it could already be split in expand(Cube* b), and PQ is not updated there
			if (b->isLeaf && b->split(epsilons))
			{

			    // check the status of the cube here
			    //assert(b->status == Cube::IN || b->status == Cube::ON);

                for (int j = 0; j < 2; ++j)
                {
                    for (int i = 0; i < 4; ++i)
                    {
                        b->pChildren[i][j]->updateStatus();
                        insertNode(b->pChildren[i][j]);
                    }			
                }
				return true;
			}			
		}
		return false;
    }

	Cube* Octree::getCube(Cube* root, Point3d &q)
	{
		if (q[0] > root->x + root->dx / 2 || q[0] < root->x - root->dx / 2
			|| q[1] > root->y + root->dy / 2 || q[1] < root->y - root->dy / 2
            || q[2] > root->z + root->dz / 2 || q[2] < root->z - root->dz / 2)
		{
			return 0;
		}

		Cube* b = root;
		while (!b->isLeaf)
		{
			double dx = q[0] - b->x;
			double dy = q[1] - b->y;
            double dz = q[2] - b->z;
			if (dx <= 0 && dy >= 0 && dz <= 0)
			{
				b = b->pChildren[0][0];
			}
			else if (dx >= 0 && dy >= 0 && dz <= 0)
			{
				b = b->pChildren[1][0];
			}
			else if (dx >= 0 && dy <= 0 && dz <= 0)
			{
				b = b->pChildren[2][0];
			}
			else if (dx <= 0 && dy <= 0 && dz <= 0)
			{
				b = b->pChildren[3][0];
			}
			if (dx <= 0 && dy >= 0 && dz >= 0)
			{
				b = b->pChildren[0][1];
			}
			else if (dx >= 0 && dy >= 0 && dz >= 0)
			{
				b = b->pChildren[1][1];
			}
			else if (dx >= 0 && dy <= 0 && dz >= 0)
			{
				b = b->pChildren[2][1];
			}
			else if (dx <= 0 && dy <= 0 && dz >= 0)
			{
				b = b->pChildren[3][1];
			}
		}
		return b;
	}

    list< pair<int,int> > Octree::get_node_connections(int nodenumber)
    {
        int nodeid = nodeid_array[nodenumber];
        list< pair<int,int> > connections;
        /* TODO: nodegroup is all Box nodes. We have Cubes...
        Cube * cube = nodegroup[nodeid];
        list<Cube*> nei[8]; //neighbors from 8 sides
        for(short layer = 0; layer < 2; layer++) {
            for(short i=0;i<4;i++) cube->getNeighbors(nei[i],i,layer);
        }

        assert(cube->isLeaf);



        //8 sides
        typedef list<Cube*>::iterator LCIT;
        for(LCIT i=nei[0].begin();i!=nei[0].end();i++) connections.push_back(make_pair((*i)->id, 2));
        for(LCIT i=nei[3].begin();i!=nei[3].end();i++) connections.push_back(make_pair((*i)->id,-1));
        for(LCIT i=nei[1].begin();i!=nei[1].end();i++) connections.push_back(make_pair((*i)->id, 1));
        for(LCIT i=nei[2].begin();i!=nei[2].end();i++) connections.push_back(make_pair((*i)->id,-2));

        // TODO: not sure what to do here to handle 3D cases?
        */
        return connections;
    }

    int * Octree::getnodeid_array()
    {
    	if(nodeid_array==NULL || nodelay_array==NULL || rebuild_nodeid_array)
    	{
    		nodegroup.rebuild_nodeid_array=true;
    		nodeid_array=nodegroup.get_nodeid_array();
    		nodelay_array=nodegroup.get_nodelay_array();
    		rebuild_nodeid_array=false;
    		rebuild_nodenumber_array=true;
    	}

    	return nodeid_array;
    }

    int * Octree::getnodenumber_array()
    {
    	if(nodenumber_array==NULL || rebuild_nodenumber_array)
    	{
	    	nodenumber_array = createnodenumber_array();
	    	rebuild_nodenumber_array=false;
    	}
    	
    	return nodenumber_array;
	}

    void Octree::initOctree()
    {
		PQ = new SequentialCubeQueue();
		assert(PQ);

		pRoot->updateStatus();

		insertNode(pRoot);
        is3d = true;

        //add the root nodes and their connections.  This creates a nodegroup
        //object and assigns it to the grid as an attribute.  Nodegroup is 
        //an array that contains all of the node objects.
        _add_root_nodes();

        //assign a variable that indicates whether or not nodeid_array needs to
        //be rebuilt because nodes have been added or removed.
        rebuild_nodeid_array = true;
        rebuild_nodenumber_array=true;
        nodeid_array=getnodeid_array();
		nodenumber_array=getnodenumber_array();

    }

    void Octree::insertNode(Cube* c)
	{
		PQ->push(c);
	}

    void Octree::subdividePhase()
    {
	    int ct = 0;
	    while(PQ->empty()==false)
	    {
	        expand();
	        ct++;
	    }
    }

    int Octree::treeHeight()
    {
        list<Cube *> l;
        pRoot->getLeaves(l);
        Cube *c = *max_element(l.begin(), l.end(), Cube::comparePtrDepths);
        return c->depth;
    }
    
} // end namespace cusg
