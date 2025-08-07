#include "GhostNode.h"
#include <fstream>
using namespace  std;

namespace cusg
{
	//check if c is between [a, b]
	bool between(double a, double b, double c)
	{
		if((a < c && c <b) || (a > c && c >b))
			return true;
		else 
			return false;
	}
	void computeGN(Box* m, Box* n, int dir, vector<GhostNode>& gns)
	{//decide the node j first

		GhostNode gn;
		gn.mNum = m->number;	gn.nNum = n->number;

		//if it is in north or south
		if(dir == 0 || dir == 2)
		{
			list<Box*> nl;
			if(between(n->x - n->dx, n->x, m->x))
			{//search the 3 direction(west)
				n->getNeighbors(nl, 3);
			}
			else if(between(n->x, n->x + n->dx, m->x))
			{//search for the 1 direction
				n->getNeighbors(nl, 1);
			}

			if(nl.empty()) return;

			if(nl.size() > 2)
			{
				cout<<"! Error : the node "<<n->id<<" has more than two ("<<nl.size()<<") neighbor box nodes in dir="<<dir<<endl;
				exit(0);
			}

			if(nl.size()==1){
				gn.j1 = gn.j2 = nl.front()->number;
				gn.apha1 = ( (n->x - m->x) / (n->x - nl.front()->x))  /  2;// fabs((n->x - m->x)/n->dx);
				gn.apha2 = gn.apha1;
			}
			else{
				assert(nl.size() == 2);
				gn.j1 = nl.front()->number;
				gn.j2 = nl.back()->number;
				//they have the same X coordinates
				gn.apha1 = ( (n->x - m->x) / (n->x - nl.front()->x)) / 2;//fabs( ( (n->x - m->x)/n->dx ) / 2 );
				gn.apha2 = gn.apha1;
			}
		}
		else if (dir == 1 || dir == 3)
		{
			list<Box*> nl;
			if(between(n->y + n->dy, n->y, m->y))
				n->getNeighbors(nl, 0);
			else if(between(n->y, n->y - n->dy, m->y))
				n->getNeighbors(nl, 2);

			if(nl.empty()) return;

			if(nl.size() > 2)
			{
				cout<<"! Error : the node "<<n->id<<" has more than two ("<<nl.size()<<") neighbor box nodes in dir="<<dir<<endl;
				exit(0);
			}

			if(nl.size() == 1){
				gn.j1 = gn.j2 = nl.front()->number;
				gn.apha1 = gn.apha2 =  ( (n->y - m->y) / (n->y - nl.front()->y)) / 2;//fabs( (n->y - m->y) / n->dy);
			}
			else 
			{
				assert(nl.size() == 2);
				gn.j1 = nl.front()->number;
				gn.j2 = nl.back()->number;
				gn.apha1 = gn.apha2 = ( (n->y - m->y) / (n->y - nl.front()->y)) / 2;//fabs( ( (n->y - m->y) / n->dy) / 2 );
			}
		}
		else
			assert(false);

		if(gn.j1==-1 || gn.j2==-1) return;//don't add it to gns if j1 or j2 is inactive...

		gns.push_back(gn);
	}
	void computeGNNeigh(Box* m, list<Box*>& nlist, int dir, vector<GhostNode>& gns)
	{
		//neighbors
		for(list<Box*>::iterator lit = nlist.begin(); lit != nlist.end(); ++lit)
		{
			Box* tb = *lit;
			
			if(tb->number == -1)
				continue;

			if(tb->depth >= m->depth)
				continue;

			//the node should be considered
			computeGN(m, tb, dir, gns);
		}
	}

	void computeGNC(QuadTree3D* qtree, vector<GhostNode>& gns)
	{
		//loop over all the leaf node
		for(NodeGroup::iterator nit = qtree->nodegroup.begin(); nit != qtree->nodegroup.end(); ++nit)
		{
			Box* b = *nit;
			if(!b->isLeaf || b->number == -1)
				continue;
			
			list<Box*> nei[4];
			for(int i = 0; i < 4; i++)
			{
				b->getNeighbors(nei[i], i);
			}

			for (int j = 0; j < 4; j++)
			{
				computeGNNeigh(b, nei[j], j, gns);
			}
		}

	}

	void writeGhostNode(string filePath, vector<GhostNode>& gns, int offset=1)
	{
		ofstream ofile;
		ofile.open(filePath.c_str(),std::ofstream::out | std::ofstream::trunc); 
		for(vector<GhostNode>::iterator git = gns.begin(); git != gns.end(); ++git)
		{
			GhostNode& ghn = *git;
			ofile<<ghn.nNum+offset<<"\t"<<ghn.mNum+offset<<"\t"<<ghn.j1+offset<<"\t"<<ghn.j2+offset<<"\t"<<ghn.apha1<<"\t"<<ghn.apha2<<endl;
		}
		ofile.close();

	}
}