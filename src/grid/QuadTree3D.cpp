/*! \file QuadTree3D.cpp

\brief QuadTree 3D grid

\author <a href="masc.cs.gmu.edu/">MASC group</a>, George Mason University,
<a href="profile.usgs.gov/langevin/">Christian Langevin</a>, USGS
\bug    No known bugs.

*/

//   $Id: $

#include "QuadTree3D.h"
#include "GridIntersection.h"
#include "ascii_grid.h"
#include <algorithm>

#include "vert_pass.h"
#include "BoxPlyClip.h"
using namespace std;

namespace cusg
{


	/*!
	Initialize the quadtree grid with a parent grid.  The parent grid should
	be a ModflowGrid object.
	*/

	QuadTree3D::QuadTree3D(ModflowGrid * mfgrid, int * chunksize)
		:nodeid_array(NULL), nodenumber_array(NULL), nodelay_array(NULL)
	{
		cout<<"- A QuadTree object was instantiated."<<endl;

		// Store the modflowgrid, too 
		m_mfgrid=mfgrid;

		//Set the chunksize, which is used to allocate node memory.
		if(chunksize ==NULL){
			this->m_chunksize = new int(m_mfgrid->nodes * 0.2);
		}
		else{
			this->m_chunksize = chunksize;
		}

		//determine if this is a 2d or 3d grid
		is3d=m_mfgrid->is3d;

		//add the root nodes and their connections.  This creates a nodegroup
		//object and assigns it to the grid as an attribute.  Nodegroup is 
		//an array that contains all of the node objects.
		_add_root_nodes();

		//assign a variable that indiactes whether or not nodeid_array needs to
		//be rebuilt because nodes have been added or removed.
		rebuild_nodeid_array = true;
		rebuild_nodenumber_array=true;
		one_based_node_numbering=true;
		nodeid_array=getnodeid_array();
		nodenumber_array=getnodenumber_array();

		//how z values of top and bot should be determined
		m_bot_layer_z_operator=m_top_layer_z_operator=vector<LAYER_Z_VALUE_OPEARTOR>(m_mfgrid->nlay,Z_VALUE_REPLICATE_OPEARTOR); //default is replicate
	
		vertical_pass_through = false;
	}

	QuadTree3D::~QuadTree3D(void)
	{

	}

	void getSpecialPointIntersection(QuadTree3D* qtree, const Point2d& pt, int max_level, int at_layer, list<Box*> & result)
	{
		int i, j, x_id, y_id, id;
		ModflowGrid* mf = qtree->getModflowGrid();
		Index2d id2d = mf->getIndex(pt);

		if(id2d[0]==INT_MAX || id2d[1]==INT_MAX|| id2d[0] ==INT_MIN || id2d[1]==INT_MIN) return ;

		for(i = -1; i <=1; i++)
		{
			y_id = id2d[0] + i;
			if(y_id < 0 || y_id >=mf->nrow)
				continue;

			for(j = -1; j <= 1; j++)
			{
				x_id = id2d[1] + j;
				if(x_id < 0 || x_id >= mf->ncol)
					continue;

				id = mf->get_nodeid(at_layer, y_id, x_id);
				Box* box = qtree->nodegroup[id];
				if(box->in(pt[0], pt[1]))
				{
					box = box->find(pt[0], pt[1]);
					assert(box);
					result.push_back(box);
				}
			}
		}
	}


	void QuadTree3D::refine(Point2d& pt, int max_level, int at_layer)
	{
		list<Box*> bl;
		getSpecialPointIntersection(this, pt, max_level, at_layer, bl);

		while(!bl.empty())
		{
			Box* b = bl.front();
			bl.pop_front();

			if(b->depth >= max_level) continue;//too deep
			bool r = b->split(0);
			assert(r);//make sure that split is successful

			for(int c = 0; c < 4; c++)
			{
				nodegroup.add(b->pChildren[c]);
				if(b->pChildren[c]->in(pt[0], pt[1]))
					bl.push_back(b->pChildren[c]);
			}

			this->rebuild_nodeid_array=true;
		}

	}

	/**************************************************************/
	//the old version
	//void QuadTree3D::refine( Point2d& pt, int max_level, int at_layer)
	//{
	//	//refined=false;
	//	bool refined=true;
	//	do
	//	{
	//		refined=false;
	//		GridIntersection gi(max_level,at_layer);
	//		bool r=gi.intersect((QuadTree3D*)this,pt);
	//		if(r)
	//		{
	//			assert(gi.getIntersections().size()==1);
	//			Box * b=gi.getIntersections().front(); //there should be only one

	//			//post-check to see if the point is on the corner or on the boundary

	//			if(b->depth>=max_level) continue; //too deep
	//			bool r=b->split(0);
	//			assert(r); //make sure that split is successful

	//			for(int c=0;c<4;c++)
	//				nodegroup.add(b->pChildren[c]);

	//			refined=true;
	//			this->rebuild_nodeid_array=true;
	//		}
	//	}
	//	while(refined);
	//}
	
	//bool lineIntersectBox(LineSeg2d line, Box* b)
	//{
	//	double a[2], b[2];
	//	double c[2]={line.pt[0][0], line.pt[0][1]};
	//	double d[2]={line.pt[1][0], line.pt[1][1]};
	//	double p[2];
	//	return true;
	//}

	//given the  intersection precondition provided by GridIntersection
	bool isReallyIntersect(Box* b, const LineSeg2d& line)
	{
		const double* s = line.pt[0].get();
		const double* e = line.pt[1].get();
		double p[2];
		//create each edge of b
		Point2d corners[4];
		corners[0].set(b->x-b->dx/2,b->y+b->dy/2); //upper left corner
		corners[1].set(b->x+b->dx/2,b->y+b->dy/2); //upper right corner
		corners[2].set(b->x+b->dx/2,b->y-b->dy/2); //lower right corner
		corners[3].set(b->x-b->dx/2,b->y-b->dy/2); //lower left corner

		bool b_intersect = false;
		int id1, id2;
		//if collinear
		for(int i = 0; i < 4; i++)
		{
			id1=i; id2 = (i+1)%4;
			char code = SegSegInt(corners[id1].get(), corners[id2].get(), s, e, p);
			if(code == '0') continue;
			else if(code == '1')
			{
				b_intersect = true; break;
			}
			else if(code == 'v')
			{
				//what vertex is it
				if(Between_strict(s, e, corners[id1].get()) || Between_strict(s, e, corners[id2].get()))
				{
					b_intersect = true; break;
				}
			}
			else if(code == 'e')
			{
				//if any strictly between/inside
				if(Between_strict(s, e, corners[id1].get()) || Between_strict(s, e, corners[id2].get()) ||
					Between_strict(corners[id1].get(), corners[id2].get(), s) || Between_strict(corners[id1].get(), corners[id2].get(), e))
				{
					b_intersect = true; break;
				}
			}
		}

		if(b)
			return true;

		//last chance, both point are in box
		if(b->in(s[0], s[1]) && b->in(e[0], e[1]))
			return true;

		return false;


	}
	/***********************************************************************************************************/
	//grid and line seg intersection for Quadtree3D refinement
#if 1
	bool specialLineIntersect(QuadTree3D * qtree, const LineSeg2d& line, int max_level, int at_layer, vector<Box*>& results)
	{
		GridIntersection gi(max_level, at_layer);
		LineSeg2d clipped_line;
		bool rclip=gi.clip_line(qtree,line,clipped_line);
		if(rclip==false) return true; //line segment is entirely outside, so there is no intersection
		if(getSqLength(clipped_line) < 1e-10)
			return true;

		gi.getIntersections().clear();
		
		list<Box*> bl1, bl2;
		getSpecialPointIntersection(qtree, line.pt[0], max_level, at_layer, bl1);
		getSpecialPointIntersection(qtree, line.pt[1], max_level, at_layer, bl2);
		
		//now, visit b2 neighbors recursively until we reaches b1
		set<Box*> bs; //, cur_res;
		bs.insert(bl1.begin(), bl1.end());
		bs.insert(bl2.begin(), bl2.end());

		for(set<Box*>::iterator sit = bs.begin(); sit != bs.end(); ++sit)
		{
			if(isReallyIntersect(*sit, line))
			{
				results.push_back(*sit);
				//cur_res.insert(*sit);
			}
			else
			{
				cerr<<"Warning: line not intersect with box actually!\n";
			}
		}

#if 0
		list<Box*> open(bs.begin(), bs.end());
		set<Box*> visited(bs.begin(), bs.end());
#else
		//list<Box*> open(cur_res.begin(), cur_res.end());
		//set<Box*> visited;//(cur_res.begin(), cur_res.end());

		list<Box*> open(results.begin(), results.end());
		set<Box*> visited(results.begin(), results.end());
#endif

		while(open.empty()==false) //not empty
		{
			Box * b=open.front();
			open.pop_front();

			//get the neighbors of the box b
			list<Box*> neighbors[4];

			for (int i = 0; i < 4; i++)
			{
				b->getNeighbors(neighbors[i], i);
			}

			//process neighbors
			for(int i=0;i<4;i++)
			{
				for(list<Box*>::iterator n=neighbors[ i ].begin();n!=neighbors[ i ].end();n++)
				{
					Box * nb=*n;

					if(visited.find(nb)!=visited.end()) continue; //visited

					if( gi.intersect(nb,line) )
					{
						if(isReallyIntersect(nb, line))
						{
							results.push_back(nb);
							open.push_back(nb);
						}
					}
					//visited before
					visited.insert(nb);
				}
			}

		}//end while

		if(results.empty() == false)
			return true;
		else
			return  false;
	}
#endif

	//void QuadTree3D::refine(LineSeg2d line, int max_level, int at_layer)
	//{
	//	bool refined = true;
	//	do
	//	{
	//		refined = false;
	//		GridIntersection gi(max_level, at_layer);
	//		bool r = gi.intersect((QuadTree3D*)this, line);
	//		//vector<Box*> results;
	//		//bool r = specialLineIntersect(this, line, max_level, at_layer, results);
	//		vector<Box*>& results = gi.getIntersections();
	//		
	//		if(r)
	//		{
	//			//debugging
	//			//vector<Box*> boxes;
	//			//specialLineIntersect(this, line, max_level, at_layer, results);

	//			for(vector<Box*>::iterator i=results.begin();i!=results.end();i++)
	//			{
	//				Box * b=*i;

	//				if(b->depth>=max_level) continue; //too deep

	//				bool r=b->split(0);
	//				assert(r); //make sure that split is successful

	//				for(int c=0;c<4;c++)
	//					nodegroup.add(b->pChildren[c]);

	//				refined=true;
	//				this->rebuild_nodeid_array=true;
	//			}
	//		}

	//	}while(refined);

	//}

	
	void QuadTree3D::refine(LineSeg2d line, int max_level, int at_layer)
	{
		bool refined = true;
		do
		{
			refined = false;
			//GridIntersection gi(max_level, at_layer);
			//bool r = gi.intersect((QuadTree3D*)this, line);
			vector<Box*> results;
			bool r = specialLineIntersect(this, line, max_level, at_layer, results);
			//vector<Box*>& results = gi.getIntersections();
			
			if(r)
			{
				//debugging
				//vector<Box*> boxes;
				//specialLineIntersect(this, line, max_level, at_layer, results);

				for(vector<Box*>::iterator i=results.begin();i!=results.end();i++)
				{
					Box * b=*i;

					if(b->depth>=max_level) continue; //too deep

					bool r=b->split(0);
					assert(r); //make sure that split is successful

					for(int c=0;c<4;c++)
						nodegroup.add(b->pChildren[c]);

					refined=true;
					this->rebuild_nodeid_array=true;
				}
			}

		}while(refined);

	}

	//void QuadTree3D::refine(LineSeg2d line, int max_level, int at_layer)
	//{
	//	//const LineSeg2d& line=*_line;

	//	bool refined=true;

	//	do
	//	{

	//		refined=false;
	//		GridIntersection gi(max_level,at_layer);
	//		bool r=gi.intersect((QuadTree3D*)this,line);

	//		if(r)
	//		{
	//			vector<Box*>& boxes=gi.getIntersections();
	//			for(vector<Box*>::iterator i=boxes.begin();i!=boxes.end();i++)
	//			{
	//				Box * b=*i;

	//				if(b->depth>=max_level) continue; //too deep

	//				//check if really intersect
	//				//todo:




	//				bool r=b->split(0);
	//				assert(r); //make sure that split is successful

	//				for(int c=0;c<4;c++)
	//					nodegroup.add(b->pChildren[c]);

	//				refined=true;
	//				this->rebuild_nodeid_array=true;
	//			}
	//		}

	//	}
	//	while(refined);
	//}

	

	void QuadTree3D::refine( c_plyline& arc, int max_level, int at_layer)
	{

		bool refined=true;

		do
		{

			refined=false;

			vector<Box*> results;

			GridIntersection gi(max_level,at_layer);

			ply_vertex * ptr=arc.getHead();

			bool r=false;
			do
			{
				ply_vertex * next=ptr->getNext();
				if(next!=NULL)
				{
					LineSeg2d line(ptr->getPos(),next->getPos());
					
#if 1
					bool myr = specialLineIntersect(this, line, max_level, at_layer, results);
					if(myr)
						r = true;
#else
					bool myr=gi.intersect((QuadTree3D*)this,line);
					if(myr)
					{
						r=true;
						vector<Box*>& cur_inter = gi.getIntersections();
						results.insert(results.end(), cur_inter.begin(), cur_inter.end());
					}
#endif
				}
				ptr=next;

				//break;
			}
			while(ptr!=NULL);

			if(r)
			{
				for(vector<Box*>::iterator i=results.begin();i!=results.end();i++)
				{
					Box * b=*i;
					if(b->depth>=max_level) continue; //too deep

					bool splitted=b->split(0);
					if(splitted)
					{
						for(int c=0;c<4;c++)
						{
							nodegroup.add(b->pChildren[c]);
						}
						refined=true;
						this->rebuild_nodeid_array=true;
					}
				}//end for i
			}

		}
		while(refined);

	}

	//void QuadTree3D::refine( c_plyline& arc, int max_level, int at_layer)
	//{

	//	bool refined=true;

	//	do
	//	{

	//		refined=false;
	//		GridIntersection gi(max_level,at_layer);

	//		ply_vertex * ptr=arc.getHead();

	//		bool r=false;
	//		do
	//		{
	//			ply_vertex * next=ptr->getNext();
	//			if(next!=NULL)
	//			{
	//				LineSeg2d line(ptr->getPos(),next->getPos());
	//				bool myr=gi.intersect((QuadTree3D*)this,line);
	//				if(myr) r=true;
	//			}
	//			ptr=next;
	//		}
	//		while(ptr!=NULL);

	//		if(r)
	//		{
	//			vector<Box*>& boxes=gi.getIntersections();
	//			for(vector<Box*>::iterator i=boxes.begin();i!=boxes.end();i++)
	//			{
	//				Box * b=*i;
	//				if(b->depth>=max_level) continue; //too deep

	//				bool r=b->split(0);
	//				if(r)
	//				{
	//					for(int c=0;c<4;c++) nodegroup.add(b->pChildren[c]);
	//					refined=true;
	//					this->rebuild_nodeid_array=true;
	//				}
	//			}//end for i
	//		}

	//	}
	//	while(refined);
	//}

	
	bool isReallyIntersect(BoxPlyClip* bpc, Box* b)
	{
		double xmin = min(b->x - b->dx/2, b->x + b->dx/2);
		double xmax = max(b->x - b->dx/2, b->x + b->dx/2);
		double ymin = min(b->y - b->dy/2, b->y + b->dy/2);
		double ymax = max(b->y - b->dy/2, b->y + b->dy/2);

		list<BoxPlyClip*> cres;
		bpc->clip(xmin, xmax, ymin, ymax, cres);

		double area(0.0);
		for (list<BoxPlyClip*>::iterator cit = cres.begin(); cit != cres.end(); ++cit)
		{
			GIS_polygon nplygon;
			convertClipToPolygon(*cit, &nplygon);
			area += nplygon.getArea();

			delete (*cit);
		}

		//delete bpc;
		if(area >= 1e-7)
			return true;
		else
			return false;
	}

	bool isReallyIntersect(c_ply& ply, Box* b)
	{
		/*****************************************************************/
		BoxPlyClip* bpc  = convertPlyToClip(ply);
		bool r = isReallyIntersect(bpc, b);
		
		bpc->destroy();
		delete bpc;
		return r;
	}
	bool isReallyIntersect(c_polygon& pgon, Box* b)
	{
		BoxPlyClip* bpc  = convertSimplePolygonToClip(pgon);
		
		bool r = isReallyIntersect(bpc, b);
		
	bpc->destroy();
	delete bpc;
		return r;
	}

	void QuadTree3D::refine( c_ply& ply, int max_level, int at_layer)
	{
		bool refined=true;

		do
		{

			refined=false;
			GridIntersection gi(max_level,at_layer);
			bool r=gi.intersect((QuadTree3D*)this,ply);

			if(r)
			{
				vector<Box*>& boxes=gi.getIntersections();

				//cout<<"add box size="<<boxes.size()*4<<endl;

				for(vector<Box*>::iterator i=boxes.begin();i!=boxes.end();i++)
				{
					Box * b=*i;
					if(b->depth>=max_level) continue; //too deep

					if(!isReallyIntersect(ply, b)) continue; //no really intersect or contain

					bool r=b->split(0);

#if DEBUG
					assert(r); //make sure that split is successful
#endif

					for(int c=0;c<4;c++)
						nodegroup.add(b->pChildren[c]);

					refined=true;
					this->rebuild_nodeid_array=true;
				}
			}

		}
		while(refined);

	}


	//refine cells enclosed by pgon
	void QuadTree3D::refine( c_polygon& pgon, int max_level, int at_layer)
	{
		bool refined=true;

		do
		{
			refined=false;
			GridIntersection gi(max_level,at_layer);
			bool r=gi.intersect((QuadTree3D*)this,pgon);

			if(r)
			{
				vector<Box*>& boxes=gi.getIntersections();

				for(vector<Box*>::iterator i=boxes.begin();i!=boxes.end();i++)
				{
					Box * b=*i;
					if(b->depth>=max_level) continue; //too deep

					if(!isReallyIntersect(pgon, b))
					{
						cout<<"warning: QuadTree3D::refine( c_polygon& pgon, int max_level, int at_layer) no intersection!\n";
						continue;//not really intersect
					}

					bool r=b->split(0);

#if DEBUG
					assert(r); //make sure that split is successful
#endif

					for(int c=0;c<4;c++)
						nodegroup.add(b->pChildren[c]);

					refined=true;
					this->rebuild_nodeid_array=true;
				}//end for i (intersected boxes)
			}

		}
		while(refined);
	}

	///return true if this is the final refinement
	bool QuadTree3D::smooth_refinement
		(Box * box, vector< pair<int,Box*> >& box_pq, int level_tol_horizontal, int level_tol_vertical)
	{
		list< pair<int,int> >  connectionlist = get_node_connections(box);

		bool refine_once=true;

		//check the neighbors
		// go through each connection
		typedef list< pair<int,int> >::iterator IT;
		for(IT ic=connectionlist.begin();ic!=connectionlist.end();ic++)
		{
			int tonodeid=ic->first;
			int dir=ic->second;
			Box * nei= nodegroup[tonodeid];

			if(nei==NULL) continue;
			assert(nei->isLeaf);

			bool refine_nei=false;

			if( abs(dir)==3 ) //vertical neighbor
			{
				if(nei->depth < box->depth-level_tol_vertical) refine_nei=true;
				if(nei->depth+1 < box->depth-level_tol_vertical)  refine_once=false;
			}
			else //horizontal neighbor
			{
				if(nei->depth < box->depth-level_tol_horizontal) refine_nei=true;
				if(nei->depth+1 < box->depth-level_tol_horizontal) refine_once=false;
			}

			//visiting all neighbors
			if(refine_nei)
			{
				bool results=nei->split(0);
				if(results) //the neighbor did split
				{
					for(int k=0;k<4;k++) //enqueue neighbors' kid
					{
						Box * nei_kid=nei->pChildren[k];

						//add to node group
						nodegroup.add(nei_kid);

						//add to the priority queue
						pair<int,Box*> tmp(nei_kid->depth,nei_kid);
						box_pq.push_back(tmp);
						push_heap(box_pq.begin(),box_pq.end());

						//flag to rebuild
						this->rebuild_nodeid_array=true;
					}//end for k
				}
				else
				{
					cerr<<"! Warning: QuadTree3D::smooth_refinement: split failed"<<endl;
				}
			}
		}//end for i

		return refine_once;
	}

	/// Smooth the refinement so no connections exceeds the level tolerance (leveltol).
	/// level_tol_horizontal is horizontal level tolerance
	/// level_tol_vertical is vertical level tolerance
	void QuadTree3D::smooth_refinement(int level_tol_horizontal, int level_tol_vertical)
	{
		//get all leaves
		list<Box*> leaves;
		int totalnodesize=nodegroup.size();
		for(int l=0;l<totalnodesize;l++)
		{
			Box * box=nodegroup[l];
			if(box->isLeaf) leaves.push_back(box);
		}//end for l

		//we use a heap here to retrieve the deepest leaves
		vector< pair<int,Box*> > box_pq;
		for(list<Box*>::iterator i=leaves.begin();i!=leaves.end();i++)
		{
			Box * box=*i;
			pair<int,Box*> tmp(box->depth,box);
			box_pq.push_back(tmp);
			push_heap(box_pq.begin(),box_pq.end());
		}

		//loop until the heap is empty
		while(box_pq.empty()==false)
		{
			pair<int,Box*> tmp=box_pq.front();
			Box * box=tmp.second;
			pop_heap(box_pq.begin(),box_pq.end());
			box_pq.pop_back();

			if(box->isLeaf==false) continue; //not a leaf anymore, already split...

			while(smooth_refinement(box, box_pq, level_tol_horizontal, level_tol_vertical)==false);

		}//end while

#if DEBUG
		//check the correctness of the result
		totalnodesize=nodegroup.size();
		for(int l=0;l<totalnodesize;l++)
		{
			Box * box=nodegroup[l];

			if(box->isLeaf)
			{
				list< pair<int,int> >  connectionlist = get_node_connections(box);


				//check the neighbors
				// go through each connection
				typedef list< pair<int,int> >::iterator IT;
				for(IT ic=connectionlist.begin();ic!=connectionlist.end();ic++)
				{
					int tonodeid=ic->first;
					int dir=ic->second;
					Box * nei= nodegroup[tonodeid];

					if(nei==NULL) continue;
					assert(nei->isLeaf);

					bool refine_nei=false;

					if( abs(dir)==3 ) //vertical neighbor
					{
						if(nei->depth < box->depth-level_tol_vertical) refine_nei=true;
					}
					else //horizontal neighbor
					{
						if(nei->depth < box->depth-level_tol_horizontal) refine_nei=true;
					}

					assert(refine_nei==false);
				}
			}
		}//end for l
#endif //DEBUG

	}

	/// Return an array of size (nodes) that contains the area for all cells.
	double * QuadTree3D::get_cell_areas()
	{
		//log the creation
		//cout<<"Getting cell areas."<<endl;
		int nodesize=nodes();

		double * area = new double[nodesize];
		assert(area);

		for(int nodenumber=0; nodenumber<nodesize; nodenumber++)
		{
			int nodeid = nodeid_array[nodenumber];
			Box * nodeobj = nodegroup[nodeid];
			area[nodenumber] = nodeobj->dx * nodeobj->dy;
		}

		return area;
	}



	/*!
	\brief Return the four vertices for the specified nodeid.
	*/    
	/*
	vector<Point2d> get_vertices(int nodeid)
	{

	}
	*/

	/*! 
	\brief Return a tuple containing two points (the lower left and the upper 
	right) in global grid coordinates.
	*/    	
	pair<Point2d,Point2d> QuadTree3D::get_extent()
	{
		return m_mfgrid->get_extent();
	}

	/*
	///Return the nodeid for the modflow grid.
	int get_nodeid(int i, int j)
	{
	//return i * ncol + j;
	}
	*/

	/*! 
	Return an array of size nodes that contains the nodeid for all
	active leaf nodes.  Thus the array contains:

	nodeid_array[nodenumber] = nodeid

	Note that this array is rebuilt, but only if the
	nodegroup.rebuild_nodeid_array flag is True.  This flag is set to
	true anytime quadtree nodes are added or removed.
	*/  
	int * QuadTree3D::getnodeid_array()
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

	///get nodelay_array
	int * QuadTree3D::get_nodelay()
	{
		if(nodeid_array==NULL || nodelay_array==NULL || rebuild_nodeid_array)
		{
			nodegroup.rebuild_nodeid_array=true;
			nodeid_array=nodegroup.get_nodeid_array();
			nodelay_array=nodegroup.get_nodelay_array();
			rebuild_nodeid_array=false;
			rebuild_nodenumber_array=true;
		}

		return nodelay_array;
	}

	/*!
	\brief Create a nodenumber_array of size nodegroup with the following:

	node number is the id of the node in CSR or Coo array

	nodenumber_array[nodeid] = nodenumber
	*/	
	int * QuadTree3D::createnodenumber_array()
	{
		int nsize = nodegroup.size();
		int * nodenumber_array=new int[nsize];
		assert(nodenumber_array);
		memset(nodenumber_array,(int)-INT_MAX,sizeof(int)*nsize);

		int leaf_size=nodes();

		for( int nodenumber=0;  nodenumber < leaf_size; nodenumber++)
		{
			int nodeid=nodeid_array[nodenumber];
			nodenumber_array[nodeid] = nodenumber;
		}

		return nodenumber_array;
	}

	/*!
	Create an array of size modflowgrid.nlay that has the number of
	nodes in each layer.
	*/
	void QuadTree3D::create_nodelay_array()
	{
		//this is not implemented in grid.py
		cerr<<"! ERROR: create_nodelay_array is not implemented"<<endl;
		return;
	}

	/*! 
	Return an array of len(nodegroup) that contains the nodenumber for all
	active leaf nodes.  Thus the array contains:

	nodenumber_array[nodeid] = nodenumber

	Note that this array is rebuilt, but only if the
	nodegroup.rebuild_nodeid_array flag is True.  This flag is set to
	true anytime quadtree nodes are added or removed.
	*/
	int * QuadTree3D::getnodenumber_array()
	{
		if(nodenumber_array==NULL || rebuild_nodenumber_array)
		{
			nodenumber_array = createnodenumber_array();
			rebuild_nodenumber_array=false;
		}

		return nodenumber_array;
	}

	///get box from id
	Box * QuadTree3D::get_nodeobj(int nodenumber)
	{
		int nodeid = nodeid_array[nodenumber];
		return nodegroup[nodeid];
	}

	void QuadTree3D::getnode(void * node, int nodenumber)
	{
		Box * box=(Box*)node;
		int nodeid = nodeid_array[nodenumber];
		*box=*(nodegroup[nodeid]);
	}

	///get number of nodes
	int QuadTree3D::get_nodes()
	{
		return nodegroup.node_size();
	}

	int QuadTree3D::get_nodenumber(int nodeid)
	{
		int nodenumber = nodenumber_array[nodeid];
		if( nodenumber == -1)
		{
			cerr<<"! ERROR: This nodeid "<<nodeid<<" does not presently have a nodenumber."<<endl;
			assert(false);
		}

		return nodenumber;
	}

	
	list<pair<int, int> > QuadTree3D::get_node_connections(int nodenumber)
	{
		int nodeid = nodeid_array[nodenumber];
		Box * box = nodegroup[nodeid];
		assert(box->isLeaf);
		return get_node_connections(box);
	}

	list<pair<int, int> > QuadTree3D::get_node_connections(Box* box)
	{
		//neighbors from 4 sides and up and down
		list<Box*> nei[6];

		//get each of the 4 neighbors
		for(short i=0;i<4;i++) box->getNeighbors(nei[i],i);

		//get all leaves of box->pUp
		if( box->pUp!=NULL )
		{
			box->refine_pUp();
			box->pUp->getLeaves(nei[4]);
		}

		//get all leaves of box->pDown
		if( box->pDown!=NULL )
		{
			box->refine_pDown();
			box->pDown->getLeaves(nei[5]);
		}

		//4 sides and up and down
		list< pair<int,int> > connections;
		for(list<Box*>::iterator i=nei[0].begin();i!=nei[0].end();i++) if((*i)->active) connections.push_back(make_pair((*i)->id, 2));
		for(list<Box*>::iterator i=nei[3].begin();i!=nei[3].end();i++) if((*i)->active) connections.push_back(make_pair((*i)->id,-1));
		for(list<Box*>::iterator i=nei[1].begin();i!=nei[1].end();i++) if((*i)->active) connections.push_back(make_pair((*i)->id, 1));
		for(list<Box*>::iterator i=nei[2].begin();i!=nei[2].end();i++) if((*i)->active) connections.push_back(make_pair((*i)->id,-2));
		for(list<Box*>::iterator i=nei[4].begin();i!=nei[4].end();i++) if((*i)->active) connections.push_back(make_pair((*i)->id,3));

		for(list<Box*>::iterator i=nei[5].begin();i!=nei[5].end();i++) if((*i)->active) connections.push_back(make_pair((*i)->id,-3));

		return connections;
	}



	///get box connections from id
	///Return a sorted list of connections, where a connection is
	///a list [tonodeid, fldir].  fldir is -1, +1, -2, +2 for -x, +x, -y, +y, -z, +z.
	list< pair<int,int> > QuadTree3D::get_node_connections_vp(int nodenumber, vector<float>* area)
	{
		int nodeid = nodeid_array[nodenumber];
		Box * box = nodegroup[nodeid];
		assert(box->isLeaf);
		return get_node_connections_vp(box, area);
	}

	list< pair<int,int> > QuadTree3D::get_node_connections_vp(Box * box, vector<float>* area)
	{
		//neighbors from 4 sides and up and down
		list<Box*> nei[6];

		//get each of the 4 neighbors
		for(short i=0;i<4;i++) box->getNeighbors(nei[i],i);

		//get all leaves of box->pUp
		if( box->pUp!=NULL )
		{
			box->refine_pUp();
			box->pUp->getLeaves(nei[4]);
		}

		//get all leaves of box->pDown
		if( box->pDown!=NULL )
		{
			box->refine_pDown();
			box->pDown->getLeaves(nei[5]);
		}

		//4 sides and up and down
		list< pair<int,int> > connections;
		for(list<Box*>::iterator i=nei[0].begin();i!=nei[0].end();i++) if((*i)->active) connections.push_back(make_pair((*i)->id, 2));
		for(list<Box*>::iterator i=nei[3].begin();i!=nei[3].end();i++) if((*i)->active) connections.push_back(make_pair((*i)->id,-1));
		for(list<Box*>::iterator i=nei[1].begin();i!=nei[1].end();i++) if((*i)->active) connections.push_back(make_pair((*i)->id, 1));
		for(list<Box*>::iterator i=nei[2].begin();i!=nei[2].end();i++) if((*i)->active) connections.push_back(make_pair((*i)->id,-2));
		for(list<Box*>::iterator i=nei[4].begin();i!=nei[4].end();i++) if((*i)->active) connections.push_back(make_pair((*i)->id,3));

		//vertical pass here
		if(area != NULL)
		{
			int csize = connections.size();
			for(int i = 0; i < csize; i++)
				area->push_back(-1.0f);
		}
		//if(addVerticalPass)
		//{
			//up direction
			vector<VPTPair> vpb1;
			getVPT(box, vpb1, false);
			for(vector<VPTPair>::iterator vpit = vpb1.begin(); vpit != vpb1.end(); ++vpit)
			{
				//cout<<"vertical pass"<<endl;
				assert(vpit->bSrc == box);
				connections.push_back(make_pair(vpit->bDst->id, 3));
				if(area != NULL) area->push_back((float)(vpit->area));
			}
		//}

		for(list<Box*>::iterator i=nei[5].begin();i!=nei[5].end();i++) if((*i)->active) connections.push_back(make_pair((*i)->id,-3));

		if(area != NULL)
		{
			int addSize = connections.size() - area->size();
			for(int i = 0; i < addSize; i++) area->push_back(-1.0f);
		}
		//add the down vertical pass
			vector<VPTPair> vpb2;
			getVPT(box, vpb2,true);
			for(vector<VPTPair>::iterator vpit = vpb2.begin(); vpit != vpb2.end(); ++vpit)
			{
				//cout<<"vertical pass"<<endl;
				assert(vpit->bSrc == box);
				connections.push_back(make_pair(vpit->bDst->id, -3));
				if(area != NULL) area->push_back((float)(vpit->area));				
			}

		return connections;
	}


	/// direction indiactor of size nja (+-1 for x, +-2 for y +-3 for z)
	int QuadTree3D::boundayNodes(vector<int>& nodeids, int dir)
	{
		//for(NodeGroup::reverse_iterator i=nodegroup.rbegin();i!=nodegroup.rend();i++)
		int size=this->nodes();
		for(int i=0;i<size;i++)
		{
			Box * box=this->get_nodeobj(i);
			list<Box*> neighbors;
			if(dir==-1){ //-x boundary
				box->getWestNeighbors(neighbors);
			}
			else if(dir==1){ //+x boundary
				box->getEastNeighbors(neighbors);
			}
			else if(dir==-2){ //-y boundary
				box->getSouthNeighbors(neighbors);
			}
			else if(dir==2){ //-y boundary
				box->getNorthNeighbors(neighbors);
			}

			if(neighbors.empty()) nodeids.push_back(box->id);
		}

		return nodeids.size();
	}


	/*! \brief Obtain CSRData

	Return a UsgCsrData object, which contains all of the necessary
	information for creating an unstructured grid model.

	If top[nodes] and bot[nodes] are provided they are used in the
	calculation of fahl.

	*/
	CSRData QuadTree3D::get_usg_csr_data()
	{
		cout<<"- Creating the UsgCsrData object"<<endl;

		if(this->rebuild_nodeid_array)
		{
			this->getnodeid_array();
		}

		if(this->rebuild_nodenumber_array)
		{
			this->getnodenumber_array();
		}

		//this is the number of leaves of the quadtree
		int box_size=this->nodes();

		//go through each cell and count diagonal position and connections
		int nja = 0;
		for( int nid=0;nid<box_size;nid++)  //nid in range(self.nodes) )
		{
			nja++;

			vector<float> areas;
			int csize= 0;
			
			if(vertical_pass_through)
				csize = get_node_connections_vp(nid, &areas).size();
			else
				csize = get_node_connections(nid).size();
			nja += csize;
		}

		//allocate the CSR arrays.
		//cout<<"Allocating ia of size "<<box_size+1<<" and ja of size "<<nja<<endl;

		int * ia = new int[box_size+1];
		int * ja = new int[nja];
		double * area = get_cell_areas();
		float * cl1 = new float[nja];
		float * cl2 = new float[nja];
		float * fahl = new float[nja];
		int * fldr = new int[nja];
		int * iac = new int[box_size]; // number of connection +1 for each cell, array of size  nodes

		assert(ia && ja && cl1 && cl2 && fahl && fldr && iac);

		//Go through each connection and fill the arrays.
		int japos = 0;
		int id_offset=(one_based_node_numbering)?1:0;
		for( int nid=0;nid<box_size;nid++)
		{
			//get an object representation of this nodeid
			//nodeobj = self.get_nodeobj(nodenumber)

			Box * nodeobj = get_nodeobj(nid);
			double dxdydzn[3] = {nodeobj->dx, nodeobj->dy, nodeobj->dz};

			if(nodeobj->dx<0 || nodeobj->dy<0 || nodeobj->dz<0)
			{
				cout<<"("<<nodeobj->dx<<","<<nodeobj->dy<<","<<nodeobj->dz<<")"<<endl;
			}

			//diagonal position
			ia[nid] = japos+id_offset;
			ja[japos] = nid+id_offset;
			cl1[japos] = 0.0;
			cl2[japos] = 0.0;
			fahl[japos] = 0.0;
			fldr[japos] = 0;
			japos++;

			vector<float> conAreas;
			list< pair<int,int> >  connectionlist;
			if(vertical_pass_through)
				connectionlist = get_node_connections_vp(nid, &conAreas);
			else
				connectionlist = get_node_connections(nid);

			iac[nid]=connectionlist.size()+1;

			if(sortconnections)
			{
				if(vertical_pass_through)
					sortareabyconid(connectionlist, conAreas);//sort the area first
				sortconnectionbynodenumber(connectionlist);
			}

			// go through each connection
			typedef list< pair<int,int> >::iterator IT;
			vector<float>::const_iterator fit = conAreas.begin();
			for(IT ic=connectionlist.begin();ic!=connectionlist.end();ic++)
			{
				int tonodeid=ic->first;
				int fldir=ic->second;

				Box * tonodeobj= nodegroup[tonodeid];

				int tonid = get_nodenumber(tonodeid); ///messed up here....

				//assert(tonid!=tonodeid);
				double tonodedxdydz[3] = {tonodeobj->dx, tonodeobj->dy, tonodeobj->dz};

				ja[japos] = tonid+id_offset;
				fldr[japos] = fldir;
				int drpos = abs(fldir) - 1;
				cl1[japos] = dxdydzn[drpos] / 2.0;
				cl2[japos] = tonodedxdydz[drpos] / 2.0;
				double dz = 1.0;

				if(is3d)
				{
					dz = 0.5 * (dxdydzn[2] + tonodedxdydz[2]);
					if(top!=NULL && bot!=NULL)
					{
						double dz1 = top[nid] - bot[nid];
						double dz2 = top[tonid] - bot[tonid];
						dz = 0.5 * (dz1 + dz2);
					}
				}

				if((!vertical_pass_through) || (*fit) <0 )
				{
					if(drpos == 0) //x
					{
						fahl[japos] = min(dxdydzn[1], tonodedxdydz[1]) * dz;
					}
					else if (drpos == 1) //y
					{
						fahl[japos] = min(dxdydzn[0], tonodedxdydz[0]) * dz;
					}
					else if (drpos == 2) //z
					{
						fahl[japos] = min(dxdydzn[0] * dxdydzn[1],tonodedxdydz[0] * tonodedxdydz[1]);
					}
				}
				else
				{
					assert(fldir == 3 || fldir == -3);
					fahl[japos] = *fit;
				}

				//cout<<"fahl["<<japos<<"]="<<fahl[japos]<<" dz="<<dz<<endl;

				japos ++;

				if(vertical_pass_through)
					++fit;

			}//end for ic

		}//end for nid


		ia[box_size] = japos+id_offset;

		//cout<<"- Returning the UsgCsrData object"<<endl;

		return CSRData(box_size, nja, area, ia, ja, cl1, cl2, fahl, fldr, iac);

	}

	void QuadTree3D::sortconnectionbynodenumber(list< pair<int,int> > & connectionlist)
	{    
		//get node numbers
		typedef list< pair<int,int> >::iterator IT;
		vector< pair<int, pair<int,int> >  > sorted_connection_list;
		for(IT ic=connectionlist.begin();ic!=connectionlist.end();ic++)
		{
			int tonid = get_nodenumber(ic->first);
			sorted_connection_list.push_back(make_pair(tonid,*ic));
		}
		//sort by node numbers
		sort(sorted_connection_list.begin(), sorted_connection_list.end());
		//put the sorted connection back
		list< pair<int,int> > tmp;
		int size=sorted_connection_list.size();
		for(int i=0; i<size; i++) tmp.push_back( sorted_connection_list[i].second );
		tmp.swap(connectionlist);
	}


	void QuadTree3D::sortareabyconid(const list<pair<int, int> >& connectionlist, vector<float>& area)
	{
		typedef list< pair<int,int> >::const_iterator CIT;
		vector<pair<int, float> > sorted_area_list;
		vector<float>::iterator fit = area.begin();
		for(CIT ic = connectionlist.begin(); ic != connectionlist.end(); ++ic, ++fit)
		{
			int tonid = get_nodenumber(ic->first);
			sorted_area_list.push_back(make_pair(tonid, *fit));
		}
		//sort by node numbers
		sort(sorted_area_list.begin(), sorted_area_list.end());
		//put the sorted connection back
		vector<float> tmp;
		int size = sorted_area_list.size();
		for(int i = 0; i < size; i++) tmp.push_back(sorted_area_list[i].second);
		tmp.swap(area);
	}




	//	void QuadTree3D::buildConnectionBetweenLayers()
	//	{
	//
	//	}

	void QuadTree3D::_add_root_nodes()
	{
		// Copy the nodegroup from the base modflow grid
		nodegroup=m_mfgrid->nodegroup;
		nodegroup.nlay=m_mfgrid->nlay;
		nodegroup.grid=this;

		// Copy the X and Y coordinates for the center of all cells
		// from base modflow grid.
		this->X = clone_array(m_mfgrid->X, m_mfgrid->ncol);
		this->Y = clone_array(m_mfgrid->Y, m_mfgrid->nrow);

		assert(X&&Y);
	}


	// Create the top array from the parent grid
	// Do not call this function if you just wanted to get the top array
	double * QuadTree3D::get_top()
	{
		if(this->top!=NULL)
		{
			delete [] this->top;
		}

		int nodesize = nodegroup.node_size();
		int mfgrid_csize=m_mfgrid->nrow * m_mfgrid->ncol;
		this->top = new double[nodesize];
		assert(this->top);

		for(int nodenumber=0;nodenumber<nodesize;nodenumber++)
		{
			Box * b=get_nodeobj(nodenumber);

			if(m_top_layer_z_operator[b->layer]==Z_VALUE_REPLICATE_OPEARTOR)
			{
				Index2d i=m_mfgrid->getIndex(Point2d(b->x,b->y));
				int id=m_mfgrid->get_nodeid(b->layer,i[0],i[1]);
			    this->top[nodenumber] = m_mfgrid->botm(b->layer,i[0],i[1]); 
			}
			else //if(m_top_layer_z_operator[b->layer]==Z_VALUE_LINEAR_INTERPOLATION_OPEARTOR)
			{
				if( m_top_layer_z_source.find(b->layer)==m_top_layer_z_source.end() ) //no name defined, interpolate from the based grid
					this->top[nodenumber] = bilinear_interpolation(b,&m_mfgrid->bot[mfgrid_csize*b->layer]);
				else //read from ascii file
				{
					bool area_weigthed = true;
					if(m_top_area_weighted_interpolate.find(b->layer) != m_top_area_weighted_interpolate.end())
						area_weigthed = m_top_area_weighted_interpolate[b->layer];

					this->top[nodenumber] = bilinear_interpolation(b,m_top_layer_z_source[b->layer], area_weigthed);
				}
			}
		}


		return this->top;
	}

	// create the bottom array from the parent grid
	// Do not call this function if you just wanted to get the bot array
	double * QuadTree3D::get_bot()
	{
		if(this->bot!=NULL)
		{
			delete this->bot;
		}

		int nodesize=nodes();
		int mfgrid_csize=m_mfgrid->nrow * m_mfgrid->ncol;
		this->bot = new double[nodesize];
		assert(this->bot);

		for(int nodenumber=0;nodenumber<nodesize;nodenumber++)
		{
			Box * b=get_nodeobj(nodenumber);

			if(m_bot_layer_z_operator[b->layer]==Z_VALUE_REPLICATE_OPEARTOR)
			{
				Index2d i=m_mfgrid->getIndex(Point2d(b->x,b->y));
			    this->bot[nodenumber] = m_mfgrid->botm(b->layer+1,i[0],i[1]); 
			}
			else //if(m_bot_layer_z_operator[b->layer]==Z_VALUE_LINEAR_INTERPOLATION_OPEARTOR)
			{
				if(  m_bot_layer_z_source.find(b->layer)==m_bot_layer_z_source.end() ) //no name defined, interpolate from the based grid
					this->bot[nodenumber] = bilinear_interpolation(b,&m_mfgrid->bot[mfgrid_csize*(b->layer+1)]);
				else
				{
					bool area_weigthed = true;
					if(m_bot_area_weighted_interpolate.find(b->layer) != m_bot_area_weighted_interpolate.end())
						area_weigthed = m_bot_area_weighted_interpolate[b->layer];

					//read from ascii file
					this->bot[nodenumber] = bilinear_interpolation(b,m_bot_layer_z_source[b->layer], area_weigthed);
				}
			}
		}

		return this->bot;
	}

	///recompute z value and dz value of each cell based on the top and bottom
	bool QuadTree3D::recompute_z_of_cells()
	{
		//check if bot and top are created
		if(bot==NULL || top==NULL) 
		{
			cerr<<"! Error: QuadTree3D::recompute_z_of_cells: Either top or bot are not initialized"<<endl;
			return false;
		}

		int nodesize=nodes();
		int mfgrid_csize=m_mfgrid->nrow * m_mfgrid->ncol;

		for(int nodenumber=0;nodenumber<nodesize;nodenumber++)
		{
			Box * b=get_nodeobj(nodenumber);
			b->dz = (top[nodenumber]-bot[nodenumber]);
			b->z = (top[nodenumber]+bot[nodenumber])/2;
		}
		return true;
	}

	///align the elevation (z) between the shared surface (bottom of top cell and top of bottom cells)
	bool QuadTree3D::align_elev()
	{
		//check if bot and top are created
		if(bot==NULL || top==NULL) 
		{
			cerr<<"! Error: QuadTree3D::align_elev: Either top or bot are not initialized"<<endl;
			return false;
		}

		int nodesize=nodes();

		//for each leaf node
		for(int nodenumber=0;nodenumber<nodesize;nodenumber++)
		{
			Box * box=get_nodeobj(nodenumber);
			assert(box->isLeaf); //make sure this node is a leaf node...

			if( box->active==false ) continue; //disabled...
			if( box->pDown==NULL ) continue; //no neighbors at lower level
			align_elev(box);
		}//end for nodenumber

		return true;
	}

	///align the elevation (z) between the shared surface (bottom of top cell and top of bottom cells) of the given box
	void QuadTree3D::align_elev(Box * top_box)
	{		
		//get all leaves of box->pDown				
		//TODO: we need to consider VP here...

		double bot_of_top_box=bot[top_box->number];

		list<Box*> lower_boxes;
		top_box->refine_pDown();
		top_box->pDown->getLeaves(lower_boxes);
						
		//align
		for(list<Box*>::iterator j=lower_boxes.begin();j!=lower_boxes.end();j++)
		{
			Box* bot_box=*j;
			if(bot_box->active==false) continue; //disabled...
			double top_of_bot_box=this->top[bot_box->number];
			if(bot_of_top_box<top_of_bot_box)
			{
				//cout<<"align = "<<bot_of_top_box<<" and "<<top_of_bot_box<<endl;
				this->top[bot_box->number]=bot_of_top_box;
			}
		}//end for j
	}

	///number the nodes (i.e. leaves) in the grid
	void QuadTree3D::number_nodes()
	{
		ModflowGrid * mfgrid=getModflowGrid();
		unsigned int number=0;
		for(int l=0;l<mfgrid->nlay;l++)
		{
			for(int r=0;r<mfgrid->nrow;r++)
			{
				for(int c=0;c<mfgrid->ncol;c++)
				{
					int bid=mfgrid->get_nodeid(l,r,c);
					Box * box=nodegroup[bid];
					number_nodes(box, number);
				}//end for c
			}//end for r
		}//end for l
	}

	void QuadTree3D::number_nodes(Box * box, unsigned int& number)
	{
		if(box->isLeaf)
		{
			if(box->active)
				box->number=number++;
			else
				box->number=-1;
		}
		else
		{
			//for(int i=0;i<4;i++) number_nodes(box->pChildren[i],number);
			number_nodes(box->pChildren[0],number);
			number_nodes(box->pChildren[1],number);
			number_nodes(box->pChildren[3],number);
			number_nodes(box->pChildren[2],number);
		}
	}

	//refresh the tree to get the new node id and number arrays
	void QuadTree3D::refresh()
	{
		rebuild_nodeid_array = true;
		nodeid_array=getnodeid_array();
		nodenumber_array=getnodenumber_array();
	}

	//print the information of the leaves
	inline void DFSprint
		(ostream& out, Box * box, unsigned int & id, const string& extra, string& location, bool one_based_numbering)
	{
		char tmp[8];

		if(box->isLeaf)
		{
			if(box->active){
				if(id!=box->number+((one_based_numbering)?1:0))
				{
					cerr<<"! Error: box number="<<box->number+((one_based_numbering)?1:0)<<" differs from "<<id<<endl;
				}
				assert(id==box->number+((one_based_numbering)?1:0));
				out<<id++<<", "<<extra<<" "<<location<<"\n";
			}
			else
				out<<box->number<<", "<<extra<<" "<<location<<"\n";
		}
		else
		{
			for(int x=0;x<4;x++)
			{
				int i=0;
				switch(x){
				case 0: i=0; break;
				case 1: i=1; break;
				case 2: i=3; break;
				case 3: i=2; break;
				}
				if(one_based_numbering) sprintf(tmp,"%d",x+1);
				else sprintf(tmp,"%d",x);
				string mylocation=location+tmp;
				DFSprint(out,box->pChildren[i],id,extra,mylocation,one_based_numbering);
			}//end for i
		}
	}

	ostream & operator<<(ostream& out, QuadTree3D &grid)
	{
		ModflowGrid * mfgrid=grid.getModflowGrid();
		unsigned int id_offset=(grid.one_based_node_numbering)?1:0;
		unsigned int id=id_offset; //starting from 1 if it's one based
		char extra[64];

		out<<grid.nodegroup.get_total_leaf_size()<<"\n";

		/* 
		   //
		   // the foloowing code is removed on 3/25/2014. 
		   // this is now output to .nod file
		   // in bool DefinitionExporter::save_node_coordinate(const string& prefix, Grid * grid)
		   // in def_exporter.cpp
		   //

		//output number of nodes and number of connections
		//this is the number of active leaves of the quadtree
		int box_size=grid.nodes();

		//go through each cell and count diagonal position and connections
		int nja = 0;
		for( int nid=0;nid<box_size;nid++)  //nid in range(self.nodes) )
		{
			nja++;
			int csize=grid.get_node_connections(nid).size();
			nja += csize;
		}

		out<<box_size<<"\n";
		out<<nja<<"\n";
		//
		*/

		for(int l=0;l<mfgrid->nlay;l++)
		{
			for(int r=0;r<mfgrid->nrow;r++)
			{
				for(int c=0;c<mfgrid->ncol;c++)
				{
					int bid=mfgrid->get_nodeid(l,r,c);
					Box * box=grid.nodegroup[bid];
					sprintf(extra,"(%d,%d,%d)",l+id_offset,r+id_offset,c+id_offset);
					string extra_str=extra;
					string location;
					DFSprint(out,box,id,extra_str,location,grid.one_based_node_numbering);
				}//end for c
			}//end for r
		}//end for l

		return out;
	}

	void QuadTree3D::saveGhostNode(string filePath)
	{
		if(m_gns.size() == 0)
			computeGNC(this, m_gns);

		writeGhostNode(filePath, m_gns, ((this->one_based_node_numbering)?1:0) );
		cout<<"- Save data to file: "<<filePath<<endl;
	}

	//linear interploation the value of the center of the box
	double QuadTree3D::bilinear_interpolation(Box * box, double * zgrid)
	{
		//create a ascii_grid from grid 
		mf_ascii_grid ag(this->getModflowGrid(),zgrid,m_mfgrid->ncol,m_mfgrid->nrow);
		Point2d pt(box->x,box->y);
		return bilinear_interpolation(pt, &ag);
	}

	//linear interploation the value of the center of the box
	double QuadTree3D::bilinear_interpolation(Box * box, const string & arcinfo_filename, bool useAreaWeighted)
	{
		ascii_grid * ag=m_str2ag[arcinfo_filename];

		if(ag==NULL) //not found, build from filename
		{
			arcinfo_ascii_grid * aag=new arcinfo_ascii_grid();
			assert(aag);
			bool br=aag->build(arcinfo_filename);
			if(br==false)
			{
				cerr<<"! Error: failed to build arcinfo_ascii_grid from file: ("<<arcinfo_filename<<")"<<endl;
				assert(false);
				exit(1);
			}

			m_str2ag[arcinfo_filename]=ag=aag;
		}

		//#if 0
		if(!useAreaWeighted)
		{
			//rotate box center if the grid is rotated...
			Point2d pt;
			{
				double cx, cy, rotation;
				getModflowGrid()->getRotatePara(cx, cy, rotation);
				double rx,ry;
				box->rotateCtr(cx,cy,rotation,rx,ry);
				pt.set(rx,ry);
			}

			//interpolate here
			return bilinear_interpolation(pt, ag);
		}
		//#endif 
		//
		//#if 1
		else
		{
			//change to use the area-weighted interpolation
			Point2d corners[4];
			corners[0].set(min(box->x - box->dx/2, box->x + box->dx/2), min(box->y - box->dy/2, box->y + box->dy/2));
			corners[1].set(max(box->x - box->dx/2, box->x + box->dx/2), min(box->y - box->dy/2, box->y + box->dy/2));
			corners[2].set(max(box->x - box->dx/2, box->x + box->dy/2), max(box->y - box->dy/2, box->y + box->dy/2));
			corners[3].set(min(box->x - box->dx/2, box->x + box->dx/2), max(box->y - box->dy/2, box->y + box->dy/2));

			int dminx = ag->getncols() + 1, dmaxx = -1, dminy = ag->getnrows() + 1, dmaxy = -1;
			double cx, cy, rotation;
			getModflowGrid()->getRotatePara(cx, cy, rotation);
			for(int i = 0; i < 4; i++)
			{
				double tx, ty;
				tx = cx + (corners[i][0] - cx) * cos(rotation) - (corners[i][1] - cy)*sin(rotation);
				ty = cy + (corners[i][1] - cy) * cos(rotation) + (corners[i][0] - cx)*sin(rotation);
				corners[i].set(tx, ty);

				Index2d idx=ag->getIndex(corners[i][0], corners[i][1]);
				if(idx[1] < dminx)
					dminx = idx[1];
				if(idx[1] > dmaxx)
					dmaxx = idx[1];
				if(idx[0] < dminy)
					dminy = idx[0];
				if(idx[0] > dmaxy)
					dmaxy = idx[0];
			}

			dminx = max(0, dminx); dmaxx = min((int)(ag->getncols() - 1), dmaxx);
			dminy = max(0, dminy); dmaxy = min((int)(ag->getnrows() - 1), dmaxy);

			double val = 0.0, totalArea = 0.0;
			for(int ix = dminx; ix <= dmaxx; ix++)
			{
				for(int iy = dminy; iy <= dmaxy; iy++)
				{
					ascii_grid_cell cell=ag->getCell(iy,ix);
					double area = intersectTwoBox(corners[0], corners[1], corners[2], corners[3], 
						min(cell.x - cell.dx/2, cell.x + cell.dx/2), max(cell.x - cell.dx/2, cell.x + cell.dx/2), 
						min(cell.y - cell.dy/2, cell.y + cell.dy/2), max(cell.y - cell.dy/2, cell.y + cell.dy/2));
					val += area * ag->getvalue(iy, ix);

					totalArea += area;
				}
			}

			//if(totalArea != fabs(box->dx * box->dy))
			//{
			//	cout<<totalArea<<" - "<<(box->dx * box->dy)<<"\n";
			//	cout<<"interpolation area error!\n";
			//}

			//if total intersection area is 0, use the bilinear
			//rotate box center if the grid is rotated...

			if(totalArea < 1e-6)
			{
				Point2d pt;
				double cx, cy, rotation;
				getModflowGrid()->getRotatePara(cx, cy, rotation);
				double rx,ry;
				box->rotateCtr(cx,cy,rotation,rx,ry);
				pt.set(rx,ry);
				//interpolate here
				return bilinear_interpolation(pt, ag);
			}
			else
				return val = val / totalArea;
		}
		//#endif

	}

	//linear interploation the value of a given point in a given layer (either top or bottom)
	double QuadTree3D::bilinear_interpolation(const Point2d& pt, int layer, bool top)
	{
		int mfgrid_csize=m_mfgrid->nrow * m_mfgrid->ncol;

		if(top)
		{
			if( m_top_layer_z_source.find(layer)==m_top_layer_z_source.end() ) //no name defined, interpolate from the based grid
				return bilinear_interpolation(pt,&m_mfgrid->bot[mfgrid_csize*layer]);
			else //read from ascii file
				return bilinear_interpolation(pt,m_top_layer_z_source[layer]);
		}
		else
		{
			if(  m_bot_layer_z_source.find(layer)==m_bot_layer_z_source.end() ) //no name defined, interpolate from the based grid
				return bilinear_interpolation(pt,&m_mfgrid->bot[mfgrid_csize*(layer+1)]);
			else //from arc ascii info file
				return bilinear_interpolation(pt,m_bot_layer_z_source[layer]);
		}
	}

	//linear interploation the value of a given point
	double QuadTree3D::bilinear_interpolation(const Point2d& pt, double * zgrid)
	{
		//create a ascii_grid from grid 
		mf_ascii_grid ag(this->getModflowGrid(),zgrid,m_mfgrid->ncol,m_mfgrid->nrow);
		return bilinear_interpolation(pt, &ag);
	}

	//linear interploation the value of a given point 
	double QuadTree3D::bilinear_interpolation(const Point2d& pt, const string & arcinfo_filename)
	{
		ascii_grid * ag=m_str2ag[arcinfo_filename];

		if(ag==NULL) //not found, build from filename
		{
			arcinfo_ascii_grid * aag=new arcinfo_ascii_grid();
			assert(aag);
			bool br=aag->build(arcinfo_filename);
			if(br==false)
			{
				cerr<<"! Error: failed to build arcinfo_ascii_grid from file: ("<<arcinfo_filename<<")"<<endl;
				assert(false);
				exit(1);
			}

			m_str2ag[arcinfo_filename]=ag=aag;
		}


		//rotate pt if the grid is rotated...
		Point2d pt2;
		{
			double cx, cy, rotation;
			getModflowGrid()->getRotatePara(cx, cy, rotation);
			double rx,ry;
			rx= cx +  (pt[0] - cx) * cos(rotation) - (pt[1] - cy)*sin(rotation);
			ry = cy + (pt[1] - cy) * cos(rotation) + (pt[0] - cx)*sin(rotation);
			pt2.set(rx,ry);
		}

		//interpolate here
		return bilinear_interpolation(pt2, ag);
	}

	//// GL: bilinear interploation 
	//double QuadTree3D::bilinear_interpolation(const Point2d& pt, ascii_grid * ag)
	//{
	//	Index2d id = ag->getIndex(pt[0], pt[1]);
	//	Index2d id1 = ag->getIndex(pt[0], pt[1]);

	//	int i=nrows-1-floor((y-yllcorner)/cellsize);
	//	int j=floor((x-xllcorner)/cellsize);

	//
	
	//	assert(false);
	//	return 0;
	//}


	
	//linear interploation the value of the center of the box
	double QuadTree3D::bilinear_interpolation(const Point2d& pt, ascii_grid * ag)
	{
		//determine the i and j to interploate from...

		//get the ag cell that contains pt
		Index2d id=ag->getIndex(pt[0],pt[1]);

#if 1
		// return nodata when the point is out of bound...

		ascii_grid_cell cell=ag->getCell(id[0],id[1]);
			

		//
		if(pt[0]<cell.x && id[1]>0 ) id[1]--;
		//if(pt[1]<=cell.y && id[0]<(int)ag->getnrows()-1 ) id[0]++;
		//if(pt[1]>=cell.y && id[0]<(int)ag->getnrows()-1 ) id[0]++;
		if(pt[1]>cell.y && id[0]>0 ) id[0]--;

		//if(pt[1]<=cell.y && id[0]>0 ) id[0]--;
		
		cell=ag->getCell(id[0],id[1]);

		//bilinear_interpolation
		double x= (pt[0]-cell.x)/cell.dx;
		double y= (cell.y-pt[1])/cell.dy;
		//double y= (pt[1]-cell.y)/cell.dy;


		if(x<0 || x>1 || y<0 || y>1)
		{
			Index2d id=ag->getIndex(pt[0],pt[1]);
			{
				if(x<0) x=0;
				else if(x>1) x=1;
				if(y<0) y=0;
				else if(y>1) y=1;
			}
			//cout<<"degenerate case: "<<x<<", "<<y<<endl;
		}


		//if(x<0 || x>1 || y<0 || y>1) 
		//{
		//	Index2d id=ag->getIndex(pt[0],pt[1]);

		//	//if(id[1]>0 && id[0]<ag->getnrows()-1 )
		//	//{
		//		//cout<<"Error:QuadTree3D::bilinear_interpolation: Index out of bound: x="<<x<<" y="<<y<<endl;
		//		//Index2d id=ag->getIndex(pt[0],pt[1]);
		//	//}
		//	//else //
		//	{
		//		if(x<0) x=0;
		//		else if(x>1) x=1;
		//		if(y<0) y=0;
		//		else if(y>1) y=1;
		//	}

		//}

		return cusg::bilinear_interpolation(ag,id[0],id[1],x,y);
#else
		// return nodata when the point is out of bound...

#if 0   // turn off this for now... 
		// (to have the following code working, getIndex(double x, double y) needs to be modified as well)
		if( id[0]>=(int)ag->getnrows() || id[1]>=(int)ag->getncols() || id[0]<0 || id[1]<0)
		{
			return ag->getNoDATA_value();
		}
#endif
		//

		ascii_grid_cell cell=ag->getCell(id[0],id[1]);

		//
		if(pt[0]<=cell.x && id[1]>0 ) id[1]--;
		if(pt[1]<=cell.y && id[0]<(int)ag->getnrows()-1 ) id[0]++;
		//if(pt[1]>=cell.y && id[0]<(int)ag->getnrows()-1 ) id[0]++;

		//if(pt[1]<=cell.y && id[0]>0 ) id[0]--;
		
		cell=ag->getCell(id[0],id[1]);

		//bilinear_interpolation
		double x= (pt[0]-cell.x)/cell.dx;
		double y= (pt[1]-cell.y)/cell.dy;

		if(x<0 || x>1 || y<0 || y>1) 
		{
			Index2d id=ag->getIndex(pt[0],pt[1]);

			//if(id[1]>0 && id[0]<ag->getnrows()-1 )
			//{
				//cout<<"Error:QuadTree3D::bilinear_interpolation: Index out of bound: x="<<x<<" y="<<y<<endl;
				//Index2d id=ag->getIndex(pt[0],pt[1]);
			//}
			//else //
			{
				if(x<0) x=0;
				else if(x>1) x=1;
				if(y<0) y=0;
				else if(y>1) y=1;
			}

		}

		return cusg::bilinear_interpolation(ag,id[0],id[1],x,y);
#endif

	}

}//end of namespace
