/*! \file GridIntersection.h

\brief a class for computing the intersection between the grid and a geometric object

\author <a href="http://masc.cs.gmu.edu">MASC group</a>, George Mason University
\author <a href="http://profile.usgs.gov/langevin/">Christian Langevin</a>, USGS
\bug    No known bugs.
*/

#pragma once

#include "ModflowGrid.h"
#include "QuadTree3D.h"
#include "intersection.h"
#include "lineseg.h"
#include "Point.h"
#include "polygon.h"
#include "polyline.h"
#include "def_struct.h"
#include "shpReader.h"
#include "BoxPlyClip.h"
#include "PolygonChecker.h"
#include <map>
using namespace std;


#define DEBUGTEST 0

namespace cusg
{

	inline double getSqLength(const LineSeg2d& line)
	{
		double d1 = (line.pt[0][0] - line.pt[1][0]) *  (line.pt[0][0] - line.pt[1][0]);
		double d2 =  (line.pt[0][1] - line.pt[1][1]) *  (line.pt[0][1] - line.pt[1][1]);
		return d1 + d2;
	}


	struct PointCellIntersectionInfo
	{
		int m_ptID;
		Box * box;
		GIS_Point2d pnt;
	};

	struct PolylineCellIntersectionInfo
	{
		int m_curveID;
		double m_len;
		double m_startGeoDist;
		double m_endGeoDist;
		Box* box;
		c_plyline* m_plyline; //the cutted polyline
	};

	struct SegmentCellIntersectionInfo
	{
		SegmentCellIntersectionInfo()
		{
			m_curveID = m_segID = -1;
			m_len = m_startGeoDist = m_endGeoDist = 0;
			box = NULL;
		}

		int m_curveID;         //id of the curve that contains this segment
		int m_segID;           //id of the segment
		double m_len;          //length of the segment inside the box
		double m_startGeoDist;
		double m_endGeoDist;
		Box * box;
		Point2d m_startPnt; // the starting point
		Point2d m_endPnt; // the ending point
	};

	struct PolygonCellIntersectionInfo
	{
		PolygonCellIntersectionInfo()
		{
			m_polyID=-1;
			m_area=0;
			box=NULL;
			m_interPolygon=NULL;
		}

		int m_polyID;
		double m_area; //intersection area
		Box * box;
		GIS_polygon*  m_interPolygon;
	};
	void destroyPlylineCell(PolylineCellIntersectionInfo& info);

	void destroyPlygonCell(PolygonCellIntersectionInfo& info);


	class GridIntersection
	{
	public:

		GridIntersection(int max_level, int at_layer)
		{
			m_max_level=max_level;
			m_layer=at_layer;
			m_b_include_poly_boundary=true;
			m_b_save_intersection_info=false;
			//m_cur_id = -1;
		}
		~GridIntersection()
		{
			for(vector<PolygonCellIntersectionInfo>::iterator pit = m_ply_interaction_info.begin();
				pit != m_ply_interaction_info.end(); ++pit)
			{
				pit->m_interPolygon->destroy();

				delete pit->m_interPolygon;
				pit->m_interPolygon = NULL;
			}
			for(vector<PolylineCellIntersectionInfo>::iterator lit = m_plyline_intersection_info.begin(); lit != m_plyline_intersection_info.end(); ++lit)
			{
				lit->m_plyline->destroy();

				delete lit->m_plyline;
				lit->m_plyline = NULL;
			}
		}

		///turn on/off m_b_include_poly_boundary to enable/disable storing
		///the intersected boxes between the grid and the given polygon
		void set_include_poly_boundary_flag(bool flag){ m_b_include_poly_boundary=flag; }

		///turn on/off m_b_save_intersection_info to enable/disable storing
		///the detailed intersection information
		void set_save_intersection_info(bool flag) { m_b_save_intersection_info=flag; }

		//build intersection from grid_intersection raw data
		void build(grid_intersection_raw_data& data);

		//save intersection information to file
		bool save(string& filename, vector<string>& attributes);
		void save(ostream& out, vector<string>& attributes);


		///get intersected boxes
		vector<Box*>& getIntersections() { return m_intersection; }

		//save the shape file
		bool saveShp(string& shpName, vector<string>& attributes);

		//save the vtk file
		bool saveVTK(string& vtkName, vector<string>& attributes);

		//
		bool intersect(QuadTree3D * grid, double x, double y)
		{
			return intersect(grid,Point2d(x,y));
		}

		bool intersect(QuadTree3D * grid, double x1, double y1, double x2, double y2)
		{
			Point2d p1(x1,y1);
			Point2d p2(x2,y2);
			LineSeg2d s(p1,p2);

			return intersect(grid,s);
		}

		bool intersect(QuadTree3D * grid, vector<double>& ply_pts)
		{
			//create ply
			c_ply ply(c_ply::POUT);
			int size=ply_pts.size();
			ply.beginPoly();
			for(int i=0;i<size;i++)
			{
				double x=ply_pts[i];
				double y=ply_pts[++i];
				ply.addVertex(x,y);
			}
			ply.endPoly();

			return intersect(grid,ply);
		}

		Box * getIntersections(int i) { return m_intersection[i]; }

		// in all intersection functions, true is returned if
		// there are intersections; otherwise false is returned

		//
		// ModflowGrid2D
		//
		bool intersect(ModflowGrid2D * grid, const Point2d& pt);


		//
		// ModflowGrid, onBound=true means ignore the boxes only intersecting at box boundary
		//
		bool intersect(ModflowGrid * grid, const Point2d& pt);//, bool noBound=false);

		//
		//QuadTree3d, onBound=true means ignore the boxes only intersecting at box boundary
		//
		bool intersect(QuadTree3D * grid, const Point2d& pt);//, bool noBound=false);


		//intersect a generic Grid with a line segment
		template<typename T> bool intersect(T * grid, const LineSeg2d& line);//, bool noBound=false);

		//intersect a generic Grid with a polygon..
		template<typename T> bool intersect(T * grid, const c_ply& ply);//, bool noBound=false);
		template<typename T> bool intersect(T * grid, const c_polygon& polygon);//, bool noBound=false);

		//intersect a generic Grid with a ply..
		template<typename T> bool intersect
			(T * grid, const c_ply& ply, list<Box*>& bd_intersection, list<Box*>& enclosed_intersection);

		//intersect a generic Grid with a polygon..
		template<typename T> bool intersect
			(T * grid, const c_polygon& polygon, list<Box*>& bd_intersection, list<Box*>& enclosed_intersection);

		///Box
		bool intersect(Box * b, const Point2d& pt);
		bool intersect(Box * b, const LineSeg2d& line);
		double intersect_len(Box * b, const LineSeg2d& line);
		//double intersect_polygon(Box* b, const GIS_polygon& ply, GIS_polygon& resultPly);
		//void intersect_polygon(Box* b, const GIS_polygon& plygon, vector<GIS_polygon*>& resultPlygons);
		//void intersect_polygon(Box* b, BoxPlyClip* bpc, vector<GIS_polygon*>& resultPlygons);
		void intersect_polygon(Box* b, const GIS_polygon& plygon,  vector<GIS_polygon*>& resultPlygons,BoxPlyClip* bpc=NULL);

		//retrive detailed information
		bool getIntersection(Box * b, const LineSeg2d& line, LineSeg2d& intseg);
		//c_ply getIntersection(Box * b, const c_ply& ply);

		//get the number of intersections...
		int getIntersectionSize() { return m_intersection.size(); }

		//intersection for GIS based primitives
		template<typename T> bool intersect(T * grid, /*const*/ GIS_Point2d& pt);
		template<typename T> bool intersect(T * grid, /*const*/ GIS_LineSeg2d& line);
		template<typename T> bool intersect(T * grid, /*const*/ GIS_polygon& ply, BoxPlyClip * bpc = NULL);
		//template<typename T> bool intersect(T* grid, BoxPlyClip* bpc);

		/// clip off a line segment if its end points are out side of the grid
		/// return false when the entire segment is outside of the grid
		bool clip_line(ModflowGrid2D * grid, const LineSeg2d& line, LineSeg2d& clipped_line);
		bool clip_line(QuadTree3D * grid, const LineSeg2d& line, LineSeg2d& clipped_line);

		//utility functions
		int checkNeighbor(Box* b1, Box* b2);
		//check box edge collinearly intersect with line
		bool checkEdgeCollinear(Box* box, int dir, const LineSeg2d& line);

	protected:

		//check if the given pt is inside the ply
		bool intersect(const Point2d& pt, const c_ply& ply);

		//check if the given pt is inside the polygon
		bool intersect(const Point2d& pt, const c_polygon& polygon);

		//create information for intersection
		//the given box should be intersecting with the given geometry already
		void buildIntersectionInfo(Box * box, const GIS_polygon& ply, vector<PolygonCellIntersectionInfo>& infoResult, BoxPlyClip* bpc = NULL);

		void buildIntersectionInfo(Box * box, const GIS_polygon& ply, PolygonCellIntersectionInfo& info);

		//void buildIntersectionInfo(Box * box, BoxPlyClip* bpc, vector<PolygonCellIntersectionInfo>& infoResult);
		SegmentCellIntersectionInfo buildIntersectionInfo(Box * box, const GIS_LineSeg2d& line);
		PointCellIntersectionInfo buildIntersectionInfo(Box * box, const GIS_Point2d& pt);

		void floodGrid(Box * sb, set<Box*>& visited, list<Box*>& flooded);

		//given a box and two points, one inside and one outside. 
		//close the gap between in and out so that their distance is smaller than
		//a given threshold
		void bracket_to_box_boundary(Box & box, Point2d& out, Point2d& in) const;

		int m_max_level;
		int m_layer;
		vector<Box*> m_intersection; /// a list of intersected boxes

		bool m_b_include_poly_boundary;

		bool m_b_save_intersection_info;

		vector<PolygonCellIntersectionInfo> m_ply_interaction_info;
		vector<SegmentCellIntersectionInfo> m_seg_interaction_info;
		vector<PointCellIntersectionInfo>   m_pt_interaction_info;
		vector<PolylineCellIntersectionInfo> m_plyline_intersection_info;

		map<Box*, list<ply_vertex *> > boxvertexMap; //mapping a box to a list of polygon vertices/edges inside the box

		//
		vector<Box*> m_onBound;


	public:
		ShapeReader m_shp_reader;

		vector<Box*> testBoxes;

		vector<PolygonCellIntersectionInfo>& getPolygonIntersectionInfo()
		{
			return m_ply_interaction_info;
		}

		vector<SegmentCellIntersectionInfo>& getSegmentIntersectionInfo()
		{
			return m_seg_interaction_info;
		}

		vector<PointCellIntersectionInfo>& getPointIntersectionInfo()
		{
			return m_pt_interaction_info;
		}

		vector<PolylineCellIntersectionInfo>& getPolylineIntersectionInfo()
		{
			return m_plyline_intersection_info;
		}

	};

	//
	//
	// implementation of templated functions
	//
	//

	//intersection a generic Grid with a polygon..
	template<typename T> bool GridIntersection::intersect(T * grid, const c_ply& ply)//, bool noBound)
	{
		list<Box*> bd_intersection;
		list<Box*> enclosed_intersection;
		bool r=intersect(grid,ply,bd_intersection,enclosed_intersection);

		if(r)
		{
			//add the newly discovered intersections
			if(m_b_include_poly_boundary)
				m_intersection.insert(m_intersection.end(),bd_intersection.begin(),bd_intersection.end());

			m_intersection.insert(m_intersection.end(),enclosed_intersection.begin(),enclosed_intersection.end());
		}

		return r;
	}


	//intersection a generic Grid with a polygon..
	template<typename T> bool GridIntersection::intersect(T * grid, const c_polygon& polygon) //, bool noBound)
	{
		list<Box*> bd_intersection;
		list<Box*> enclosed_intersection;
		bool r=intersect(grid,polygon,bd_intersection,enclosed_intersection);
		if(r)
		{
			//add the newly discovered intersections
			if(m_b_include_poly_boundary)
				m_intersection.insert(m_intersection.end(),bd_intersection.begin(),bd_intersection.end());
			m_intersection.insert(m_intersection.end(),enclosed_intersection.begin(),enclosed_intersection.end());
		}
		return r;
	}

	//intersection a Grid with a ply..
	template<typename T> bool GridIntersection::intersect
		(T * grid, const c_ply& ply, list<Box*>& bd_intersection, list<Box*>& enclosed_intersection)
	{
		vector<Box*> backup;

		backup.swap(m_intersection);
		set<Box*> visited;

		//get all boundary cells
		ply_vertex * ptr=ply.getHead();
		do{
			ply_vertex * next=ptr->getNext();

			LineSeg2d line(ptr->getPos(),next->getPos());
			bool r=intersect(grid, line);
			assert(r);

			//record the relationship between the boxes and this segment
			if(m_b_save_intersection_info)
			{
				for(vector<Box*>::iterator i=m_intersection.begin();i!=m_intersection.end();i++)
					boxvertexMap[*i].push_back(ptr);
			}

			visited.insert(m_intersection.begin(),m_intersection.end());
			m_intersection.clear();
			ptr=next;
		}
		while(ptr!=ply.getHead());

		if(visited.empty()) return false;

		//border intersections
		bd_intersection.insert(bd_intersection.end(),visited.begin(),visited.end());

		//get all enclosed cells
		list<Box*> open=bd_intersection;
		while(open.empty()==false) //not empty
		{
			Box * b=open.front();
			open.pop_front();

			//collect neighboring cells
			list<Box*> neighbors[4];

			for (int i = 0; i < 4; i++)
			{
				b->getNeighbors(neighbors[i], i);
			}

			//check if cells are inside the polygon
			for(int i=0;i<4;i++)
			{
				for(list<Box*>::iterator n=neighbors[i].begin(); n!=neighbors[i].end(); n++)
				{
					Box * nb=*n;

					if(visited.find(nb)!=visited.end()) continue; //visited
					Point2d b_o(nb->x,nb->y); //center of the box
					if( intersect(b_o,ply) ) //if the center of the box is inside...
					{
						open.push_back(nb);
						visited.insert(nb);
						enclosed_intersection.push_back(nb);

						/*
						//now we should just flood the entire grid with this interior cell
						list<Box*> flooded;
						floodGrid(nb,visited,flooded);

						open.insert(open.end(),flooded.begin(),flooded.end());
						visited.insert(flooded.begin(),flooded.end());
						enclosed_intersection.insert(enclosed_intersection.end(),flooded.begin(),flooded.end());
						*/
					}
				}//end for n
			}//end for i

		}//end while


		//restore
		backup.swap(m_intersection);

		return true;
	}


	//intersection a Grid with a polygon..
	template<typename T> bool GridIntersection::intersect
		(T * grid, const c_polygon& polygon, list<Box*>& bd_intersection, list<Box*>& enclosed_intersection)
	{
		vector<Box*> backup;

		backup.swap(m_intersection);
		set<Box*> visited;

		//get all boundary cells
		for(c_polygon::const_iterator i=polygon.begin(); i!=polygon.end();i++)
		{
			const c_ply& ply=*i;

			ply_vertex * ptr=ply.getHead();

			do{
				ply_vertex * next=ptr->getNext(); //ptr is the current vertex and next is the next vertex in this ply

				LineSeg2d line(ptr->getPos(),next->getPos());
				bool r=intersect(grid, line);

#if DEBUGTEST
				if (r == false)
				{
					intersect(grid, line);
				}
#endif

				assert(r);

				//record the relationship between the boxes and this segment
				if(m_b_save_intersection_info)
				{
					for(vector<Box*>::iterator ib=m_intersection.begin();ib!=m_intersection.end();ib++)
						boxvertexMap[*ib].push_back(ptr);
				}

				visited.insert(m_intersection.begin(),m_intersection.end());
				m_intersection.clear();
				ptr=next;
			}
			while(ptr!=ply.getHead());
		}

		if(visited.empty()) return false;

		//border intersections
		bd_intersection.insert(bd_intersection.end(),visited.begin(),visited.end());

		//get all enclosed cells
		list<Box*> open=bd_intersection;
		while(open.empty()==false) //not empty
		{
			Box * b=open.front();
			open.pop_front();

			//collect neighboring cells
			list<Box*> neighbors[4];

			for (int i = 0; i < 4; i++)
			{
				b->getNeighbors(neighbors[i], i);
			}

			//check if cells are inside the polygon
			for(int i=0;i<4;i++)
			{
				for(list<Box*>::iterator n=neighbors[i].begin(); n!=neighbors[i].end(); n++)
				{
					Box * nb=*n;

					if(visited.find(nb)!=visited.end()) continue; //visited
					Point2d b_o(nb->x,nb->y); //center of the box
					if( intersect(b_o,polygon) ) //if the center of the box is inside...
					{ 
						open.push_back(nb);
						visited.insert(nb);
						enclosed_intersection.push_back(nb);

						/*
						//now we should just flood the entire grid with this interior cell
						list<Box*> flooded;
						floodGrid(nb,visited,flooded);

						open.insert(open.end(),flooded.begin(),flooded.end());
						visited.insert(flooded.begin(),flooded.end());
						enclosed_intersection.insert(enclosed_intersection.end(),flooded.begin(),flooded.end());
						*/
					}
				}//end for n
			}//end for i

		}//end while


		//restore
		backup.swap(m_intersection);

		return true;

	}


	//modify here to solve the problem of non-symmetrical refinement

	template<typename T>
	bool GridIntersection::intersect(T * grid, const LineSeg2d& line)//, bool noBoud)
	{
		LineSeg2d clipped_line;
		bool rclip=clip_line(grid,line,clipped_line);
		if(rclip==false) return true; //line segment is entirely outside, so there is no intersection
		if(getSqLength(clipped_line) < 1e-10)
			return true;

		//get the first box that encloses one of the end pt
		bool r=intersect(grid,clipped_line.pt[0]);
		if(r==false){
			cerr<<"! Error: GridIntersection::intersect(ModflowGrid2D * grid, const LineSeg2d& line) #1: Pt="<<line.pt[0]<<endl;
			return false;
		}

		//first box
		Box * b1=m_intersection.back();

		//shpArcId.insert(make_pair(b1, m_cur_id));

		r=intersect(grid,clipped_line.pt[1]);
		if(r==false){
			cerr<<"! Error: GridIntersection::intersect(ModflowGrid2D * grid, const LineSeg2d& line) #2: Pt="<<line.pt[1]<<endl;
			return false;
		}

		//last box
		Box * b2=m_intersection.back();

		//the line segment is short and contained in the same box
		if(b1==b2){
			m_intersection.pop_back();
			return true; //done
		}

		//check if b1 and b2 share faces and the line is just colinear with that
		int neiID = checkNeighbor(b1, b2);
		if(neiID != -1)
		{
			//check if the box edge collinearly intersect with the line
			if(checkEdgeCollinear(b1, neiID, line))
			{
				m_intersection.pop_back();
				return true; //done
			}
		}

		//now, visit b2 neighbors recursively until we reaches b1
		list<Box*> open;
		open.push_back(b1);
		open.push_back(b2);

		set<Box*> visited;
		visited.insert(b1);
		visited.insert(b2);

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

					if( intersect(nb,line) )
					{
						//////////////////////////////////////////////////////////////////////////
						//add by Guilin
						//if collinear intersection
						if(checkEdgeCollinear(b, i, line))
						{
							cout<<"warning: edge collinear with the line\n";
							continue;
						}

						double epsilon=min(nb->dx,nb->dy)/1000;
						if(intersect_len(nb,line)>epsilon) //long enough...
						{

							m_intersection.push_back(nb);
						}

						open.push_back(nb);
						visited.insert(nb);
					}
				}
			}

		}//end while

		return true;
	}





	/**************************************************************************************************/
	//intersection for GIS based primitives
	template<typename T> bool GridIntersection::intersect(T * grid, /*const*/ GIS_Point2d& pt)
	{
		////rotate
		//ModflowGrid* mfgrid = dynamic_cast<ModflowGrid*>(grid);
		//if(mfgrid==NULL)
		//{
		//	QuadTree3D* quadgrid = dynamic_cast<QuadTree3D*>(grid);
		//	if(quadgrid!=NULL)
		//		mfgrid = quadgrid->getModflowGrid();
		//}
		//if(mfgrid == NULL)
		//{
		//	LGR_Grid* lgrgrid = dynamic_cast<LGR_Grid*>(grid);
		//	if(lgrgrid != NULL)
		//		mfgrid = lgrgrid->getModflowGrid();
		//}
		//if(mfgrid == NULL)
		//{
		//	cerr<<"error: abnormal grid !"<<endl; 
		//}
		//double cx, cy, angle;
		//mfgrid->getRotatePara(cx, cy, angle);
		////rotate
		//pt.rotate(cx, cy, -angle);


		vector<Box*> tmp;
		m_intersection.swap(tmp);
		bool r=intersect(grid,(const Point2d&)pt);

		//process information
		if(this->m_b_save_intersection_info && r)
		{
			//store intersection info here...
			for(vector<Box*>::iterator i=m_intersection.begin();i!=m_intersection.end();i++)
			{
				m_pt_interaction_info.push_back(buildIntersectionInfo(*i,pt));
			}
		}

		tmp.insert(tmp.end(),m_intersection.begin(),m_intersection.end());
		m_intersection.swap(tmp);
		return r;
	}

	template<typename T> bool GridIntersection::intersect(T * grid, /*const*/ GIS_LineSeg2d& line)
	{
		////rotate
		//ModflowGrid* mfgrid = dynamic_cast<ModflowGrid*>(grid);
		//if(mfgrid==NULL)
		//{
		//	QuadTree3D* quadgrid = dynamic_cast<QuadTree3D*>(grid);
		//	if(quadgrid!=NULL)
		//		mfgrid = quadgrid->getModflowGrid();
		//}
		//if(mfgrid == NULL)
		//{
		//	LGR_Grid* lgrgrid = dynamic_cast<LGR_Grid*>(grid);
		//	if(lgrgrid != NULL)
		//		mfgrid = lgrgrid->getModflowGrid();
		//}
		//if(mfgrid == NULL)
		//{
		//	cerr<<"error: abnormal grid !"<<endl; 
		//}
		//double cx, cy, angle;
		//mfgrid->getRotatePara(cx, cy, angle);
		////rotate
		//line.rotate(cx, cy, -angle);

		vector<Box*> tmp;
		m_intersection.swap(tmp);
		bool r=intersect(grid,(const LineSeg2d&)line);

		//process information
		if(this->m_b_save_intersection_info && r)
		{
			//store intersection info here...
			for(vector<Box*>::iterator i=m_intersection.begin();i!=m_intersection.end();i++)
			{
				m_seg_interaction_info.push_back( buildIntersectionInfo(*i,line) );
			}//end for i

		}//end if

		tmp.insert(tmp.end(),m_intersection.begin(),m_intersection.end());
		m_intersection.swap(tmp);
		return r;
	}

	template<typename T> bool GridIntersection::intersect(T * grid, /*const*/ GIS_polygon& ply, BoxPlyClip* bpc)
	{

		list<Box*> bd_intersection;
		list<Box*> enclosed_intersection;
		bool r=intersect(grid,(const c_ply&)ply.front(),bd_intersection,enclosed_intersection);

		if(r)
		{
			//handle the enclosed boxes
			m_intersection.insert(m_intersection.end(),enclosed_intersection.begin(),enclosed_intersection.end());

			//process intersection information
			if(this->m_b_save_intersection_info)
			{
				//store intersection info here...
				for(list<Box*>::iterator i=enclosed_intersection.begin();i!=enclosed_intersection.end();i++)
				{ 
					Box * box=*i;

#if DEBUGTEST
					if(box->number != 8408)
						continue;
#endif

					//store intersection info here...
#if 0
					PolygonCellIntersectionInfo info;
					buildIntersectionInfo(box, ply, info);
					m_ply_interaction_info.push_back(info);
#else
					vector<PolygonCellIntersectionInfo> infos;
					buildIntersectionInfo(box, ply, infos, bpc);
					m_ply_interaction_info.insert(m_ply_interaction_info.end(), infos.begin(), infos.end());
#endif
				}
			}

			//handle the boundary boxes if needed
			if(m_b_include_poly_boundary)
			{
				m_intersection.insert(m_intersection.end(),bd_intersection.begin(),bd_intersection.end());
				//process information
				if(this->m_b_save_intersection_info)
				{
					//store intersection info here...
					for(list<Box*>::iterator i=bd_intersection.begin();i!=bd_intersection.end();i++)
					{
						Box * box=*i;

#if DEBUGTEST
						if(box->number != 8408)
							continue;
#endif

						//store intersection info here...
#if 0
						PolygonCellIntersectionInfo info;
						buildIntersectionInfo(box, ply, info);
						m_ply_interaction_info.push_back(info);
#else
						vector<PolygonCellIntersectionInfo> infos;
						buildIntersectionInfo(*i, ply, infos, bpc);
						m_ply_interaction_info.insert(m_ply_interaction_info.end(), infos.begin(), infos.end());
#endif

					}
				}
			}//end if(m_b_include_poly_boundary)

		}

		return r;
	}

}//end of namespace
