/*! \file GridIntersection.cpp

\brief a class for computing the intersection between the grid and a geometric object

\author <a href="http://masc.cs.gmu.edu">MASC group</a>, George Mason University
\author <a href="http://profile.usgs.gov/langevin/">Christian Langevin</a>, USGS
\bug    No known bugs.
*/

#include "GridIntersection.h"
#include "shapefile/shpReader.h"
#include <algorithm>
#include <functional>
#include <vector>
#include <string>
#include "intersecReadWrite.h"
#include "BoxPlyClip.h"
#include "PolygonChecker.h"
using namespace std;


#define DEBUG 00

namespace cusg
{

	void destroyPlylineCell(PolylineCellIntersectionInfo& info)
	{
		info.m_plyline->destroy();
	}
	void destroyPlygonCell(PolygonCellIntersectionInfo& info)
	{		
		for(GIS_polygon::iterator git =info. m_interPolygon->begin(); git != info.m_interPolygon->end(); ++git)
		{
			git->destroy();
		}
	}

	//GridIntersection

	void GridIntersection::build(grid_intersection_raw_data& data)
	{
		//check the refinement rules
		map<string,ShapeReader> shp_readers;

		//cast
		QuadTree3D* qtgrid=dynamic_cast<QuadTree3D*>(data.grid);
		ModflowGrid* mfgrid=dynamic_cast<ModflowGrid*>(data.grid); 
		ModflowGrid2D* mf2dqtgrid=dynamic_cast<ModflowGrid2D*>(data.grid);

		//as least one of the grids needs to be valid
		assert(qtgrid||qtgrid||mf2dqtgrid);

			//rotate
		ModflowGrid* rmfgrid = mfgrid;
		if(rmfgrid==NULL)
		{
			if(qtgrid!=NULL)
				rmfgrid = qtgrid->getModflowGrid();
		}
		if(rmfgrid == NULL)
		{
			cerr<<"error: abnormal grid !"<<endl; 
		}
		double cx, cy, angle;
		rmfgrid->getRotatePara(cx, cy, angle);

		//read from shapefile
		//ShapeReader shp_reader;
		m_shp_reader.destroy();
		m_shp_reader.read(data.shapefile);

		//only support
		if(data.feature_type=="line")
		{
			vector<GIS_plyline>& arcs=m_shp_reader.getArcList();
			//subdivide the qtree using the features
			int arcssize=arcs.size();
			for(int j=0;j<arcssize;j++)
			{
				int previousSegInterNum  = m_seg_interaction_info.size();

				int seg_size=arcs[j].getSize()-1; //getSize returns the number of vertices
				for(int is=0;is<seg_size;is++)
				{
					GIS_LineSeg2d line=arcs[j].getLineSegment(is);
					//rotate
					if(angle!=0)
						line.rotate(cx, cy, -angle);

					//cout<<"m_start_geo_dist["<<j<<"]["<<is<<"]="<<line.m_start_geo_dist<<","<<line.m_end_geo_dist<<endl;

					if(qtgrid!=NULL) intersect(qtgrid,line);
					else if(mfgrid!=NULL) intersect(mfgrid,line);
					else if(mf2dqtgrid!=NULL) intersect(mf2dqtgrid,line);
				}//end for i

				//merge the intersection segments that are inside the same box
				int afterSegInterNum = m_seg_interaction_info.size();
				int k = previousSegInterNum;
				while(k < afterSegInterNum)
				{
					PolylineCellIntersectionInfo plyInfo;
					plyInfo.box = m_seg_interaction_info[k].box;
					plyInfo.m_curveID = m_seg_interaction_info[k].m_curveID;
					plyInfo.m_startGeoDist = m_seg_interaction_info[k].m_startGeoDist;
					plyInfo.m_endGeoDist = m_seg_interaction_info[k].m_endGeoDist;
					plyInfo.m_len = m_seg_interaction_info[k].m_len;
					plyInfo.m_plyline = new c_plyline();
					plyInfo.m_plyline->beginPoly();
					
					Point2d& tSPnt = m_seg_interaction_info[k].m_startPnt;
					Point2d& tEPnt = m_seg_interaction_info[k].m_endPnt;
					if(angle!=0)
					{
						tSPnt.rotate(cx ,cy, angle);  tEPnt.rotate(cx, cy, angle);//rotate back
					}
					plyInfo.m_plyline->addVertex(tSPnt[0], tSPnt[1]);
					plyInfo.m_plyline->addVertex(tEPnt[0], tEPnt[1]);

					while(k < (afterSegInterNum - 1) && m_seg_interaction_info[k].m_curveID==m_seg_interaction_info[k+1].m_curveID 
						&& m_seg_interaction_info[k].box == m_seg_interaction_info[k+1].box 
						&&m_seg_interaction_info[k].m_endPnt.almost_equ(m_seg_interaction_info[k+1].m_startPnt))
					{
						tEPnt = m_seg_interaction_info[k+1].m_endPnt;
						if(angle!=0)
						{
							tEPnt.rotate(cx, cy, angle);//rotate back
						}
						plyInfo.m_plyline->addVertex(tEPnt[0], tEPnt[1]);
						plyInfo.m_endGeoDist = m_seg_interaction_info[k+1].m_endGeoDist;
						plyInfo.m_len += m_seg_interaction_info[k+1].m_len;
						k++;
					}

					plyInfo.m_plyline->endPoly();
					m_plyline_intersection_info.push_back(plyInfo);

					k++;
				}
				//arcs[j].destroy();

			}//end for j
		}//end if
		else if(data.feature_type=="polygon")
		{
			int prePlygonSize = m_ply_interaction_info.size();
			vector<GIS_polygon>& polygons=m_shp_reader.getPolyList();
			//subdivide the qtree using the features
			int psize=polygons.size();
			for(int j=0;j<psize;j++)
			{
				//if(j != 3728)
				//	continue;
				//cout<<"poly "<<j<<endl;

				if(j%500==0)
					cout<<"- Processing polygon: "<<j<<"    ......"<<endl;

				//rotate
				if(angle!=0)
					polygons[j].rotate(cx, cy, -angle);
				
				vector<GIS_polygon*> tvplygons;
				bool allInhole =checkHole(polygons[j], tvplygons);  

				//cout<<"begin intersection"<<endl;
				//int k = 0;
				//loop over all the polygons
				for(vector<GIS_polygon*>::iterator git = tvplygons.begin(); git != tvplygons.end(); ++git)
				{
					//cout<<k<<endl; k++;
					GIS_polygon gp = *(*git);

					BoxPlyClip * bpc = convertCPolygonToClip(*(*git));

					if(qtgrid!=NULL) intersect(qtgrid, gp, bpc);	//polygons[j]);
					else if(mfgrid!=NULL) intersect(mfgrid, gp,  bpc);		//polygons[j]);
					else if(mf2dqtgrid!=NULL) intersect(mf2dqtgrid, gp,  bpc);		//polygons[j]);

					bpc->destroy();
					delete bpc;
				}
				//tvplygons.pop_back();
				for(vector<GIS_polygon*>::iterator vit = tvplygons.begin(); vit != tvplygons.end(); ++vit)
				{
					GIS_polygon * gp = *vit;
					gp->destroy();
					delete gp;
				}
				//polygons[j].destroy();//delete the object
			}
			int afterSize =  m_ply_interaction_info.size();
			if(angle != 0)
			{
				for(int k = prePlygonSize; k < afterSize; k++)
				{
					m_ply_interaction_info[k].m_interPolygon->rotate(cx, cy, angle);
				}
			}

		}
		else if(data.feature_type=="point")
		{
			int prePlygonSize = m_pt_interaction_info.size();
			vector<GIS_Point2d>& pts=m_shp_reader.getPointList();
			//subdivide the qtree using the features
			int psize=pts.size();
			for(int j=0;j<psize;j++)
			{

				//rotate
				if(angle!=0)
					pts[j].rotate(cx, cy, -angle);

				if(qtgrid!=NULL) intersect(qtgrid,pts[j]);
				else if(mfgrid!=NULL) intersect(mfgrid,pts[j]);
				else if(mf2dqtgrid!=NULL) intersect(mf2dqtgrid,pts[j]);
			}
			//rotate back
			int afterSize =  this->m_pt_interaction_info.size();
			if(angle != 0)
			{
				for(int k = prePlygonSize; k < afterSize; k++)
				{
					m_pt_interaction_info[k].pnt.rotate(cx, cy, angle);
				}
			}

		}//end if
	}

	//shape file
	bool GridIntersection::saveShp(string& shpName, vector<string>& attributes)
	{
		if(shpName.size()!=0)
		{

			if(m_pt_interaction_info.size() != 0){
				writePntIntersect(this, shpName, attributes, WRITE_SHP);
			}
			else if(m_seg_interaction_info.size() != 0){
				writeLineIntersect(this, shpName, attributes,WRITE_SHP);
			}
			else if(m_ply_interaction_info.size() != 0){
				writePolyIntersect(this, shpName, attributes,WRITE_SHP);
			}
			else{
				cerr<<"! Warning: no data to write!\n";
				return true;
			}
			cout<<"- Save shape data to file: "<<shpName<<endl;
		}

		return true;
	}

	//save txt file
	bool GridIntersection::save(string& filename, vector<string>& attributes)
	{
		if(filename.size()!=0)
		{

			if(m_pt_interaction_info.size() != 0){
				writePntIntersect(this, filename, attributes, WRITE_TXT);
			}
			else if(m_seg_interaction_info.size() != 0){
				writeLineIntersect(this, filename, attributes, WRITE_TXT);
			}
			else if(m_ply_interaction_info.size() != 0){
				writePolyIntersect(this, filename, attributes, WRITE_TXT);
			}
			else{
				cerr<<"no data to write!\n";
				return true;
			}
			cout<<"- Save txt data to file: "<<filename<<endl;
		}
		else{
			cerr<<"output_file is not defined!"<<endl;
		}

		return true;
	}

	//save vtk file
	bool GridIntersection::saveVTK(string& filename, vector<string>& attributes)
	{
		//if(filename.size()!=0)
		//{
		//	if(m_pt_interaction_info.size() != 0){
		//		writePntIntersect2VTK(this, filename, attributes);
		//	}
		//	else if(m_seg_interaction_info.size() != 0){
		//		writeLineIntersect2VTK(this, filename, attributes);
		//	}
		//	else if(m_ply_interaction_info.size() != 0){
		//		writePolyIntersect2VTK(this, filename, attributes);
		//	}
		//	else{
		//		cerr<<"no data to write!\n";
		//		return true;
		//	}
		//	cout<<"- Savevtk data to file: "<<filename<<endl;
		//}
		//else{
		//	cerr<<"output_file is not defined!"<<endl;
		//}
		cerr<<"Saving intersection file to VTK is not suported yet!\n";
		return true;
	}


	//
	//
	// ModflowGrid2D
	//
	//


	bool GridIntersection::intersect(ModflowGrid2D * grid, const Point2d& pt)
	{
		Index2d id2d=grid->getIndex(pt);

		if(id2d[0]==INT_MAX || id2d[1]==INT_MAX || id2d[0] ==INT_MIN || id2d[1]==INT_MIN) return false;

		//get the box and put it in m_intersection
		int id=grid->get_nodeid(id2d[0],id2d[1]);
		Box * box=grid->nodegroup[id];
		m_intersection.push_back(box);

		////shpPntId.insert(make_pair(box, m_cur_id));

		return true;
	}


	bool GridIntersection::clip_line(ModflowGrid2D * grid, const LineSeg2d& line, LineSeg2d& clipped_line)
	{
		double cx = grid->xoffset + grid->width / 2;
		double cy = grid->yoffset + grid->height / 2;

		Box box(cx, cy, grid->width, grid->height);

		bool r = intersect(&box, line);
		if (r == false) return false; //there is no intersection...

		r = getIntersection(&box, line, clipped_line);
		if (r == false) return false; //there is no intersection...

		if (box.in(clipped_line.pt[0][0], clipped_line.pt[0][1]) == false)
		{
			//move clipped_line.pt[0] toward clipped_line.pt[1] until the point is inside the box
			Point2d mid((clipped_line.pt[0][0] + clipped_line.pt[1][0]) / 2, (clipped_line.pt[0][1] + clipped_line.pt[1][1]) / 2);
			bracket_to_box_boundary(box, clipped_line.pt[0], mid);
			clipped_line.pt[0]=mid;
		}

		if (box.in(clipped_line.pt[1][0], clipped_line.pt[1][1]) == false)
		{
			//move clipped_line.pt[1] toward clipped_line.pt[0]
			Point2d mid((clipped_line.pt[0][0] + clipped_line.pt[1][0]) / 2, (clipped_line.pt[0][1] + clipped_line.pt[1][1]) / 2);
			bracket_to_box_boundary(box, clipped_line.pt[1], mid);
			clipped_line.pt[1] = mid;
		}

		return true;
	}

	bool GridIntersection::clip_line(QuadTree3D * grid, const LineSeg2d& line, LineSeg2d& clipped_line)
	{
		return clip_line((ModflowGrid2D *)grid->getModflowGrid(),line,clipped_line);
	}

	void GridIntersection::bracket_to_box_boundary(Box & box, Point2d& out, Point2d& in) const
	{
		//make sure that in point is in and out point is out of box
		if (box.in(in[0], in[1]) == false)
		{
			cout << "holy" << endl;
		}

		assert(box.in(in[0], in[1]));
		assert(box.in(out[0], out[1]) == false);
		//
		while ( true )
		{
			Point2d mid((in[0] + out[0]) / 2, (in[1] + out[1]) / 2);
			if (Equal(mid.get(), in.get())) break;
			if (Equal(mid.get(), out.get())) break;

			if (box.in(mid[0], mid[1]))
			{
				in = mid;
			}
			else out = mid;
		}//end while
	}

	//
	//
	// ModflowGrid
	//
	//
	bool GridIntersection::intersect(ModflowGrid * grid, const Point2d& pt)	//, bool noBound)
	{
        Index2d id2d=grid->getIndex(pt);
        if(id2d[0]==INT_MAX || id2d[1]==INT_MAX|| id2d[0] ==INT_MIN || id2d[1]==INT_MIN) return false;

        //get the box and put it in m_intersection
        int id=grid->get_nodeid(this->m_layer,id2d[0],id2d[1]);

        Box * box=grid->nodegroup[id];
        m_intersection.push_back(box);

        return true;
	}

	//
	//
	//QuadTree3d
	//
	//


	bool GridIntersection::intersect(QuadTree3D * grid, const Point2d& pt) //, bool noBoud)
	{
		//vector<Box*> tmp;
		//tmp.swap(m_intersection);
		ModflowGrid * mf=grid->getModflowGrid();

		//cout<<"Check: mf->nrow="<<mf->nrow<<" mf->ncol="<<mf->ncol<<endl;

		Index2d id2d=mf->getIndex(pt);

		//cout<<"pt="<<pt<<endl;
		//cout<<"id2d="<<id2d[0]<<","<<id2d[1]<<endl;
		if(id2d[0]==INT_MAX || id2d[1]==INT_MAX|| id2d[0] ==INT_MIN || id2d[1]==INT_MIN) return false;

		//get the box and put it in m_intersection
		int id=mf->get_nodeid(m_layer,id2d[0],id2d[1]);
		Box * box=grid->nodegroup[id];

#if DEBUG
		bool inside=box->in(pt[0],pt[1]);
		if(inside==false)
		{
			Index2d id2d = mf->getIndex(pt);
			mf->get_nodeid(m_layer, id2d[0], id2d[1]);
			box->in(pt[0], pt[1]);
			cout<<"id="<<id<<" m_layer="<<m_layer<<endl;
			cout<<"modflowGrid="<<mf<<endl;
			cout<<"qtree_grid="<<grid<<endl;
			cout<<"Point "<<pt<<" is not inside "<<*box<<" index="<<id2d<<endl;
		}
#endif

		box=box->find(pt[0],pt[1]); //find a leave that contains pt

		//make sure box is valid
		assert(box);

		//remember
		m_intersection.push_back(box);

		////shpPntId.insert(make_pair(box, m_cur_id));//add by guilin

		//remember the boxes that are exactly on the boundary of the polygon
		if(box->onBound(pt[0], pt[1]))
			this->m_onBound.push_back(box);

		return true;
	}


	//Box
	bool GridIntersection::intersect(Box * b, const Point2d& pt)
	{
		return b->in(pt[0],pt[1]);
	}

	//check intersection between a box and a line
	//TODO: this can be improved
	bool GridIntersection::intersect(Box * b, const LineSeg2d& line)
	{
		//check if the end points are in
		if(b->in(line.pt[0][0],line.pt[0][1])) return true;
		if(b->in(line.pt[1][0],line.pt[1][1])) return true;

		//create each edge of b
		const Point2d ul(b->x-b->dx/2,b->y+b->dy/2); //upper left corner
		const Point2d ur(b->x+b->dx/2,b->y+b->dy/2); //upper right corner
		const Point2d lr(b->x+b->dx/2,b->y-b->dy/2); //lower right corner
		const Point2d ll(b->x-b->dx/2,b->y-b->dy/2); //lower left corner

		const double * s=line.pt[0].get();
		const double * t=line.pt[1].get();

		//now, check intersections between the line segment and the edges
		if( SegSegInt(ul.get(), ur.get(), s, t) ) return true;
		if( SegSegInt(ur.get(), lr.get(), s, t) ) return true;
		if( SegSegInt(lr.get(), ll.get(), s, t) ) return true;
		if( SegSegInt(ll.get(), ul.get(), s, t) ) return true;

		//JML: I am not sure if the code below is really necessary.
		if (AreaSign(ul.get(), ur.get(), s) == 0) if (Between(ul.get(), ur.get(), s)) return true;
		if (AreaSign(ur.get(), lr.get(), s) == 0) if (Between(ur.get(), lr.get(), s)) return true;
		if (AreaSign(lr.get(), ll.get(), s) == 0) if (Between(lr.get(), ll.get(), s)) return true;
		if (AreaSign(ll.get(), ul.get(), s) == 0) if (Between(ll.get(), ul.get(), s)) return true;
		
		if (AreaSign(ul.get(), ur.get(), t) == 0) if (Between(ul.get(), ur.get(), t)) return true;
		if (AreaSign(ur.get(), lr.get(), t) == 0) if (Between(ur.get(), lr.get(), t)) return true;
		if (AreaSign(lr.get(), ll.get(), t) == 0) if (Between(lr.get(), ll.get(), t)) return true;
		if (AreaSign(ll.get(), ul.get(), t) == 0) if (Between(ll.get(), ul.get(), t)) return true;

		//otherwise
		return false;
	}

	void GridIntersection::buildIntersectionInfo(Box * box, const GIS_polygon& ply, PolygonCellIntersectionInfo& info)
	{
		info.box=box;
		info.m_area=box->dx*box->dy*2;
		info.m_polyID = ply.m_id;
		info.m_interPolygon=new GIS_polygon();

		Point2d UL, UR, LR, LL;
		box->getCorners(UL,UR,LR,LL);

		c_ply cply(c_ply::POUT);
		cply.beginPoly();
		cply.addVertex(UL[0],UL[1]);
		cply.addVertex(UR[0],UR[1]);
		cply.addVertex(LR[0],LR[1]);
		cply.addVertex(LL[0],LL[1]);
		cply.endPoly();
		info.m_interPolygon->push_back(cply);
	}

	void GridIntersection::buildIntersectionInfo
	(Box * box, const GIS_polygon& ply,  vector<PolygonCellIntersectionInfo>& infoResult, BoxPlyClip* bpc)
	{
		vector<GIS_polygon*> resPlygons;
		intersect_polygon(box, ply, resPlygons,bpc);
		for(vector<GIS_polygon*>::iterator pit = resPlygons.begin(); pit != resPlygons.end(); ++pit)
		{
			GIS_polygon* gp = *pit;
			PolygonCellIntersectionInfo info;
			info.m_interPolygon = gp;
			info.box = box;
			info.m_polyID = ply.m_id;
			info.m_area = info.m_interPolygon->getArea();

#if DEBUGTEST
			info.m_polyID = gp->m_id;
#endif

			infoResult.push_back(info);
		}

	}

	SegmentCellIntersectionInfo GridIntersection::buildIntersectionInfo(Box * box, const GIS_LineSeg2d& line)
	{
		SegmentCellIntersectionInfo info;
		info.box=box;
		info.m_segID=line.m_seg_id;
		info.m_curveID=line.m_curve_id;

		LineSeg2d intseg;
		bool r=getIntersection(box, (LineSeg2d)line, intseg);

		if (r == false)
		{
			cerr << "! Warning: GridIntersection::buildIntersectionInfo(Box,const GIS_LineSeg2d): Cannot find intersection information" << endl;
			return info;
		}

		info.m_len=(intseg.pt[0]-intseg.pt[1]).norm();
		info.m_startGeoDist=(intseg.pt[0]-line.pt[0]).norm()+line.m_start_geo_dist;
		info.m_endGeoDist=(intseg.pt[1]-line.pt[0]).norm()+line.m_start_geo_dist;

		if(info.m_startGeoDist>info.m_endGeoDist)
		{
				swap(info.m_endGeoDist,info.m_startGeoDist);
				info.m_startPnt = intseg.pt[1];
				info.m_endPnt = intseg.pt[0];
		}
		else
		{
			info.m_startPnt = intseg.pt[0];
			info.m_endPnt = intseg.pt[1];
		}

		return info;
	}

	PointCellIntersectionInfo GridIntersection::buildIntersectionInfo(Box * box, const GIS_Point2d& pt)
	{
		PointCellIntersectionInfo info;
		info.box=box;
		info.m_ptID=pt.m_id;
		info.pnt = pt;
		return info;
	}


	bool GridIntersection::getIntersection(Box * b, const LineSeg2d& line, LineSeg2d& intseg)
	{

		//check if the end points are in
		bool vinbox[2]={false};
		if(b->in(line.pt[0][0],line.pt[0][1])) vinbox[0]=true;
		if(b->in(line.pt[1][0],line.pt[1][1])) vinbox[1]=true;

		if(vinbox[0] && vinbox[1]){
			intseg=line;
			return true;
		}

		//create each edge of b
		Point2d ul(b->x-b->dx/2,b->y+b->dy/2); //upper left corner
		Point2d ur(b->x+b->dx/2,b->y+b->dy/2); //upper right corner
		Point2d lr(b->x+b->dx/2,b->y-b->dy/2); //lower right corner
		Point2d ll(b->x-b->dx/2,b->y-b->dy/2); //lower left corner

		//now, check intersections between the line segment and the edges
		char r[4];
		Point2d p[4];

		r[0]=SegSegInt(ul.get(), ur.get(), line.pt[0].get(),line.pt[1].get(),p[0].get());
		r[1]=SegSegInt(ur.get(), lr.get(), line.pt[0].get(),line.pt[1].get(),p[1].get());
		r[2]=SegSegInt(lr.get(), ll.get(), line.pt[0].get(),line.pt[1].get(),p[2].get());
		r[3]=SegSegInt(ll.get(), ul.get(), line.pt[0].get(),line.pt[1].get(),p[3].get());

		//no intersection....
		if ((r[0] || r[1] || r[2] || r[3]) == false)
			return false;

		if(vinbox[0] || vinbox[1])
		{
			float max_len=0;
			for(int i=0;i<4;i++)
			{
				if(r[i]=='0') continue;
				if(r[i]=='e')
				{
					//compute the union between the edge of box and the line
					Point2d seg[2];
					const double * u[2]={NULL,NULL};

					switch(i){
					case 0: seg[0]=ul; seg[1]=ur; break;
					case 1: seg[0]=ur; seg[1]=lr; break;
					case 2: seg[0]=lr; seg[1]=ll; break;
					case 3: seg[0]=ll; seg[1]=ul; break;
					}

					bool r=Union(seg[0].get(),seg[1].get(),line.pt[0].get(),line.pt[1].get(),u);
					if(r)
					{
						intseg.pt[0].set(u[0][0],u[0][1]);
						intseg.pt[1].set(u[1][0],u[1][1]);
					}

					return r;
				}
				//found the intersection...
				const Point2d& pt1=p[i];
				const Point2d& pt2=(vinbox[0])?line.pt[0]:line.pt[1];

				double len=(pt1-pt2).normsqr();
				if(len>max_len)
				{
					intseg.pt[0]=pt1;
					intseg.pt[1]=pt2;
					max_len=len;
				}

			}
			//done

			return true;
		}

		//count
		float max_len=0;
		for(int i=0;i<4;i++)
		{
			if(r[i]=='0') continue;

			if(r[i]=='e') continue; //this can't happen, it should be handled in the previous loop

			for(int j=i;j<4;j++)
			{
				if(r[j]=='0') continue;
				if(r[j]=='e') continue;

				//found the intersection...
				double len=(p[i]-p[j]).normsqr();
				if(len>max_len)
				{
					intseg.pt[0]=p[i];
					intseg.pt[1]=p[j];
					max_len=len;
				}
			}
		}

		return true;
	}

	double GridIntersection::intersect_len(Box * b, const LineSeg2d& line)
	{
		//check if the end points are in
		LineSeg2d intseg;
		bool r = getIntersection(b, line, intseg);
		if (r)
			return (intseg.pt[0]-intseg.pt[1]).norm();
		return 0; //no intersection...
	}

	//c_ply GridIntersection::getIntersection(Box * b, const c_ply& ply)
	//{
	//	c_ply intply(c_ply::POUT);
	//	return intply;
	//}

	//bool gridintersection::intersect(box * b, const c_ply& ply)
	//{

	//    return true;
	//}

	//
	// check if a point is inside a ply
	// ray shooting is slow. Can be improved.
	//
	bool GridIntersection::intersect(const Point2d& pt, const c_ply& ply)
	{

		static int totoal_call=0;
		static int totoal_failed_call=0;
		totoal_call++;

		int max_trial=200;

		Vector2d offset( ply.getHead()->getPos().get() );
		float d=ply.getRadius()*4;
		Point2d near=pt-offset;
		Point2d far=near;

		while(max_trial>0)
		{
			//get a random point
			Vector2d vec(drand48(),drand48());
			if(drand48()<0.5) vec[0]=-vec[0];
			if(drand48()<0.5) vec[1]=-vec[1];
			far=pt-offset+vec.normalize()*d;

			//
			ply_vertex * ptr=ply.getHead();
			bool degenerate=false;
			int int_count=0;

			do{
				ply_vertex * next=ptr->getNext();
				double x[2];
				Point2d c=ptr->getPos()-offset;
				Point2d d=next->getPos()-offset;

				char r=SegSegInt(near.get(),far.get(), c.get(), d.get() , x);

				if(r=='1') int_count++;
				else if(r!='0')
				{
					degenerate=true;
					break;
				}

				ptr=next;
			}
			while(ptr!=ply.getHead());

			//cout<<"int_count="<<int_count<<endl;

			if(degenerate){  max_trial--; continue; }
			return (int_count%2)==1;
		}

		totoal_failed_call++;

		//cout<<"Max iteration reached.. near="<<near<<" far="<<far<<" ratio="<<totoal_failed_call<<"/"<<totoal_call<<endl;

		return false;
	}

	//
	// check if a point is inside a polygon
	// ray shooting is slow. Can be improved.
	//
	bool GridIntersection::intersect(const Point2d& pt, const c_polygon& pgon)
	{

		static int totoal_call=0;
		static int totoal_failed_call=0;
		totoal_call++;

		int max_trial=200;

		Vector2d offset( pgon.front().getHead()->getPos().get() );
		const double * pgon_box=pgon.getBBox();
		float d=max(pgon_box[1]-pgon_box[0], pgon_box[3]-pgon_box[2])*4;

		Point2d near=pt-offset;
		Point2d far=near;

		while(max_trial>0)
		{
			//get a random point
			Vector2d vec(drand48(),drand48());
			if(drand48()<0.5) vec[0]=-vec[0];
			if(drand48()<0.5) vec[1]=-vec[1];
			far=pt-offset+vec.normalize()*d;

			//count the number of intersection
			bool degenerate=false; //flag raised if degenerate intersection is found.
			int int_count=0;       //# of intersections

			//loop through each ply
			for(c_polygon::const_iterator i=pgon.begin();i!=pgon.end();i++)
			{
				const c_ply& ply=*i;
				ply_vertex * ptr=ply.getHead();

				do{
					ply_vertex * next=ptr->getNext();
					double x[2];
					Point2d c=ptr->getPos()-offset;
					Point2d d=next->getPos()-offset;

					char r=SegSegInt(near.get(),far.get(), c.get(), d.get() , x);

					if(r=='1') int_count++;
					else if(r!='0')
					{
						degenerate=true;
						break;
					}

					ptr=next;
				}
				while(ptr!=ply.getHead());

				if(degenerate) break; //something is wrong
			}

			if(degenerate){  max_trial--; continue; }

#if DEBUG
			cout<<"- GridIntersection::intersect: int_count="<<int_count<<endl;
#endif

			return (int_count%2)==1;
		}

		totoal_failed_call++;

#if DEBUG
		cerr<<"! Warning: Max iteration reached.. near="<<near<<" far="<<far<<" ratio="<<totoal_failed_call<<"/"<<totoal_call<<endl;
#endif

		return false;
	}

	//void GridIntersection::intersect_polygon(Box* b, const GIS_polygon& plygon, vector<GIS_polygon*>& resultPlygons)
	void GridIntersection::intersect_polygon(Box* b, const GIS_polygon& plygon,  vector<GIS_polygon*>& resultPlygons,BoxPlyClip* bpcInit)
	{	
		//if(b->number != 16443)
		//	return;

		//ofstream ofile("C:\\Guilin\\matlab files\\cor.txt");
		//GIS_polygon gpn;
		//gpn.copy(plygon);
		//int vsize = gpn.getSize();
		//for(int i = 0; i< vsize; i++)
		//{
		//	ofile<<gpn[i]->getPos()[0]<<"\t"<<gpn[i]->getPos()[1]<<endl;
		//}
		//ofile.close();
		//cout<<"num "<<b->number<<endl;


		double xmin = min(b->x - b->dx/2, b->x + b->dx/2);
		double xmax = max(b->x - b->dx/2, b->x + b->dx/2);
		double ymin = min(b->y - b->dy/2, b->y + b->dy/2);
		double ymax = max(b->y - b->dy/2, b->y + b->dy/2);
		//cout<<"current id : "<<b->id<<endl;

		/*****************************************************************/
		BoxPlyClip* bpc = bpcInit;
		if(bpcInit==NULL)
		{
			bpc = convertCPolygonToClip(plygon);
		}

		list<BoxPlyClip*> cres;
		bpc->clip(xmin, xmax, ymin, ymax, cres);

		for (list<BoxPlyClip*>::iterator cit = cres.begin(); cit != cres.end(); ++cit)
		{
			GIS_polygon* nplygon = new GIS_polygon();
			convertClipToPolygon(*cit, nplygon);
			resultPlygons.push_back(nplygon);

#if DEBUGTEST
			nplygon->m_id = plygon.m_id;
#endif
			
			(*cit)->destroy();
			delete *cit;
		}

		bpc->resetStatus();


		//the holes 
		for(list<VNode*>::iterator vit = bpc->plys.begin(); vit != bpc->plys.end(); ++vit)
		{
			VNode* vn = *vit;
			if(vn == bpc->first)
				continue;
			VNode* nh = copyOnePly(vn, false);
			BoxPlyClip* tbp = new BoxPlyClip();
			tbp->first = nh;
			tbp->plys.push_back(tbp->first);

			list<BoxPlyClip*> hcres;
			tbp->clip(xmin, xmax, ymin, ymax, hcres);
			for (list<BoxPlyClip*>::iterator cit = hcres.begin(); cit != hcres.end(); ++cit)
			{
				GIS_polygon* nplygon = new GIS_polygon();
				convertClipToPolygon(*cit, nplygon);
				nplygon->front().set(c_ply::POUT, nplygon->front().getHead());
				resultPlygons.push_back(nplygon);

				(*cit)->destroy();
				delete *cit;
			}
			tbp->destroy();

			delete tbp;
		}


		if(bpcInit==NULL)
		{
			bpc->destroy();
			delete bpc;
		}
		else
		{
			bpc->resetStatus();
		}

		//bpc->destroy();
		//delete bpc;
	}

	//check if b1 and b2 are neighbors that sharing a same edge/face, return -1 if nor neighbors
	int GridIntersection::checkNeighbor(Box* b1, Box* b2)
	{
		for(int i = 0; i < 4; i++)
		{
			list<Box*> bl;
			b1->getNeighbors(bl, i);
			if(std::find(bl.begin(), bl.end(), b2) != bl.end())
				return i;
		}

		return -1;
	}

	//check box edge collinearly intersect with line
	bool GridIntersection::checkEdgeCollinear(Box* box, int dir, const LineSeg2d& line)
	{
		double a[2], b[2];
		double c[2]={line.pt[0][0], line.pt[0][1]};
		double d[2]={line.pt[1][0], line.pt[1][1]};
		double p[2];
		if(dir == 0)
		{
			a[0] = box->x - box->dx/2;		b[0] = box->x + box->dx/2;
			a[1] = box->y + box->dy/2;		b[1] = box->y + box->dy/2;
		}
		else if(dir ==1)
		{
			a[0] = box->x + box->dx/2;		b[0] = box->x + box->dx/2;
			a[1] = box->y + box->dy/2;		b[1] = box->y - box->dy/2;
		}
		else if(dir ==2)
		{
			a[0] = box->x + box->dx/2;		b[0] = box->x - box->dx/2;
			a[1] = box->y - box->dy/2;		b[1] = box->y - box->dy/2;
		}
		else
		{
			assert(dir==3);
			a[0] = box->x - box->dx/2;		b[0] = box->x - box->dx/2;
			a[1] = box->y + box->dy/2;		b[1] = box->y - box->dy/2;
		}
		//check collinear
		if(SegSegInt(a, b, c, d, p) == 'e')
			return true;

		return false;
	}

	
	///collect all cells in the grid from b 
	void GridIntersection::floodGrid(Box * sb, set<Box*>& visited, list<Box*>& flooded)
	{
		list<Box*> open;
		open.push_back(sb);
	
		while(open.empty()==false) //not empty
		{
			Box* b=open.front();
			open.pop_front();
			flooded.push_back(b);

			//collect neighboring cells
			list<Box*> neighbors[4];

			for (int i = 0; i < 4; i++)
			{
				b->getNeighbors(neighbors[i], i);
			}

			//get all unvisited cells
			for(int i=0;i<4;i++)
			{
				for(list<Box*>::iterator n=neighbors[i].begin(); n!=neighbors[i].end(); n++)
				{
					Box * nb=*n;
					if(visited.find(nb)!=visited.end()) continue; //visited
					visited.insert(nb);
					open.push_back(nb);
				}//end for n
			}//end for i
		}//end while

	}//end GridIntersection::floodGrid
}

