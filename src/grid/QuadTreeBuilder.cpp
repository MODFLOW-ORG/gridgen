/*! \file QuadTreeBuilder.cpp

\brief QuadTree Builder

\author <a href="masc.cs.gmu.edu/">MASC group</a>, George Mason University,
        <a href="profile.usgs.gov/langevin/">Christian Langevin</a>, USGS
\bug    No known bugs.

 */

//   $Id: $
#include <locale>
#include <string>
#include "QuadTreeBuilder.h"
#include "shapefile/shpReader.h"
#include "GridDomain.h"

namespace cusg
{

inline string tolower(const string& s)
{
    string lower;
    for(string::const_iterator i=s.begin();i!=s.end();i++)
    {
        char c=::tolower(*i);
        lower.push_back(c);
    }
    return lower;
}

inline int getField(string fname, ShapeReader& shpReader)
{
	int id = 0;
	vector<ShapeIntField>& intFlds = shpReader.getIntFields();
	for(vector<ShapeIntField>::iterator sit = intFlds.begin(); sit != intFlds.end(); ++sit, ++id)
	{
		string tstr = sit->name;
		//std::transform(sit->name.begin(), sit->name.end(), tstr.begin(), ::tolower);
		if(fname == tstr.c_str())
			return id;
	}

	return -1;
}



void Quadtree_Builder::build(QuadTree3D * qtree, quadtree_builder_raw_data& rawdata)
{
    //rotate qtree
	ModflowGrid* rmfgrid = qtree->getModflowGrid();
	double cx, cy, angle;
	rmfgrid->getRotatePara(cx, cy, angle);

	//set z-value interploation methods
	for(int i=0;i<rawdata.modflow_data->nlay;i++)
	{
		string operation_top = tolower(rawdata.top[i].first);

		if(operation_top=="replicate") qtree->setTopOperators(i,QuadTree3D::Z_VALUE_REPLICATE_OPEARTOR);
		else if(operation_top=="interpolate") qtree->setTopOperators(i,QuadTree3D::Z_VALUE_LINEAR_INTERPOLATION_OPEARTOR);
		else if(operation_top=="asciigrid") //ASCII file
		{
			qtree->setTopOperators(i,QuadTree3D::Z_VALUE_LINEAR_INTERPOLATION_OPEARTOR);
			qtree->setTopSource(i,rawdata.top[i].second);

			if(rawdata.top_area_weighted.find(i) != rawdata.top_area_weighted.end())
			{
				qtree->setTopAreaWeighted(i, rawdata.top_area_weighted[i]);
			}
		}
		else
		{
			cerr<<"! Error: Unknown top "<<i<<" elevation operator: "<<rawdata.top[i].first<<endl;
			exit(1);
		}

		string operation_bot = tolower(rawdata.bot[i].first);
		if(operation_bot=="replicate") qtree->setBotOperators(i,QuadTree3D::Z_VALUE_REPLICATE_OPEARTOR);
		else if(operation_bot=="interpolate") qtree->setBotOperators(i,QuadTree3D::Z_VALUE_LINEAR_INTERPOLATION_OPEARTOR);
		else if(operation_bot=="asciigrid") //ASCII file
		{
			qtree->setBotOperators(i,QuadTree3D::Z_VALUE_LINEAR_INTERPOLATION_OPEARTOR);
			qtree->setBotSource(i,rawdata.bot[i].second);

			if(rawdata.bot_area_weighted.find(i) != rawdata.bot_area_weighted.end())
			{
				qtree->setBotAreaWeighted(i, rawdata.bot_area_weighted[i]);
			}
		}
		else
		{
			cerr<<"! Error: Unknown bottom "<<i<<" elevation operator: ("<<rawdata.bot[i].first<<")"<<endl;
			exit(1);
		}
	}

    //check the refinement rules
    map<string,ShapeReader> shp_readers;
    for(int i=0;i<rawdata.modflow_data->nlay;i++)
    {
		cout<<"  * Refine grid at layer "<<i+1<<endl;

        list<refinement_raw_data*>& rF=rawdata.refine_by_features[i];
        for(list<refinement_raw_data*>::iterator f=rF.begin();f!=rF.end();f++)
        {
            refinement_raw_data* rf=*f;
            ShapeReader& shp_reader=shp_readers[rf->shapefile];
            shp_reader.read(rf->shapefile);

			//prepare for attribute based refinement and store the specified using fid 
			int fid = -1; //index of refinement_level
			bool bByAttribute  = (rf->refinement_level_by_attribute.empty()==false);
			vector<ShapeIntField>& intFlds = shp_reader.getIntFields();
			
			if(bByAttribute)
			{
				//cerr<<"refine by attribute"<<endl;
				fid = getField(rf->refinement_level_by_attribute, shp_reader);
				if(fid == -1)
				{
					cerr<<"! Error: The field for refinement_level_by_attribute : "
						<<rf->refinement_level_by_attribute<<"is not found in shape file: "
						<<rf->shapefile<<endl;
					continue;
				}
			}

            // refine using the given feature (rf)
			string rftype=tolower(rf->featuretype);
			if(rftype=="line")
            {
                vector<GIS_plyline>& arcs=shp_reader.getArcList();
                //subdivide the qtree using the features
                int arcssize=arcs.size();
                for(int j=0;j<arcssize;j++)
				{
					//rotate
					if(angle!=0) arcs[j].rotate(cx, cy, -angle);

					int rf_level=(bByAttribute)?intFlds[fid].vals[j]:rf->refinement_level;
					qtree->refine( (c_plyline&)arcs[j], rf_level, i);
					
					//rotate back
					if(angle!=0) arcs[j].rotate(cx, cy, angle);
					//arcs[j].destroy();
				}
            }//end if
            else if(rftype=="polygon")
            {
                vector<GIS_polygon>& polygons=shp_reader.getPolyList();
                //subdivide the qtree using the features
                int psize=polygons.size();

                for(int j=0;j<psize;j++)
				{
					//make sure that bounding box is computed
			        polygons[j].buildBoxAndCenter();

					//rotate
					if(angle!=0) polygons[j].rotate(cx, cy, -angle);

					int rf_level=(bByAttribute)?intFlds[fid].vals[j]:rf->refinement_level;
					qtree->refine(polygons[j], rf_level,i); //i is the level index

					//rotate back
					if(angle!=0) polygons[j].rotate(cx, cy, angle);
					//polygons[j].destroy();
				}
            }
            else if(rftype=="point")
            {
                vector<GIS_Point2d>& pts=shp_reader.getPointList();
                //subdivide the qtree using the features
                int psize=pts.size();
                for(int j=0;j<psize;j++)
				{
					//rotate
					if(angle!=0) pts[j].rotate(cx, cy, -angle);

					int rf_level=(bByAttribute)?intFlds[fid].vals[j]:rf->refinement_level;
					qtree->refine(pts[j], rf_level, i);
					
					//rotate back
					if(angle!=0) pts[j].rotate(cx, cy, angle);
				}
            }
			else
			{
				cerr<<" ! Error: Unknown feature type: "<<rf->featuretype<<endl;
				exit(1);
			}
        }//end for f

    }//end for i

    if(rawdata.smoothing)
    {
		cout<<"  * Smooth grid: horizontal="<<rawdata.smoothing_level_horizontal<<" vertical="<<rawdata.smoothing_level_vertical<<endl;
		qtree->smooth_refinement(rawdata.smoothing_level_horizontal,rawdata.smoothing_level_vertical);
    }

    //if active domain is specified
    GridDomain gd;
    for(int i=0;i<rawdata.modflow_data->nlay;i++)
    {
        list<active_domain_raw_data*>& rF=rawdata.active_domain_data[i];
        if(rF.empty()) continue;

		cout<<"  * Determine the active domain for layer "<<i+1<<endl;

        qtree->nodegroup.disable_nodes_in_layer(i);
        for(list<active_domain_raw_data*>::iterator f=rF.begin();f!=rF.end();f++)
        {
            gd.build(qtree, i, *f);
        }
    }

	//
	// ask the grid to renumber nodes
	// node number are used to retrieve the leaf nodoes
	//
    qtree->set_one_based_node_numbering(rawdata.one_based_node_numbering);
	qtree->number_nodes();
	qtree->refresh();

	//
	// build surface elevation
	//
	cout<<"  * Build surface elevation"<<endl;
	double * top_elev=qtree->get_top();
	double * bot_elev=qtree->get_bot();
	qtree->recompute_z_of_cells();

	if(rawdata.auto_alignment)
	{
		//perform auto alignment here
		qtree->align_elev();
	}

	// count number of active cells with negative thickness
	analyze_elevation(qtree,top_elev, bot_elev);

}// void Quadtree_Builder::build

void Quadtree_Builder::analyze_elevation(QuadTree3D * qtree, double * top_elev, double * bot_elev)
{
	int nodesize=qtree->nodes();
	int neg_thinkness_count=0;
	int max_neg_thinkness_node=-1;
	double max_neg_thinkness=-FLT_MAX;

	//for each leaf node
	for(int nodenumber=0;nodenumber<nodesize;nodenumber++)
	{
		Box * box=qtree->get_nodeobj(nodenumber);
		assert(box->isLeaf); //make sure this node is a leaf node...
		if( box->active==false ) continue; //disabled...

		//check to top and bottom
		if(top_elev[nodenumber]<bot_elev[nodenumber]) //negative thinkness..
		{
			double negz=bot_elev[nodenumber]-top_elev[nodenumber];
			neg_thinkness_count++;
			if(negz>max_neg_thinkness) max_neg_thinkness=negz;
			max_neg_thinkness_node=nodenumber;
		}
	}//end for nodenumber

	if(neg_thinkness_count>0)
	{
		cerr<<"! Warning: There are "<<neg_thinkness_count<<" cells with negative thinkness. Node #"
			<<max_neg_thinkness_node<<" has the largest negative thinkness="<<-max_neg_thinkness<<endl;
	}
	//

}
//end Quadtree_Builder::analyze_elevation()


} //namespace cusg

