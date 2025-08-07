
#include "def_struct.h"
#include "def_parser.h"
#include "def_exporter.h"
#include "ModflowGrid.h"
#include "QuadTree3D.h"
#include "shapefile/shpReader.h"

#include "GridIntersection.h"
#include "GridToShp.h"
#include "GridToVTK.h"

#include "vert_pass.h"
#include "GridToVTK.h"

#if GL_GUI
#include "glDraw.h"
void fillinFeatures(cusgGridGUI& gui, quadtree_builder_raw_data * qtree_data);
#endif

using namespace cusg;

#if 1

#define GRIDGEN_VERSION "1.0"
#define GRIDGEN_SUB_VERSION "02"
#define GRIDGEN_DATE "January 6, 2017"
#define GRIDGEN_COPYRIGHT "This program is public domain and is released on the condition that neither the U.S. Geological Survey nor the United States Government may be held liable for any damages resulting from their authorized or unauthorized use."

int main(int argc, char ** argv)
{
	cout<<"GRIDGEN Version "<<GRIDGEN_VERSION<<"."<<GRIDGEN_SUB_VERSION<<" "<<GRIDGEN_DATE<<"\n";
	cout<<"A program for generating unstructured grids.\n";
	cout<<"\n"<<GRIDGEN_COPYRIGHT<<"\n"<<endl;
	
	if(argc<3)
	{
		cerr<<"usage: "<<argv[0]
#if GL_GUI	
		    <<" [-g]"
#endif		    
		    <<" action *.dfn"<<endl;
		return 1;
	}
	

	string action;
	string dfn_file;

	bool gl_display = false;

    for(int i=1;i<argc;i++){
        if( strcmp(argv[i],"-g")==0 ) gl_display=true;
        else if(action.empty()) action=argv[i];
        else dfn_file=argv[i];
    }

	DefinitionParser def_parser;
	bool parse_suc=def_parser.parse(dfn_file);
	
	block_raw_data * data=def_parser.getBlockData(action);

	if(data==NULL)
	{
		cerr<<"! Error: There is no block named: "<<action<<endl;
		return 1;
	}

	bool rinit=data->initialize( def_parser.getBlockData() );
	if(rinit==false)
	{
		cerr<<"! Error: Failed to initialize block: "<<action<<endl;
		return 1;
	}


	Grid * grid=NULL;

	//try to case the data to various blocks
	grid_raw_data * grid_data=dynamic_cast<grid_raw_data*>(data);
	quadtree_builder_raw_data * qtree_builder_data=dynamic_cast<quadtree_builder_raw_data*>(data);
    gridtoshape_raw_data * grid2shp = dynamic_cast<gridtoshape_raw_data*>(data);
    grid_to_usgdata_raw_data * grid2usg = dynamic_cast<grid_to_usgdata_raw_data*>(data);
    grid_intersection_raw_data * gridintersection = dynamic_cast<grid_intersection_raw_data*>(data);
	gridtovtk_raw_data* grid2vtk = dynamic_cast<gridtovtk_raw_data*>(data);

	//check if this is a type of grid
    if(grid_data!=NULL)
    {
        grid=grid_data->build();
        if(grid==NULL)
        {
            cerr<<"! Error: Failed to build block: "<<action<<endl;
            return 1;
        }
    }
    else if(grid2shp!=NULL)
    {
        grid=grid2shp->grid;
	}
    else if(grid2usg!=NULL)
    {
        grid=grid2usg->grid;
    }
    else if(gridintersection!=NULL)
    {
        grid=gridintersection->grid;
    }
	else if(grid2vtk != NULL)
	{
		grid = grid2vtk->grid;
	}


	QuadTree3D * qtree_grid=dynamic_cast<QuadTree3D*>(grid);
	ModflowGrid * modflow_grid=dynamic_cast<ModflowGrid*>(grid);
	
	//check if save files
	if(qtree_builder_data!=NULL)
	{
        DefinitionExporter def_exporter;
        const string& outputfile=qtree_builder_data->grid_definition_file;

        def_exporter.save(outputfile,*qtree_builder_data);
	}
	else if(grid2shp!=NULL)
	{
	    cout<<"- Saving grid: "<<grid2shp->grid_name<<" to shapefile: "<<grid2shp->shapefile<<endl;
		//save to shapefile...
	    if(qtree_grid!=NULL)
	    {
			writeQuadtree2Shp(qtree_grid, grid2shp->shapefile, grid2shp->feature_type, grid2shp->concide);
	    }
	    else if(modflow_grid!=NULL)
	    {
			writeModflowGrid2Shp(modflow_grid, grid2shp->shapefile, grid2shp->feature_type, grid2shp->without_inactive, grid2shp->one_based_node_numbering, grid2shp->concide);
	    }
	    else{
	    	cerr<<"! Error: Unknown Grid type"<<endl;
	    }	    
	}
	else if(grid2vtk!=NULL)
	{
		cout<<"- Saving grid: "<<grid2vtk->grid_name<<" to vtkfile: "<<grid2vtk->vtkfile<<endl;
		//save to shapefile...
	    if(qtree_grid!=NULL)
	    {
			if(grid2vtk->b_share_vertex)
			{
				writeQuadtreeGrid2VTK_shareV(qtree_grid, grid2vtk->vtkfile);
			}
			else
			{
				writeQuadtreeGrid2VTK(qtree_grid, grid2vtk->vtkfile);
			}
	    }
	    else if(modflow_grid!=NULL)
	    {
			if(grid2vtk->b_share_vertex)
			{
				writeModflowGrid2VTK_shareV(modflow_grid, grid2vtk->vtkfile, grid2vtk->one_based_node_numbering);
			}
			else
			{
				writeModflowGrid2VTK(modflow_grid, grid2vtk->vtkfile, grid2vtk->one_based_node_numbering);
			}
	    }
	    else{
	    	cerr<<"! Error: Unknown Grid type"<<endl;
	    }
	}
	else if(grid2usg!=NULL)
	{
        DefinitionExporter def_exporter;
        def_exporter.save(*grid2usg);
        def_exporter.save_node_coordinate(grid2usg->usg_data_prefix,grid2usg->grid);
	}
	else if(gridintersection!=NULL)
	{
        GridIntersection gi(0,gridintersection->layer);
        gi.set_save_intersection_info(true);
        gi.build(*gridintersection);
        //save results
        gi.save(gridintersection->output_file, gridintersection->attributes);
		gi.saveShp(gridintersection->output_shapefile, gridintersection->attributes);
	}

#if GL_GUI

	if( (qtree_grid!=NULL || modflow_grid!=NULL || lgr_grid!=NULL) && gl_display)
	{
		//show rendering
        cusgGridGUI gui(dfn_file);

        //populate refinement features
        if(qtree_grid!=NULL && qtree_builder_data!=NULL)
        {
            fillinFeatures(gui,qtree_builder_data);
        }

        //Retrieve the basic modflow grid
		ModflowGrid * mfg = NULL;
		if(modflow_grid!=NULL) mfg = modflow_grid;
		else if(qtree_grid!=NULL) mfg = qtree_grid->getModflowGrid();
		else if(lgr_grid!=NULL) mfg = lgr_grid->getModflowGrid();
        
        double ratio=mfg->height/mfg->width;
        double winwidth=0, winheigh=0;
        int long_side=700;
        if(ratio<1)
        {
            winwidth=long_side;
            winheigh=long_side; //winwidth*ratio;
            gui.setScale(mfg->width);
        }
        else{
            winheigh=long_side;
            winwidth=long_side; //winheigh/ratio;
            gui.setScale(mfg->height);
        }

        Point2d offset(-mfg->xoffset-mfg->width/2,-mfg->yoffset-mfg->height/2);
        gui.setOffset(offset);
        gui.setWindowSize(winwidth,winheigh);

        drawGridBase * renderer=NULL;
        if(qtree_grid!=NULL)
        {
            renderer=new drawQuadtree3D(qtree_grid);

        }
		else if(lgr_grid!=NULL)
		{
			renderer = new drawLGRGrid(lgr_grid);
		}
        else if(modflow_grid!=NULL)
        {
            renderer=new drawModFlowGrid(modflow_grid);
        }

        assert(renderer);
        gui.setGridRenderer(renderer);

        gui.setNumberofNodes(grid->nodes());
        gui.setNumberofLayers(mfg->nlay);
        CSRData csr=grid->get_usg_csr_data();
        gui.setNumberofConnections(csr.nja);
        csr.destory();
        gui.startGUI();
	}

#endif

	return 0;
}

#endif

#if GL_GUI
void fillinFeatures(cusgGridGUI& gui, quadtree_builder_raw_data * qtree_data)
{
    //check the refinement rules
    map<string,ShapeReader> shp_readers;
    for(int i=0;i<qtree_data->modflow_data->nlay;i++)
    {
        vector<c_plyline> all_arcs;
        vector<c_polygon> all_plys;
        vector<Point2d> all_pts;

        list<refinement_raw_data*> rF=qtree_data->refine_by_features[i];
        for(list<refinement_raw_data*>::iterator f=rF.begin();f!=rF.end();f++)
        {
            refinement_raw_data* rf=*f;
            ShapeReader& shp_reader=shp_readers[rf->shapefile];
            shp_reader.read(rf->shapefile);
            //only support
            if(rf->featuretype=="line")
            {
                vector<GIS_plyline>& arcs=shp_reader.getArcList();
                for(vector<GIS_plyline>::iterator i=arcs.begin();i!=arcs.end();i++)
                    all_arcs.push_back( (c_plyline)*i );
            }//end if
            else if(rf->featuretype=="polygon")
            {
                vector<GIS_polygon>& polygons=shp_reader.getPolyList();
                //all_plys.insert(all_plys.end(),polygons.begin(),polygons.end());
                for(vector<GIS_polygon>::iterator i=polygons.begin();i!=polygons.end();i++)
                    all_plys.push_back( (c_polygon)*i );
            }
            else if(rf->featuretype=="point")
            {
                vector<GIS_Point2d>& pts=shp_reader.getPointList();
                //all_pts.insert(all_pts.end(),pts.begin(),pts.end());
                for(vector<GIS_Point2d>::iterator i=pts.begin();i!=pts.end();i++)
                    all_pts.push_back( (Point2d)*i );
            }

        }//end for f

        cusgRefinementFeatures features;
        features.feature_arcs.swap(all_arcs);
        features.feature_polygons.swap(all_plys);
        features.feature_points.swap(all_pts);

        gui.addFeature(features,i);
    }


}

#endif




