
#include "def_exporter.h"
#include "ModflowGrid.h"
#include "QuadTree3D.h"
#include "QuadTreeBuilder.h"
#include <sstream>
#include <algorithm>
#include <vector>

#define PRECISION 40
namespace cusg
{

inline void rotateOneBox_exporter(Box* b, double cx, double cy, double rotation, double* padfX, double* padfY, double* padfZ, double* padfM, double* ctrX, double* ctrY)
{
	double ldX, ldY, luX, luY, rdX, rdY, ruX, ruY, ctrx, ctry;
	b->rotate(cx, cy, rotation, ldX, ldY, luX, luY, rdX, rdY, ruX, ruY);
	b->rotateCtr(cx, cy, rotation, ctrx, ctry);

	padfX[0] = luX; padfX[1] = ruX;
	padfX[2] = rdX; padfX[3] = ldX;
	padfY[0] = luY; padfY[1] = ruY;
	padfY[2] = rdY; padfY[3] = ldY;
	ctrX[0] = ctrx; ctrY[0] = ctry;
}

inline ostream& operator<<(ostream& out, list<string>& strlist)
{
    typedef list<string>::iterator IT;
    for(IT i=strlist.begin();i!=strlist.end();i++)
    {
        out<<" "<<*i;
    }
    return out;
}

bool DefinitionExporter::save(ostream& out, grid_raw_data& data)
{
    if( dynamic_cast<modflow_grid_raw_data*>(&data)!=NULL )
        return save(out,(modflow_grid_raw_data&)data);
    else if( dynamic_cast<quadtree_builder_raw_data*>(&data)!=NULL )
        return save(out,(quadtree_builder_raw_data&)data);
    return false;
}

bool DefinitionExporter::save(ostream& out, modflow_grid_raw_data& data)
{
    assert(data.type=="modflow_grid");
	out<<"BEGIN MODFLOW_GRID "<<data.name<<"\n";
	out<<"\tROTATION_ANGLE = "<<setprecision(PRECISION)<<data.rotation_angle * 360 / (2*3.1415926)<<"\n";
	out<<"\tX_OFFSET = "<<setprecision(PRECISION)<<data.x_offset<<"\n";
	out<<"\tY_OFFSET = "<<setprecision(PRECISION)<<data.y_offset<<"\n";
	out<<"\tLENGTH_UNIT = "<<setprecision(PRECISION)<<data.length_unit<<"\n";
	out<<"\tNLAY = "<<data.nlay<<"\n";
	out<<"\tNROW = "<<data.nrow<<"\n";
	out<<"\tNCOL = "<<data.ncol<<"\n";
	out<<"\tDELR = "<<data.delr_str<<"\n";
	out<<"\tDELC = "<<data.delc_str<<"\n";
	out<<"\tTOP = "<<data.top_str<<"\n";

	typedef map<int, list<string> >::iterator IT;
	for(IT i=data.bot_str.begin();i!=data.bot_str.end();i++)
	{
	    out<<"\tBOTTOM LAYER "<<i->first+1<<" = "<<i->second<<"\n";
	}

	out<<"END MODFLOW_GRID\n";

	return true;
}

bool DefinitionExporter::save(ostream& out, ModflowGrid * grid)
{
    const string& name=grid->name;
    string top_data_name=name+".top.dat";
    string delc_data_name=name+".delc.dat";
    string delr_data_name=name+".delr.dat";

    out<<"BEGIN MODFLOW_GRID "<<name<<"\n";
    out<<"\tROTATION_ANGLE = "<<setprecision(PRECISION)<<grid->rotation * 360 / (2*3.1415926)<<"\n";
    out<<"\tX_OFFSET = "<<setprecision(PRECISION)<<grid->xoffset<<"\n";
    out<<"\tY_OFFSET = "<<setprecision(PRECISION)<<grid->yoffset<<"\n";
    //out<<"length_unit = "<<grid-><<"\n"; //not available...
    out<<"\tNLAY = "<<grid->nlay<<"\n";
    out<<"\tNROW = "<<grid->nrow<<"\n";
    out<<"\tNCOL = "<<grid->ncol<<"\n";
    out<<"\tDELR = OPEN/CLOSE "<<delr_data_name<<"\n";
    out<<"\tDELC = OPEN/CLOSE "<<delc_data_name<<"\n";
    out<<"\tTOP = OPEN/CLOSE "<<top_data_name<<"\n";

    save(delc_data_name,grid->delc,grid->nrow);
    save(delc_data_name,grid->delr,grid->ncol);
    save(top_data_name,grid->top,grid->ncol*grid->nrow);

    int csize=grid->nrow*grid->ncol;
    for(int i=0;i<grid->nlay;i++)
    {
        stringstream ss;
        ss<<name<<".bot"<<i+1<<".dat";
        out<<"\tBOTTOM LAYER "<<i+1<<" = OPEN/CLOSE "<<ss.str() <<"\n";
        double * ptr=&grid->bot[(i+1)*csize];
        save(ss.str(),ptr,csize);
    }

    out<<"END MODFLOW_GRID\n";

    return true;
}


bool DefinitionExporter::save(ostream& out, QuadTree3D * grid)
{
	ModflowGrid * mflowgrid=grid->getModflowGrid();
	//filenames...
	const string& name=grid->name;

	string ghostnode_name=name+".gnc.dat";//ghnost nodes are saved in get_usg_csr_data()
    
	//CSRData csr=grid->get_usg_csr_data();

    string tree_structure_file= name+".tsf";
    string nodesperlay_data_name=name+".nodesperlay.dat";

    //start to generate the block
    out<<"BEGIN QUADTREE "<<name<<"\n";
    out<<"\tMODFLOW_GRID = "<<mflowgrid->name<<"\n";
    out<<"\tSTRUCTURE_File  = OPEN/CLOSE "<<tree_structure_file<<"\n";


    //save quadtree structure...
    {
        ofstream fout(tree_structure_file.c_str());
        if(fout.good()==false){
            cerr<<"! Error: Cannot open file:"<<tree_structure_file<<endl;
            return false;
        }
        fout<<(*grid);
        fout.close();
        cout<<"- Save Quadtree to file: "<<tree_structure_file<<endl;
    }

	//
	// Note that the following code works because nodes in quadtree are numbered layer by layer
	//

    //output tops
    double * top_ptr=grid->top;
    for(int i=0;i<mflowgrid->nlay;i++)
    {
        int nsize_per_layer=grid->nodegroup.nodelay_array[i];
        stringstream ss;
        ss<<name<<".top"<<i+1<<".dat";
        out<<"\tTOP LAYER "<<i+1<<" = OPEN/CLOSE "<<ss.str() <<"\n";
        save(ss.str(),top_ptr,nsize_per_layer);
        top_ptr+=nsize_per_layer;
    }

    //output bottoms
    double * bot_ptr=grid->bot;
    for(int i=0;i<mflowgrid->nlay;i++)
    {
        int nsize_per_layer=grid->nodegroup.nodelay_array[i];
        stringstream ss;
        ss<<name<<".bot"<<i+1<<".dat";
        out<<"\tBOTTOM LAYER "<<i+1<<" = OPEN/CLOSE "<<ss.str() <<"\n";
        save(ss.str(),bot_ptr,nsize_per_layer);
        bot_ptr+=nsize_per_layer;
    }

    out<<"END QUADTREE\n";

    //free memory
    //csr.destory();

    return true;
}

//save csr data into a set of files with "prefix"
bool DefinitionExporter::save(const string& prefix, CSRData& csr, int *nperlay)
{
    //filenames...
    string ia_data_name=prefix+".ia.dat";
    string ja_data_name=prefix+".ja.dat";
    string connection_direction_array = prefix+".fldr.dat";
    string cell_area_array = prefix+".area.dat";
    string connection_area_array = prefix+".fahl.dat";
    string connection_length_array = prefix+".c1.dat";
    string connection_length_array_reversed = prefix+".c2.dat";
    string connection_count_array = prefix+".iac.dat";
	
	//needed information about number of nodes per layer here
	

    //save csr data
    this->set_record_per_line(200);
    save(ia_data_name,csr.ia,csr.nodes+1);
    save(ja_data_name,csr.ja,csr.nja,csr.ia);
    save(connection_direction_array,csr.fldr,csr.nja,csr.ia);
    save(cell_area_array,csr.area,csr.nodes,nperlay);
    save(connection_area_array,csr.fahl,csr.nja,csr.ia);
    save(connection_length_array,csr.cl1,csr.nja,csr.ia);
    save(connection_length_array_reversed,csr.cl2,csr.nja,csr.ia);
    if(csr.iac!=NULL) save(connection_count_array,csr.iac,csr.nodes);

    return true;
}

bool DefinitionExporter::save(ostream& out, quadtree_builder_raw_data& data)
{
	bool r=save(out,*data.modflow_data);
	if(r==false) return false;
	out<<"\n";
	return save(out,data.qtree_grid);
}

bool DefinitionExporter::save(grid_to_usgdata_raw_data& data)
{
	QuadTree3D * quadGrid = dynamic_cast<QuadTree3D*>(data.grid);

	//needed information about number of nodes per layer here
	int * nperlay=NULL;
	if(quadGrid != NULL)
	{
		nperlay=quadGrid->get_nodelay();
	}
	else //get number of nodes per lay onfly
	{
		ModflowGrid * mflowgrid=dynamic_cast<ModflowGrid*>(data.grid);
		if(mflowgrid!=NULL)
		{
			nperlay=new int[mflowgrid->nlay];
			assert(nperlay);
			for(int i=0;i<mflowgrid->nlay;i++) nperlay[i]=mflowgrid->ncol*mflowgrid->nrow;
		}
	}

	if(data.vertical_pass_through)
		quadGrid->vertical_pass_through = true;

    //save the rest
    assert(data.grid);
    CSRData csr=data.grid->get_usg_csr_data();
    save(data.usg_data_prefix,csr,nperlay);
    csr.destory(); //free memory

    //save extra infor for irragular map
	if(quadGrid != NULL)
	{
		// Save ghost nodes to file
		string ghostnode_name = data.usg_data_prefix + ".gnc.dat";
		quadGrid->saveGhostNode(ghostnode_name);

		// save nodes per layer
		string nodesperlay_data_name=data.usg_data_prefix+".nodesperlay.dat";
		ModflowGrid * mflowgrid=quadGrid->getModflowGrid();
		save(nodesperlay_data_name,nperlay,mflowgrid->nlay);
	}
	else
	{
		delete [] nperlay;	
	}

    return true;
}

bool DefinitionExporter::save_node_coordinate(const string& prefix, Grid * grid)
{
    int id_offset=1;
    string filename=prefix+".nod";
    ofstream fout(filename.c_str());
    if(fout.good()==false)
    {
        cerr<<"! Error: Failed to open file: "<<filename<<endl;
        return false;
    }

	///////////////////////////////////////////////////////
	QuadTree3D * qtree_grid=dynamic_cast<QuadTree3D*>(grid);
	ModflowGrid * modflow_grid=dynamic_cast<ModflowGrid*>(grid);

	//count the total number of cells needed to be saved
	int total = 0, conNum = 0;
	for(NodeGroup::iterator nit = grid->nodegroup.begin(); nit != grid->nodegroup.end(); nit++)
	{
		Box* b = *nit;
		if(b->active && b->isLeaf)
		{
			total++;
			if(qtree_grid != NULL)
			{
				if(qtree_grid->vertical_pass_through)
				{
					vector<float> conAreas;
					conNum += qtree_grid->get_node_connections_vp(b, &conAreas).size();
				}
				else
					conNum += qtree_grid->get_node_connections(b).size();

				//conNum += qtree_grid->get_node_connections(b).size();
			}
			else if(modflow_grid != NULL)
			{
				conNum += modflow_grid->get_node_connections(b->id).size();
			}
		}
	}

	//save the NACTIVECELLS and NCONNECTIONS  to the file
	fout<<total<<"\t"<<conNum<<endl;

	vector< pair<long, Box *> > sorted_boxes;
	sorted_boxes.reserve(total);

    for(NodeGroup::iterator i=grid->nodegroup.begin();i!=grid->nodegroup.end();i++)
    {
        Box * box=*i;
        if(box->active==false || box->isLeaf==false) continue;
		sorted_boxes.push_back(make_pair(box->number+id_offset,box));
	}

	std::sort(sorted_boxes.begin(),sorted_boxes.end());

	// add rotation here: Mar 15th 2016
	double cx, cy, rotation;
	ModflowGrid * mfgrid=qtree_grid->getModflowGrid();
	mfgrid->getRotatePara(cx, cy, rotation);
	double padfX[5] = {0}, padfY[5] = {0}, padfZ[4] = {0}, padfM[4] = {0}, ctrX[1] = {0}, ctrY[1] = {0};

	//
	for(vector< pair<long, Box *> >::iterator i=sorted_boxes.begin();i!=sorted_boxes.end();i++)
	{
		Box * box=i->second;

		//rotate box
		rotateOneBox_exporter(box, cx, cy, rotation, padfX, padfY, padfZ, padfM, ctrX, ctrY);

		//fout<<i->first<<" "<<box->layer+id_offset
        //    <<" "<<setprecision(PRECISION)<<box->x<<" "<<box->y<<" "<<box->z<<" "
        //    <<box->dx<<" "<<box->dy<<" "<<box->dz<<"\n";
		fout<<i->first<<" "<<box->layer+id_offset
            <<" "<<setprecision(PRECISION)<<ctrX[0]<<" "<<ctrY[0]<<" "<<box->z<<" "
            <<box->dx<<" "<<box->dy<<" "<<box->dz<<"\n";
    }



    fout.close();
    cout<<"- Save data to file: "<<filename<<endl;
    return true;
}

bool DefinitionExporter::save(grid_intersection_raw_data& data)
{
    return true;
}



}//end namespace 
