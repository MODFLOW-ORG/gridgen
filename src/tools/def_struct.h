#pragma once
#ifndef CUSG_DEF_STRUCT_H
#define CUSG_DEF_STRUCT_H

#include <iostream>
#include <list>
#include <vector>
#include <cstdlib>
#include <cstring>
#include <map>
#include <cassert>
#include <fstream>
#include "Basic.h"
#include "reader_treestruct.h"
//#include "LGR_Grid.h"
//#include "LGR_Child_Grid.h"

using namespace std;

namespace cusg
{

//-----------------------------------------------------------------------------
class block_raw_data;         //defined later
class active_domain_raw_data; //defined later in this file
class Grid;                   //defined in Grid.h
class ModflowGrid;            //defined in ModflowGrid.h
class QuadTree3D;             //defined in QuadTree3D.h
class LGR_Child_Grid;
class LGR_Grid;

//-----------------------------------------------------------------------------
//This class parses definition files
//and also sets up and initialize the simulation

class BaseParser 
{
protected:
    
    virtual ~BaseParser(){}

    //read a block
    virtual bool readBlock(istream& in, block_raw_data& block);

    //convert a string to a list of tokens
    virtual list<string> tokenize(char * tmp, const char * ignore);
    
    // split a string into name and value by "="
    pair<char*,char*> split(char * tmp);
    
    //check if the contents of name and tmp are the same
    bool equal(const list<string>& name, string tmp);

    //merge a list of strings in to a single string connected by del
    string merge(const list<string>& name, string del);

    //convert s to lower case string
    string tolower(const string& s);

    template<class T> T * copyArray(T* from, int size)
    {
        T * A=new T[size];
        assert(A);
        for(int i=0;i<size;i++) A[i]=from[i];
        return A;
    }

    //
    template<class T> T * createArray(list<string>& spec, int size, string folderPath="")
    {
        if(spec.empty()) return NULL;

        if(tolower(spec.front())=="constant")
        {
            T * A=new T[size];
            assert(A);
            T value=(T)atof(spec.back().c_str());
            for(int i=0;i<size;i++) A[i]=value;
            return A;
        }
        else if(tolower(spec.front())=="open/close")
        {
            T * A=new T[size];
            assert(A);
			string filePath = folderPath;
			filePath+=spec.back().c_str();
			//replace the '\' with '/'
			int pathsize = filePath.size();
			for (int i = 0; i < pathsize; i++)
			{
				if (filePath[i] == '\\')
				{
					filePath[i] = '/';
				}
			}

			char tmp[1024]; 
            ifstream fin(filePath.c_str());
            if(!fin.good()){
                cerr<<"! Error: Cannot open file:"<<filePath<<endl;
                delete [] A;
            }//end good
            else
            {
                if(fin.peek()=='#'){ fin.getline(tmp, 1024); }
                for(int i=0;i<size;i++) fin>>A[i];
            }

            fin.close();
            return A;
        }
        else 
        {
            //throw spec.back(); //don't know what to do in this function
			cerr<<"! Error: Unknown operator: "<<spec.front()<<" (in BaseParser::createArray)"<<endl;
			exit(1);
        }

        return NULL;
    }
};


//-----------------------------------------------------------------------------
struct entry_raw_data
{
	string name;
	list<string> value;
};

//-----------------------------------------------------------------------------
class block_raw_data : public list<entry_raw_data>, public BaseParser
{
public:

    block_raw_data(){ folderPath = ""; m_initizlied=false; }

	block_raw_data(string fpath){folderPath = fpath; m_initizlied = false;}

    block_raw_data(block_raw_data& other) : BaseParser(), list<entry_raw_data>()
    {
        m_initizlied=false;
        assign(other.begin(),other.end());
        name=other.name;
        type=other.type;
		folderPath = other.folderPath;
    }

    virtual ~block_raw_data() {}

	virtual bool parse(){ return true; } //do nothing by default
    virtual bool initialize(map<string, block_raw_data *>& blocks){ return true; } //do nothing by default

    bool isInitialized() const { return m_initizlied; }
	string getFolderPath() 
	{ 
		if(folderPath.empty() || folderPath == "")
			return "";
		else
			return folderPath + "/"; 
	}

    //data
	string name;
	string type;

	//related folder path
	string folderPath;

protected:

	bool m_initizlied;
};

//-----------------------------------------------------------------------------

class grid_raw_data : public block_raw_data
{
public:
    grid_raw_data(){}
    grid_raw_data(block_raw_data& block):block_raw_data(block){}
    virtual Grid * build()=0;
};

//-----------------------------------------------------------------------------
class modflow_grid_raw_data : public grid_raw_data
{
public:

	modflow_grid_raw_data(block_raw_data& block):grid_raw_data(block)
	{ 
		rotation_angle=x_offset=y_offset=0;
		nlay=nrow=ncol=0;
		delr=delc=NULL;
		top=NULL;
		mfgrid=NULL;
		length_unit="undefined";
	}
	
	virtual bool parse();
	virtual bool initialize(map<string, block_raw_data *>& blocks);
	virtual Grid * build();

	REAL rotation_angle;
	REAL x_offset;
	REAL y_offset;
    string length_unit;
    int nlay;
    int nrow;
    int ncol;

    //raw string data
    list<string> delr_str; //size of ncol
    list<string> delc_str; //size of nrow
    list<string> top_str;  //size of ncol*nrow
    map<int, list<string> > bot_str;

    //numerical data
    ModflowGrid * mfgrid; //the grid built by this class

    REAL * delr; //size of ncol
    REAL * delc; //size of nrow
    REAL * top;  //size of ncol*nrow
    vector<REAL *> bot;
    

    vector<REAL **> layers;
    vector<string> layer_filenames;
};

//-----------------------------------------------------------------------------
class quadtree_grid_raw_data : public grid_raw_data
{
public:

    quadtree_grid_raw_data(block_raw_data& block):grid_raw_data(block)
    {
//        ncon=nodes=0;
//        nodesperlay=ia=ja=connection_direction_array=NULL;
        modflow_data=NULL;
        qtree_grid=NULL;
    }

    virtual bool parse();
    virtual bool initialize(map<string, block_raw_data *>& blocks);
    virtual Grid * build();

    string modflow_grid;
//    int nodes;
//    long ncon;

    //raw string data
    string tree_structure_file;
//    list<string> nodesperlay_str;
//    list<string> ia_str;
//    list<string> ja_str;
//    list<string> connection_direction_array_str;
    map< int, list<string> > top_str;
    map< int, list<string> > bot_str;

    //numerical data
    modflow_grid_raw_data * modflow_data;
    QuadTree3D * qtree_grid; //the quadtree built by this raw data

//    int * nodesperlay; //size of nlay
//    int * ia; //size of nodes+1
//    int * ja; //size of ncon
//    int * connection_direction_array; //size of ncon
    //vector< REAL* > top;
    //vector< REAL* > bot;

	//add by Guilin
	vector<TreeStructNode> allNodes;
	//use the allNodes to fill the QuadTree qtree_grid
	void fillQuadTree3D();
};

//-----------------------------------------------------------------------------
class refinement_raw_data; //defined below
class quadtree_builder_raw_data : public grid_raw_data
{
public:

    quadtree_builder_raw_data(block_raw_data& block):grid_raw_data(block)
    {
        smoothing_level_horizontal=smoothing_level_vertical=1;
        modflow_data=NULL;
        qtree_grid=NULL;
        smoothing=false;
        one_based_node_numbering=true;
		auto_alignment=true;
    }

    virtual bool parse();
    virtual bool initialize(map<string, block_raw_data *>& blocks);
    virtual Grid * build();

    int smoothing_level_horizontal;
	int smoothing_level_vertical;

    string modflow_grid; //name of the modflow_grid
    string grid_definition_file;
    bool smoothing;
    bool one_based_node_numbering;
	bool auto_alignment; //for surface between two layers

    //active domains
    map< int, list<string> > active_domain_str;

    //string node_numbering;
    map< int, list<string> > refine_by_features_str; //one string for each layer
    map< int, list<string> > refine_by_array_str; //one array of size (nrow,ncol) for each layer
    map< int, list<string> > top_str;
    map< int, list<string> > bot_str;

    //numerical data
    modflow_grid_raw_data * modflow_data;

    //
    vector< list<active_domain_raw_data *> >  active_domain_data;

    //the quadtree created by this class
    QuadTree3D * qtree_grid;

    vector< list<refinement_raw_data*> > refine_by_features;
    vector<int*> refine_by_array; //one array of size (nrow,ncol) for each layer
    vector< pair<string,string> > top; //first string is the operation and the second string is the source name
	                                   //the operation can be "REPLICATE", "INTERPOLATE" or "ASCIIGRID"

	vector< pair<string,string> > bot; //first string is the operation and the second string is the source name
	
	//bool areaWeightedInterpolate;
	map<int, bool> top_area_weighted; //the mode of interpolation: "AREA_WEIGHTED"
	map<int, bool> bot_area_weighted; //the mode of interpolation: "AREA_WEIGHTED"

};


//-----------------------------------------------------------------------------
class refinement_raw_data : public block_raw_data
{
public:

    refinement_raw_data(block_raw_data& block):block_raw_data(block)
    {
        refinement_level=1;
    }

    virtual bool parse();

    string shapefile;
    string featuretype;
    string refinement_level_by_attribute;
    int refinement_level;
};

//-----------------------------------------------------------------------------

class surface_interpolation_raw_data : public block_raw_data
{
public:

    surface_interpolation_raw_data(block_raw_data& block):block_raw_data(block)
    {
    }

    virtual bool parse();

    string filename;
    string surface_type;
    string modflow_grid;
    string interpolation_type;
};

//-----------------------------------------------------------------------------

class grid_intersection_raw_data : public block_raw_data
{
public:

    grid_intersection_raw_data(block_raw_data& block):block_raw_data(block)
    {
        grid=NULL;
        layer=0;
    }

    virtual bool parse();
    virtual bool initialize(map<string, block_raw_data *>& blocks);

    string grid_name;  //name of the grid to convert usg_data from
    string shapefile;
    string feature_type;
    string output_file;
	string output_shapefile;
	vector<string> attributes;

	//add on Aug. 22nd
	string output_vtkfile;

    //actual data...
    Grid * grid;
    int layer;
};

//-----------------------------------------------------------------------------

class grid_to_usgdata_raw_data : public block_raw_data
{
public:

    grid_to_usgdata_raw_data(block_raw_data& block):block_raw_data(block)
    {
        grid=NULL;
		vertical_pass_through = false;
    }

    virtual bool parse();
    virtual bool initialize(map<string, block_raw_data *>& blocks);

    string grid_name;  //name of the grid to convert usg_data from
    string usg_data_prefix; //name of the usg_data file
	bool vertical_pass_through;

    //actual data...
    Grid * grid;
};

//----------------------------------------------------------------------------
class gridtoshape_raw_data : public block_raw_data
{
public:

    gridtoshape_raw_data(block_raw_data& block):block_raw_data(block)
	{
        grid=NULL;
	    feature_type="polygon";
		one_based_node_numbering = true;
		without_inactive = true;
		concide = true;
	}

	virtual bool parse();
	virtual bool initialize(map<string, block_raw_data *>& blocks);

	string grid_name;
	//int feature_type;
	string feature_type;

	string shapefile;
	bool one_based_node_numbering;
	bool without_inactive;

    //actual data...
    Grid * grid;

	//first and last vertex concide
	bool concide;
};

//----------------------------------------------------------------------------
class active_domain_raw_data : public block_raw_data
{
public:

    active_domain_raw_data(block_raw_data& block):block_raw_data(block)
    {
        include_boundary=true;
    }

    virtual bool parse();
    bool initialize(map<string, block_raw_data *>& blocks);

    string feature_type;
    string shapefile;

    //actual data...
    bool include_boundary;
};


//----------------------------------------------------------------------------
//added on Aug. 21st
class gridtovtk_raw_data : public block_raw_data
{
public:
	gridtovtk_raw_data(block_raw_data& block):block_raw_data(block)
	{
        grid=NULL;
		b_share_vertex = true;
		one_based_node_numbering = true;
		without_inactive = true;
	}

	virtual bool parse();
	virtual bool initialize(map<string, block_raw_data *>& blocks);

	string grid_name;

	string vtkfile;
	bool one_based_node_numbering;
	bool without_inactive;

	bool b_share_vertex;
    //actual data...
    Grid * grid;
};



}//end namespace

#endif //CUSG_DEF_STRUCT_H
