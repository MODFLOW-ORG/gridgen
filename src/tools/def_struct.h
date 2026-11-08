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
#include "GridgenError.h"
#include "reader_treestruct.h"
//#include "LGR_Grid.h"
//#include "LGR_Child_Grid.h"

using namespace std;

namespace cusg
{

//-----------------------------------------------------------------------------
class DefinitionBlock;   //defined later
class ActiveDomainBlock; //defined later in this file
class ModflowGridBlock;  //defined later in this file
class Grid;              //defined in Grid.h
class ModflowGrid;       //defined in ModflowGrid.h
class QuadTree3D;        //defined in QuadTree3D.h
class LGR_Child_Grid;
class LGR_Grid;

//-----------------------------------------------------------------------------
//This class parses definition files
//and also sets up and initialize the simulation

class BaseParser
{
protected:
    virtual ~BaseParser() {}

    //read a block
    virtual bool readBlock(istream & in, DefinitionBlock & block);

    //convert a string to a list of tokens
    virtual list<string> tokenize(char * line, const char * ignore);

    // split a string into name and value by "="
    pair<char *, char *> split(char * line);

    //merge a list of strings in to a single string connected by del
    string merge(const list<string> & name, string del);

    //convert s to lower case string
    string tolower(const string & s);

    //convert s to upper case string
    string toupper(const string & s);

    template <class T> T * copyArray(T * from, int size)
    {
        T * A = new T[size];
        assert(A);
        for (int i = 0; i < size; i++) A[i] = from[i];
        return A;
    }

    //
    template <class T> T * createArray(list<string> & spec, int size, string folderPath = "")
    {
        if (spec.empty()) return nullptr;

        if (tolower(spec.front()) == "constant")
        {
            T * A = new T[size];
            assert(A);
            T value = (T)atof(spec.back().c_str());
            for (int i = 0; i < size; i++) A[i] = value;
            return A;
        }
        else if (tolower(spec.front()) == "open/close")
        {
            string filePath = folderPath;
            filePath += spec.back().c_str();
            //replace the '\' with '/'
            int pathsize = filePath.size();
            for (int i = 0; i < pathsize; i++)
            {
                if (filePath[i] == '\\')
                {
                    filePath[i] = '/';
                }
            }

            char line[1024];
            ifstream fin(filePath.c_str());
            if (!fin.good()) throw GridgenError("cannot open array file " + filePath);

            T * A = new T[size];
            if (fin.peek() == '#')
            {
                fin.getline(line, 1024);
            }
            for (int i = 0; i < size; i++) fin >> A[i];
            if (fin.fail())
            {
                delete[] A;
                throw GridgenError("file " + filePath + " has fewer than " + to_string(size) +
                                   " numbers, or a value that is not a number");
            }

            fin.close();
            return A;
        }
        else
        {
            throw GridgenError("array input " + spec.front() + " is not CONSTANT or OPEN/CLOSE");
        }
    }
};

//-----------------------------------------------------------------------------
struct BlockRecord
{
    string name;
    list<string> value;
};

//-----------------------------------------------------------------------------
class DefinitionBlock : public list<BlockRecord>, public BaseParser
{
public:
    DefinitionBlock()
    {
        folderPath = "";
        m_initialized = false;
    }

    DefinitionBlock(string fpath)
    {
        folderPath = fpath;
        m_initialized = false;
    }

    DefinitionBlock(DefinitionBlock & other) : BaseParser(), list<BlockRecord>()
    {
        m_initialized = false;
        assign(other.begin(), other.end());
        name = other.name;
        type = other.type;
        folderPath = other.folderPath;
    }

    ~DefinitionBlock() override {}

    virtual bool parse() { return true; } //do nothing by default
    virtual bool initialize(map<string, DefinitionBlock *> & blocks)
    {
        return true;
    } //do nothing by default

    bool isInitialized() const { return m_initialized; }
    string getFolderPath()
    {
        if (folderPath.empty() || folderPath == "")
            return "";
        else
            return folderPath + "/";
    }

    //data
    string name;
    string type;

    //related folder path
    string folderPath;

    //whether the block has a record of this name
    bool hasRecord(const string & rname) const
    {
        for (const auto & e : *this)
            if (e.name == rname) return true;
        return false;
    }

    //warnings reported when the block was initialized
    const list<string> & getWarnings() const { return m_warnings; }

    //the block as it is written in the definition file, e.g. "MODFLOW_GRID block mfgrid"
    string label() { return toupper(type) + " block " + name; }

    //whether running the block writes output, as an action block does
    virtual bool isAction() const { return false; }

    //warn that running the block writes nothing, unless it is an action
    void warnIfNotAction();

protected:
    //report that a record given in this block is ignored, and why
    void warnIgnored(const string & rname, const string & reason)
    {
        string w = "record " + toupper(rname) + " in " + label() + " is ignored; " + reason;
        cerr << "! Warning: " << w << endl;
        m_warnings.push_back(w);
    }

    //the error for a record the block requires but does not have
    GridgenError missingRecord(const string & rname)
    {
        return GridgenError(label() + " has no " + toupper(rname) + " record");
    }

    //the error for a record that names a block that is not defined
    GridgenError undefinedBlock(const string & rname, const string & block_name)
    {
        return GridgenError("record " + toupper(rname) + " in " + label() + " names block " +
                            block_name + ", which is not defined");
    }

    //the block named by record rname of this block
    DefinitionBlock * findBlock(map<string, DefinitionBlock *> & blocks, const string & rname,
                                const string & block_name);

    //the error for a record that names a block of the wrong type
    GridgenError wrongBlockType(const string & rname, DefinitionBlock * block,
                                const string & expected)
    {
        return GridgenError("record " + toupper(rname) + " in " + label() + " names " +
                            block->label() + ", which is not " + expected);
    }

    //the MODFLOW_GRID block named by record MODFLOW_GRID of this block, initialized
    ModflowGridBlock * findModflowGrid(map<string, DefinitionBlock *> & blocks,
                                       const string & grid_name);

    //the grid named by record rname of this block, initialized and built
    Grid * buildReferencedGrid(map<string, DefinitionBlock *> & blocks, const string & rname,
                               const string & grid_name);

    //report NODE_NUMBERING and WITHOUT_INACTIVE as ignored for a quadtree grid
    void warnQuadtreeExportRecords(Grid * grid);

    bool m_initialized;
    list<string> m_warnings;
};

//-----------------------------------------------------------------------------

class GridBlock : public DefinitionBlock
{
public:
    GridBlock() {}
    GridBlock(DefinitionBlock & block) : DefinitionBlock(block) {}
    virtual Grid * build() = 0;
};

//-----------------------------------------------------------------------------
class ModflowGridBlock : public GridBlock
{
public:
    ModflowGridBlock(DefinitionBlock & block) : GridBlock(block)
    {
        rotation_angle = rotation_degrees = x_offset = y_offset = 0;
        nlay = nrow = ncol = 0;
        delr = delc = nullptr;
        top = nullptr;
        mfgrid = nullptr;
        length_unit = "undefined";
    }

    ~ModflowGridBlock() override;

    bool parse() override;
    bool initialize(map<string, DefinitionBlock *> & blocks) override;
    Grid * build() override;

    REAL rotation_angle;   //radians
    REAL rotation_degrees; //as given, so it is written back unchanged
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
    map<int, list<string>> bot_str;

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
class QuadtreeBlock : public GridBlock
{
public:
    QuadtreeBlock(DefinitionBlock & block) : GridBlock(block)
    {
        modflow_data = nullptr;
        qtree_grid = nullptr;
    }

    ~QuadtreeBlock() override;

    bool parse() override;
    bool initialize(map<string, DefinitionBlock *> & blocks) override;
    Grid * build() override;

    string modflow_grid;

    //raw string data
    string tree_structure_file;
    map<int, list<string>> top_str;
    map<int, list<string>> bot_str;

    //numerical data
    ModflowGridBlock * modflow_data;
    QuadTree3D * qtree_grid; //the quadtree built by this raw data

    //add by Guilin
    vector<TreeStructNode> allNodes;
    //use the allNodes to fill the QuadTree qtree_grid
    void build_from_structure();
};

//-----------------------------------------------------------------------------
class RefinementFeaturesBlock; //defined below
class QuadtreeBuilderBlock : public GridBlock
{
public:
    bool isAction() const override { return true; }

    QuadtreeBuilderBlock(DefinitionBlock & block) : GridBlock(block)
    {
        smoothing_level_horizontal = smoothing_level_vertical = 1;
        modflow_data = nullptr;
        qtree_grid = nullptr;
        smoothing = false;
        one_based_node_numbering = true;
        auto_alignment = true;
    }

    ~QuadtreeBuilderBlock() override;

    bool parse() override;
    bool initialize(map<string, DefinitionBlock *> & blocks) override;
    Grid * build() override;

    int smoothing_level_horizontal;
    int smoothing_level_vertical;

    string modflow_grid; //name of the modflow_grid
    string grid_definition_file;
    bool smoothing;
    bool one_based_node_numbering;
    bool auto_alignment; //for surface between two layers

    //active domains
    map<int, list<string>> active_domain_str;

    map<int, list<string>> refine_by_features_str; //one string for each layer
    map<int, list<string>> top_str;
    map<int, list<string>> bot_str;

    //numerical data
    ModflowGridBlock * modflow_data;

    //
    vector<list<ActiveDomainBlock *>> active_domain_data;

    //the quadtree created by this class
    QuadTree3D * qtree_grid;

    vector<list<RefinementFeaturesBlock *>> refine_by_features;
    vector<pair<string, string>>
        top; //first string is the operation and the second string is the source name
             //the operation can be "REPLICATE", "INTERPOLATE" or "ASCIIGRID"

    vector<pair<string, string>>
        bot; //first string is the operation and the second string is the source name

    map<int, bool> top_area_weighted; //the mode of interpolation: "AREA_WEIGHTED"
    map<int, bool> bot_area_weighted; //the mode of interpolation: "AREA_WEIGHTED"
};

//-----------------------------------------------------------------------------
class RefinementFeaturesBlock : public DefinitionBlock
{
public:
    RefinementFeaturesBlock(DefinitionBlock & block) : DefinitionBlock(block)
    {
        refinement_level = 1;
    }

    bool parse() override;

    string shapefile;
    string featuretype;
    string refinement_level_by_attribute;
    int refinement_level;
};

//-----------------------------------------------------------------------------

class SurfaceInterpolationBlock : public DefinitionBlock
{
public:
    SurfaceInterpolationBlock(DefinitionBlock & block) : DefinitionBlock(block) {}

    string filename;
    string surface_type;
    string modflow_grid;
    string interpolation_type;
};

//-----------------------------------------------------------------------------

class GridIntersectionBlock : public DefinitionBlock
{
public:
    bool isAction() const override { return true; }

    GridIntersectionBlock(DefinitionBlock & block) : DefinitionBlock(block)
    {
        grid = nullptr;
        layer = 0;
    }

    bool parse() override;
    bool initialize(map<string, DefinitionBlock *> & blocks) override;

    string grid_name; //name of the grid to convert usg_data from
    string shapefile;
    string feature_type;
    string output_file;
    string output_shapefile;
    vector<string> attributes;

    //actual data...
    Grid * grid;
    int layer;
};

//-----------------------------------------------------------------------------

class GridToUsgdataBlock : public DefinitionBlock
{
public:
    bool isAction() const override { return true; }

    GridToUsgdataBlock(DefinitionBlock & block) : DefinitionBlock(block)
    {
        grid = nullptr;
        vertical_pass_through = false;
    }

    bool parse() override;
    bool initialize(map<string, DefinitionBlock *> & blocks) override;

    string grid_name;       //name of the grid to convert usg_data from
    string usg_data_prefix; //name of the usg_data file
    bool vertical_pass_through;

    //actual data...
    Grid * grid;
};

//----------------------------------------------------------------------------
class GridToShapefileBlock : public DefinitionBlock
{
public:
    bool isAction() const override { return true; }

    GridToShapefileBlock(DefinitionBlock & block) : DefinitionBlock(block)
    {
        grid = nullptr;
        feature_type = "polygon";
        one_based_node_numbering = true;
        without_inactive = true;
        close_polygons = true;
    }

    bool parse() override;
    bool initialize(map<string, DefinitionBlock *> & blocks) override;

    string grid_name;
    string feature_type;

    string shapefile;
    bool one_based_node_numbering;
    bool without_inactive;

    //actual data...
    Grid * grid;

    //repeat the first vertex of each cell polygon at its end
    bool close_polygons;
};

//----------------------------------------------------------------------------
class ActiveDomainBlock : public DefinitionBlock
{
public:
    ActiveDomainBlock(DefinitionBlock & block) : DefinitionBlock(block) { include_boundary = true; }

    bool parse() override;
    bool initialize(map<string, DefinitionBlock *> & blocks) override;

    string feature_type;
    string shapefile;

    //actual data...
    bool include_boundary;
};

//----------------------------------------------------------------------------
//added on Aug. 21st
class GridToVtkfileBlock : public DefinitionBlock
{
public:
    bool isAction() const override { return true; }

    GridToVtkfileBlock(DefinitionBlock & block) : DefinitionBlock(block)
    {
        grid = nullptr;
        b_share_vertex = true;
        one_based_node_numbering = true;
        without_inactive = true;
    }

    bool parse() override;
    bool initialize(map<string, DefinitionBlock *> & blocks) override;

    string grid_name;

    string vtkfile;
    bool one_based_node_numbering;
    bool without_inactive;

    bool b_share_vertex;
    //actual data...
    Grid * grid;
};

} //end namespace

#endif //CUSG_DEF_STRUCT_H
