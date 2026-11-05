#include "def_parser.h"
#include "Warnings.h"

#include <cctype>
#include <memory>
#include <set>

namespace cusg
{

//records accepted in each block type; a record given for a layer is listed
//as "<name> layer", and both the report and the code spellings are listed
static const map<string, set<string>> & block_records()
{
    static map<string, set<string>> records;
    if (!records.empty()) return records;

    const char * modflow_grid[] = {"rotation_angle", "x_offset", "y_offset",     "length_unit",
                                   "nlay",           "nrow",     "ncol",         "delr",
                                   "delc",           "top",      "bottom layer", nullptr};
    const char * quadtree[] = {"modflow_grid", "structure_file", "top layer", "bottom layer",
                               nullptr};
    const char * quadtree_builder[] = {"modflow_grid",
                                       "grid_definition_file",
                                       "refinement_features layer",
                                       "refine_by_features layer",
                                       "active_domain layer",
                                       "smoothing",
                                       "horizontal_smoothing_level",
                                       "smoothing_level_horizontal",
                                       "vertical_smoothing_level",
                                       "smoothing_level_vertical",
                                       "max_smoothing_level_target",
                                       "top layer",
                                       "bottom layer",
                                       "autoalignment",
                                       "node_numbering",
                                       nullptr};
    const char * refinement_features[] = {"shapefile",
                                          "feature_type",
                                          "refinement_level",
                                          "refine_level_by_attribute",
                                          "refinement_level_by_attribute",
                                          nullptr};
    const char * active_domain[] = {"shapefile", "feature_type", "include_boundary", nullptr};
    const char * grid_to_shapefile[] = {
        "grid",           "shapefile", "feature_type", "node_numbering", "without_inactive",
        "close_polygons", "concide",   nullptr};
    const char * grid_to_usgdata[] = {"grid", "usg_data_prefix", "vertical_pass_through", nullptr};
    const char * grid_to_vtkfile[] = {"grid",           "vtkfile",          "share_vertex",
                                      "node_numbering", "without_inactive", nullptr};
    const char * grid_intersection[] = {"grid",         "layer",       "shapefile",
                                        "feature_type", "output_file", "output_shapefile",
                                        "attributes",   nullptr};

    struct
    {
        const char * type;
        const char ** names;
    } blocks[] = {
        {"modflow_grid", modflow_grid},          {"quadtree", quadtree},
        {"quadtree_builder", quadtree_builder},  {"refinement_features", refinement_features},
        {"active_domain", active_domain},        {"grid_to_shapefile", grid_to_shapefile},
        {"grid_to_usgdata", grid_to_usgdata},    {"grid_to_vtkfile", grid_to_vtkfile},
        {"grid_intersection", grid_intersection}};
    for (auto & block : blocks)
    {
        set<string> & known = records[block.type];
        for (const char ** n = block.names; *n != nullptr; n++) known.insert(*n);
    }
    return records;
}

//records a block type reads but does not use, with the reason given in the warning
static const map<string, map<string, string>> & ignored_records()
{
    static map<string, map<string, string>> records;
    if (!records.empty()) return records;
    const char * no_array = "GRIDGEN does not refine by array; use REFINEMENT_FEATURES";
    records["quadtree_builder"]["refinement_array layer"] = no_array;
    records["quadtree_builder"]["refine_by_array layer"] = no_array;
    records["grid_intersection"]["output_vtkfile"] =
        "GRIDGEN does not write a VTK file of an intersection";
    return records;
}

static string upper(string s)
{
    for (char & i : s) i = (char)toupper((unsigned char)i);
    return s;
}

//a warning for each record in block that its block type does not accept
static list<string> unrecognized_records(const block_raw_data & block)
{
    list<string> warnings;
    const map<string, set<string>> & records = block_records();
    map<string, set<string>>::const_iterator known = records.find(block.type);
    if (known == records.end()) return warnings;

    for (block_raw_data::const_iterator e = block.begin(); e != block.end(); e++)
    {
        //a record given for a layer is matched without its layer number
        string name = e->name;
        size_t pos = name.rfind(" layer ");
        if (pos != string::npos) name = name.substr(0, pos + 6);
        if (known->second.count(name)) continue;
        map<string, map<string, string>>::const_iterator ignored =
            ignored_records().find(block.type);
        if (ignored != ignored_records().end() && ignored->second.count(name))
        {
            warnings.push_back("record " + upper(e->name) + " in " + upper(block.type) + " block " +
                               block.name + " is ignored; " + ignored->second.at(name));
            continue;
        }
        warnings.push_back("record " + upper(e->name) + " is not recognized in " +
                           upper(block.type) + " block " + block.name + " and is ignored");
    }
    return warnings;
}

//template<typename T> block_raw_data* createBlockType(block_raw_data& block) { return new T(block); }

//-----------------------------------------------------------------------------
//
//
// DefinitionParser
//
//
//-----------------------------------------------------------------------------

//constructor
DefinitionParser::DefinitionParser()
{
    m_block_type_map["modflow_grid"] = &createBlockType<modflow_grid_raw_data>;
    m_block_type_map["quadtree"] = &createBlockType<quadtree_grid_raw_data>;
    m_block_type_map["quadtree_builder"] = &createBlockType<quadtree_builder_raw_data>;
    m_block_type_map["refinement_features"] = &createBlockType<refinement_raw_data>;
    m_block_type_map["active_domain"] = &createBlockType<active_domain_raw_data>;
    m_block_type_map["grid_to_shapefile"] = &createBlockType<gridtoshape_raw_data>;
    m_block_type_map["grid_to_usgdata"] = &createBlockType<grid_to_usgdata_raw_data>;
    m_block_type_map["grid_intersection"] = &createBlockType<grid_intersection_raw_data>;
    m_block_type_map["grid_to_vtkfile"] = &createBlockType<gridtovtk_raw_data>;
}

//destructor
DefinitionParser::~DefinitionParser()
{
    for (auto & entry : m_blocks)
    {
        delete entry.second;
    }

    m_blocks.clear();
    m_block_type_map.clear();
}

bool DefinitionParser::parse(string filename)
{
    ifstream fin(filename.c_str());
    if (!fin.good()) throw GridgenError("cannot open definition file " + filename);

    bool r = parse(fin);
    fin.close();
    return r;
}

//parsing the ws file
bool DefinitionParser::parse(istream & in, string folder)
{
    while (!in.eof())
    {
        //parse raw data
        block_raw_data block;
        bool done = readBlock(in, block);
        if (done == false) break;

        block.folderPath = folder; //whether it is a relative path

        if (m_blocks.find(block.name) != m_blocks.end()) //name exists
        {
            warn("Block " + block.name +
                 " is defined more than once; only the last definition is used.");
            delete m_blocks[block.name];
            m_blocks[block.name] = nullptr;
        }

        //load the related path's .dfn file
        if (block.type == "load")
        {
            string dfnPath = folder;
            dfnPath += "/";
            if (folder.empty() || folder == "") dfnPath = "";

            dfnPath += block.name;

            int pathsize = dfnPath.size();
            for (int i = 0; i < pathsize; i++)
            {
                if (dfnPath[i] == '\\') dfnPath[i] = '/';
            }

            ifstream relfin(dfnPath.c_str());
            if (!relfin.good())
                throw GridgenError("cannot open included definition file " + dfnPath);

            int dotIdx = dfnPath.find_last_of('/');
            string relFolder = "";
            if (dotIdx != -1) relFolder = dfnPath.substr(0, dotIdx);

            parse(relfin, relFolder);
            relfin.close();

            continue;
        }

        //create block data with specific type
        if (m_block_type_map.find(block.type) != m_block_type_map.end())
        {
            list<string> warnings = unrecognized_records(block);
            for (auto & warning : warnings) warn(warning);
            m_warnings.splice(m_warnings.end(), warnings);
            //held until parsed, so a block whose records are in error is freed
            unique_ptr<block_raw_data> ptr_block(m_block_type_map[block.type](block));
            if (ptr_block != nullptr)
            {
                bool r = ptr_block->parse();
                if (r)
                    m_blocks[block.name] = ptr_block.release();
                else
                    throw GridgenError("block " + block.name + " cannot be read");
            } //end if
        }
        else
        {

            warn("Block " + block.name + " has the unrecognized type " + block.type +
                 " and is ignored.");

            continue;
        } //end if

    } //end while

    return true;
}

} //end of namespace
