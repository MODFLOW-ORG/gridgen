#include <ctype.h>
#include <climits>
#include <set>
#include <sstream>
#include "def_struct.h"
#include "Warnings.h"
#include "ModflowGrid.h"
#include "QuadTree3D.h"
#include "QuadTreeBuilder.h"
#include "reader_treestruct.h"
#include "GridToShp.h"

namespace cusg
{

//indices (zero based) of the layers named in "... layer n" records
static set<int> named_layers(map<string, list<string>> & name_value_map)
{
    set<int> layers;
    for (auto & r : name_value_map)
    {
        size_t pos = r.first.rfind(" layer ");
        if (pos == string::npos) continue;
        long layer = strtol(r.first.c_str() + pos + 7, nullptr, 10);
        if (layer > 0 && layer <= INT_MAX) layers.insert((int)layer - 1);
    }
    return layers;
}

//value of the "<name> layer <i+1>" record, or an empty list if there is none
static list<string> layer_record(map<string, list<string>> & name_value_map, const char * name,
                                 int i)
{
    char layer_name[64];
    sprintf(layer_name, "%s layer %d", name, i + 1);
    map<string, list<string>>::iterator found = name_value_map.find(layer_name);
    if (found == name_value_map.end()) return list<string>();
    return found->second;
}

//value of the record spelled as in the report, or else as in the code
static list<string> record(map<string, list<string>> & name_value_map, const char * report_name,
                           const char * code_name)
{
    map<string, list<string>>::iterator found = name_value_map.find(report_name);
    if (found != name_value_map.end() && !found->second.empty()) return found->second;
    found = name_value_map.find(code_name);
    if (found != name_value_map.end()) return found->second;
    return list<string>();
}

/* ********************************************************************** */
// skip blanks, tabs, line breaks and comment lines,
//  leaving us at the beginning of a token (or EOF)
//  (This code is taken from CORE lib)
inline int skip_comment_line(std::ifstream & in)
{
    int c;

    do
    {
        c = in.get();
        while (c == '#')
        {
            do
            { // ignore the rest of this line
                c = in.get();
            } while (c != '\n');
            c = in.get(); // now, reach the beginning of the next line
        }
    } while (c == ' ' || c == '\t' || c == '\n' || c == '\r'); //ignore white spaces and newlines

    if (c == EOF) warn("The definition file ends unexpectedly.");

    in.putback(c); // this is non-white and non-comment char!
    return c;

} //skip_comment_line

// skips '\' followed by '\n'
//  NOTE: this assumes a very special file format (e.g., our BigInt File format)
//  in which the only legitimate appearance of '\' is when it is followed
//  by '\n' immediately!
//  (This code is taken from CORE lib)
inline int skip_backslash_new_line(std::istream & in)
{
    int c = in.get();

    while (c == '\\')
    {
        c = in.get();

        if (c == '\n' || c == '\r')
            c = in.get();
        else // assuming the very special file format noted above!
            warn("A continuation character \\ in the definition file is not at the end of a "
                 "line and is ignored.");
    } //while
    return c;
} //skip_backslash_new_line

//-----------------------------------------------------------------------------
//
//
// BaseParser
//
//
//-----------------------------------------------------------------------------

bool BaseParser::read_block(istream & in, DefinitionBlock & block)
{
    const int size = 1024;
    const char * ignore = " \t[]()<>,!{}\r\n";
    vector<char> buffer(size);
    char * line = buffer.data();

    //
    while (!in.eof())
    {
        in.getline(line, size);

        if (strlen(line) == 0) continue;

        pair<char *, char *> name_value = split(line);

        list<string> name = tokenize(name_value.first, ignore);

        if (name.empty())
        {
            continue;
        }

        string label = name.front();
        if (label[0] == '#') continue; //comment

        label = tolower(label);

        /**********************************************************************/
        //here we simply regard "load" as a type
        if (label == "load")
        {
            block.type = "load";
            block.name = name.back();
            return true;
        }

        if (label == "begin")
        {
            name.pop_front();
            block.type = tolower(name.front());
            block.name = name.back();
        }
        else if (label == "end")
        {
            name.pop_front();
            string end_type = tolower(name.front());
            if (block.type != end_type)
                throw GridgenError("block " + block.name + " begins as " + toupper(block.type) +
                                   " but ends as " + toupper(end_type));

            break; //we are done for this block
        }
        else
        {
            BlockRecord data;
            data.name = tolower(merge(name, " "));

            data.value = tokenize(name_value.second, ignore);
            block.push_back(data);
        }
    } //end while

    if (block.name.empty()) return false; //read nothing...

    return true;
}

//convert a string to a list of tokens
list<string> BaseParser::tokenize(char * line, const char * ignore)
{
    list<string> tokens;
    char * token = strtok(line, ignore);
    while (token != nullptr)
    {
        tokens.push_back(token);
        token = strtok(nullptr, ignore);
    }
    return tokens;
}

// split a string into name and value by "="
pair<char *, char *> BaseParser::split(char * line)
{
    char * pos = strchr(line, '=');
    if (pos != nullptr)
    {
        *pos = '\0';
        pos++;
    }

    return make_pair(line, pos);
}

string BaseParser::merge(const list<string> & name, string del)
{
    string merged;
    using IT = list<string>::const_iterator;
    IT last = --name.end();
    for (IT i = name.begin(); i != name.end(); i++)
    {
        merged += (*i);
        if (i != last) merged += del;
    }

    return merged;
}

string BaseParser::tolower(const string & s)
{
    string lower;
    for (char i : s)
    {
        char c = ::tolower(i);
        lower.push_back(c);
    }
    return lower;
}

string BaseParser::toupper(const string & s)
{
    string upper;
    for (char i : s) upper.push_back(::toupper(i));
    return upper;
}

//-----------------------------------------------------------------------------
//
//
// ModflowGridBlock
//
//
//-----------------------------------------------------------------------------

bool ModflowGridBlock::parse()
{
    map<string, list<string>> name_value_map;

    for (auto & i : *this)
    {
        string & name = i.name;
        list<string> & value = i.value;

        if (name.empty() || value.empty())
            throw GridgenError("a record in " + label() + " has no name or no value");

        name_value_map[name] = value;
    }

    //check the names
    list<string> & value = name_value_map["rotation_angle"];
    if (!value.empty())
    {
        rotation_degrees = atof(value.front().c_str());
        rotation_angle = rotation_degrees * PI / 180;
    }
    else
    {
        throw missing_record("rotation_angle");
    }

    value = name_value_map["x_offset"];
    if (!value.empty())
        x_offset = atof(value.front().c_str());
    else
    {
        throw missing_record("x_offset");
    }

    value = name_value_map["y_offset"];
    if (!value.empty())
        y_offset = atof(value.front().c_str());
    else
    {
        throw missing_record("y_offset");
    }

    value = name_value_map["length_unit"];
    if (!value.empty())
        length_unit = value.front();
    else
    {
    }

    value = name_value_map["nlay"];
    if (!value.empty())
        nlay = atoi(value.front().c_str());
    else
    {
        throw missing_record("nlay");
    }

    value = name_value_map["nrow"];
    if (!value.empty())
        nrow = atoi(value.front().c_str());
    else
    {
        throw missing_record("nrow");
    }

    value = name_value_map["ncol"];
    if (!value.empty())
        ncol = atoi(value.front().c_str());
    else
    {
        throw missing_record("ncol");
    }

    value = name_value_map["delr"];
    if (!value.empty())
    {
        delr_str = value;
    }
    else
    {
        throw missing_record("delr");
    }

    value = name_value_map["delc"];
    if (!value.empty())
    {
        delc_str = value;
    }
    else
    {
        throw missing_record("delc");
    }

    value = name_value_map["top"];
    if (!value.empty())
    {
        top_str = value;
    }
    else
    {
        throw missing_record("top");
    }

    for (int i = 0; i < nlay; i++)
    {
        char layer_name[36];
        sprintf(layer_name, "bottom layer %d", i + 1);
        value = name_value_map[layer_name];
        if (!value.empty())
        {
            bot_str[i] = value;
            //relative path
        }
        else
        {
            throw missing_record(layer_name);
        }
    }

    return true;
}

ModflowGridBlock::~ModflowGridBlock()
{
    delete mfgrid;
    delete[] delr;
    delete[] delc;
    //bot[0] is top
    for (REAL * b : bot) delete[] b;
}

bool ModflowGridBlock::initialize(map<string, DefinitionBlock *> & blocks)
{
    if (m_initialized) return true;

    //kept only when every array is read, so an error leaves none
    unique_ptr<REAL[]> r(create_array<REAL>(delr_str, ncol, get_folder_path()));
    unique_ptr<REAL[]> c(create_array<REAL>(delc_str, nrow, get_folder_path()));
    vector<unique_ptr<REAL[]>> b;
    b.emplace_back(create_array<REAL>(top_str, ncol * nrow, get_folder_path()));
    for (int i = 0; i < nlay; i++)
    {
        b.emplace_back(create_array<REAL>(bot_str[i], ncol * nrow, get_folder_path()));
    }

    delr = r.release();
    delc = c.release();
    bot.reserve(nlay + 1);
    for (auto & a : b) bot.push_back(a.release());
    top = bot[0];

    m_initialized = true;

    return true;
}

Grid * ModflowGridBlock::build()
{
    if (mfgrid == nullptr)
    {
        mfgrid = new ModflowGrid(*this);
        assert(mfgrid);
        mfgrid->name = this->name;
    }

    return (Grid *)mfgrid;
}

//-----------------------------------------------------------------------------
//
//
// QuadtreeBlock
//
//
//-----------------------------------------------------------------------------

bool QuadtreeBlock::parse()
{
    map<string, list<string>> name_value_map;

    for (auto & i : *this)
    {
        string & name = i.name;
        list<string> & value = i.value;

        if (name.empty() || value.empty())
            throw GridgenError("a record in " + label() + " has no name or no value");

        name_value_map[name] = value;
    }

    //check the names
    list<string> & value = name_value_map["modflow_grid"];
    if (!value.empty())
        modflow_grid = value.front();
    else
    {
        throw missing_record("modflow_grid");
    }

    value = name_value_map["structure_file"];
    if (!value.empty())
        tree_structure_file = value.back(); //value.front();
    else
    {
        throw missing_record("structure_file");
    }

    //read the surfaces of each layer named, independently of the other layers
    set<int> layers = named_layers(name_value_map);
    for (int i : layers)
    {
        list<string> values = layer_record(name_value_map, "top", i);
        if (!values.empty()) top_str[i] = values;
        values = layer_record(name_value_map, "bottom", i);
        if (!values.empty()) bot_str[i] = values;
    }

    return true;
}

bool QuadtreeBlock::initialize(map<string, DefinitionBlock *> & blocks)
{
    if (m_initialized) return true;

    this->modflow_data = find_modflow_grid(blocks, modflow_grid);

    //read/open this->tree_structure_file here
    parse_tree_struct_file(tree_structure_file, all_nodes, get_folder_path());

    m_initialized = true;
    return true;
}

QuadtreeBlock::~QuadtreeBlock()
{
    delete qtree_grid;
}

Grid * QuadtreeBlock::build()
{
    //
    if (this->modflow_data == nullptr)
        throw GridgenError(label() + " is built before the blocks it names are read");

    if (qtree_grid == nullptr) //if it's not already built
    {
        ModflowGrid * base_grid = dynamic_cast<ModflowGrid *>(modflow_data->build());
        assert(base_grid);

        qtree_grid = new QuadTree3D(base_grid);
        try
        {
            //fill the QuadTree3d
            build_from_structure();

            //copy surface data from top/bottom to those in qtree_grid
            int nodesize = qtree_grid->nodes();
            if (qtree_grid->top != nullptr) delete[] qtree_grid->top;
            if (qtree_grid->bot != nullptr) delete[] qtree_grid->bot;
            qtree_grid->top = new double[nodesize];
            qtree_grid->bot = new double[nodesize];
            int index = 0;
            for (int layer = 0; layer < base_grid->nlay; layer++)
            {
                int node_per_layer = qtree_grid->node_group.layer_cell_counts[layer];

                unique_ptr<REAL[]> top_tmp(
                    create_array<REAL>(top_str[layer], node_per_layer, get_folder_path()));
                memcpy(&qtree_grid->top[index], top_tmp.get(), sizeof(double) * node_per_layer);

                unique_ptr<REAL[]> bot_tmp(
                    create_array<REAL>(bot_str[layer], node_per_layer, get_folder_path()));
                memcpy(&qtree_grid->bot[index], bot_tmp.get(), sizeof(double) * node_per_layer);

                index += node_per_layer;
            } //end for layer

            //done
            qtree_grid->name = this->name;
        }
        catch (...)
        {
            //a grid left half built would be returned by the next call
            delete qtree_grid;
            qtree_grid = nullptr;
            throw;
        }
    }

    return (Grid *)qtree_grid;
}

void QuadtreeBlock::build_from_structure()
{
    int count = 0;
    bool zero_based = false;
    for (vector<TreeStructNode>::iterator nite = all_nodes.begin(); nite != all_nodes.end();
         ++nite, count++)
    {

        TreeStructNode & tree_node = *nite;
        //an active cell numbered 0 means the grid was built with zero-based numbering
        if (tree_node.number == 0) zero_based = true;
        int layer_num = tree_node.layer - 1;
        int row = tree_node.row - 1;
        int col = tree_node.col - 1;

        int box_index = qtree_grid->get_modflow_grid()->get_box_index(layer_num, row, col);
        Box * nodeobj = qtree_grid->node_group[box_index];

        if ((strcmp(tree_node.code.c_str(), "") == 0) || tree_node.code.empty())
        {
            nodeobj->number = tree_node.number;
            nodeobj->is_leaf = true;
            if (nodeobj->number == -1) nodeobj->active = false;
        }
        else
        {
            Box * node_box = nodeobj;
            int code_depth = tree_node.code.size();
            int depth_index = 0;
            while (depth_index < code_depth)
            {
                if (node_box->is_leaf)
                {
                    bool was_split = node_box->split(0);
                    assert(was_split);

                    for (int c = 0; c < 4; c++)
                    {
                        int id = c;
                        if (c == 3) id = 2;
                        if (c == 2) id = 3;
                        node_box->p_children[id]->is_leaf = true;
                        node_box->p_children[id]->update_status();
                        qtree_grid->node_group.add(node_box->p_children[id]);
                    }
                    node_box->is_leaf = false;
                }
                char t = tree_node.code[depth_index];
                switch (t)
                {
                    case '1':
                    {
                        node_box = node_box->p_children[0];
                        break;
                    }
                    case '2':
                    {
                        node_box = node_box->p_children[1];
                        break;
                    }
                    case '3':
                    {
                        node_box = node_box->p_children[3];
                        break;
                    }
                    case '4':
                    {
                        node_box = node_box->p_children[2];
                        break;
                    }
                    default:
                    {
                        throw GridgenError("structure file " + tree_structure_file +
                                           " has the quadrant code " + tree_node.code +
                                           ", which is not made of the digits 1 to 4");
                    }
                }
                depth_index++;
            }

            node_box->number = tree_node.number;
            node_box->is_leaf = true;
            if (node_box->number == -1) node_box->active = false;
        }
    }

    qtree_grid->set_one_based_node_numbering(!zero_based);

    //refresh the tree to get the new node id array
    qtree_grid->refresh();

    /********************************************************************/
    //check whether the quadtree is built successfully
}
//-----------------------------------------------------------------------------
//
//
// QuadtreeBuilderBlock
//
//
//-----------------------------------------------------------------------------

bool QuadtreeBuilderBlock::parse()
{
    map<string, list<string>> name_value_map;

    for (auto & i : *this)
    {
        string & name = i.name;
        list<string> & value = i.value;

        if (name.empty() || value.empty())
            throw GridgenError("a record in " + label() + " has no name or no value");

        name_value_map[name] = value;
    }

    //check the names
    list<string> & value = name_value_map["max_smoothing_level_target"];
    if (!value.empty())
        smoothing_level_vertical = smoothing_level_horizontal = atoi(value.front().c_str());

    list<string> level =
        record(name_value_map, "horizontal_smoothing_level", "smoothing_level_horizontal");
    if (!level.empty()) smoothing_level_horizontal = atoi(level.front().c_str());

    level = record(name_value_map, "vertical_smoothing_level", "smoothing_level_vertical");
    if (!level.empty()) smoothing_level_vertical = atoi(level.front().c_str());

    value = name_value_map["modflow_grid"];
    if (!value.empty())
        modflow_grid = value.front();
    else
    {
        throw missing_record("modflow_grid");
    }

    value = name_value_map["grid_definition_file"];
    if (!value.empty())
        grid_definition_file = value.front();
    else
    {
        throw missing_record("grid_definition_file");
    }

    value = name_value_map["smoothing"];
    if (!value.empty()) smoothing = (tolower(value.front()) == "full") ? true : false;

    value = name_value_map["node_numbering"];
    if (!value.empty())
    {
        if (tolower(value.front()) == "zero-based") one_based_node_numbering = false;
    }

    value = name_value_map["autoalignment"];
    if (!value.empty()) auto_alignment = (tolower(value.front()) == "true") ? true : false;

    //read the records of each layer named, independently of the other layers
    set<int> layers = named_layers(name_value_map);
    for (int i : layers)
    {
        list<string> values;

        values = layer_record(name_value_map, "active_domain", i);
        if (!values.empty()) active_domain_str[i] = values;

        //refinement features, also allowed as refinement_features
        values = layer_record(name_value_map, "refine_by_features", i);
        if (values.empty()) values = layer_record(name_value_map, "refinement_features", i);
        if (!values.empty()) refine_by_features_str[i] = values;

        //surface arrays
        values = layer_record(name_value_map, "top", i);
        if (!values.empty()) top_str[i] = values;
        values = layer_record(name_value_map, "bottom", i);
        if (!values.empty()) bot_str[i] = values;
    }

    return true;
}

bool QuadtreeBuilderBlock::initialize(map<string, DefinitionBlock *> & blocks)
{
    if (m_initialized) return true;

    this->modflow_data = find_modflow_grid(blocks, modflow_grid);

    active_domain_data.reserve(modflow_data->nlay);
    for (int i = 0; i < modflow_data->nlay; i++)
    {
        list<ActiveDomainBlock *> domains;
        list<string> & domain_names = active_domain_str[i];
        string rname = "active_domain layer " + to_string(i + 1);
        for (auto & domain_name : domain_names)
        {
            DefinitionBlock * block = find_block(blocks, rname, domain_name);
            auto * domain = dynamic_cast<ActiveDomainBlock *>(block);
            if (domain == nullptr) throw wrong_block_type(rname, block, "an ACTIVE_DOMAIN block");
            domain->initialize(blocks);
            domains.push_back(domain);
        } // end for domain_name
        active_domain_data.push_back(domains);
    }

    // initialize refinement features
    refine_by_features.reserve(modflow_data->nlay);
    top.reserve(modflow_data->nlay);
    bot.reserve(modflow_data->nlay);

    //fill in missing info about the depth values of each layer
    for (int i = 0; i < modflow_data->nlay; i++)
    {
        if (top_str[i].empty())
        {
            top_str[i].push_back("replicate");
            top_str[i].push_back(modflow_data->name);
        }
        if (bot_str[i].empty())
        {
            bot_str[i].push_back("replicate");
            bot_str[i].push_back(modflow_data->name);
        }
    }

    //process each layer now
    for (int i = 0; i < modflow_data->nlay; i++)
    {
        //refine by features
        list<RefinementFeaturesBlock *> refinements;
        list<string> & refinement_names = refine_by_features_str[i];
        string rname = "refinement_features layer " + to_string(i + 1);
        for (auto & refinement_name : refinement_names)
        {
            DefinitionBlock * block = find_block(blocks, rname, refinement_name);
            auto * refinement = dynamic_cast<RefinementFeaturesBlock *>(block);
            if (refinement == nullptr)
                throw wrong_block_type(rname, block, "a REFINEMENT_FEATURES block");
            refinements.push_back(refinement);
        }

        refine_by_features.push_back(refinements);

        //surface data
        top.push_back(make_pair(tolower(top_str[i].front()), top_str[i].back()));
        bot.push_back(make_pair(tolower(bot_str[i].front()), bot_str[i].back()));

        //interpolation mode
        if (top_str[i].size() == 3)
        {
            vector<string> surface(top_str[i].begin(), top_str[i].end());
            if (tolower(surface[1]) == "area_weighted")
                top_area_weighted.insert(make_pair(i, true));
            else if (tolower(surface[1]) == "no_area_weighted")
                top_area_weighted.insert(make_pair(i, false));
        }
        if (bot_str[i].size() == 3)
        {
            vector<string> surface(bot_str[i].begin(), bot_str[i].end());
            if (tolower(surface[1]) == "area_weighted")
                bot_area_weighted.insert(make_pair(i, true));
            else if (tolower(surface[1]) == "no_area_weighted")
                bot_area_weighted.insert(make_pair(i, false));
        }
    }

    m_initialized = true;
    return true;
}

QuadtreeBuilderBlock::~QuadtreeBuilderBlock()
{
    delete qtree_grid;
}

Grid * QuadtreeBuilderBlock::build()
{
    //
    if (this->modflow_data == nullptr)
        throw GridgenError(label() + " is built before the blocks it names are read");

    if (qtree_grid == nullptr) //if it's not already built
    {
        Grid * base_grid = modflow_data->build();

        qtree_grid = new QuadTree3D((ModflowGrid *)base_grid);
        assert(qtree_grid);
        if (this->grid_definition_file.empty()) //no name defined...
            qtree_grid->name = this->name;
        else
        {
            int pos = this->grid_definition_file.find_last_of('.');
            qtree_grid->name = this->grid_definition_file.substr(0, pos);
        }

        QuadtreeBuilder builder;
        try
        {
            builder.build(qtree_grid, *this);
        }
        catch (...)
        {
            //a grid left half built would be returned by the next call
            delete qtree_grid;
            qtree_grid = nullptr;
            throw;
        }
    }

    return (Grid *)qtree_grid;
}

//-----------------------------------------------------------------------------
//
//
// SurfaceInterpolationBlock
//
//
//-----------------------------------------------------------------------------

//-----------------------------------------------------------------------------
//
//
// GridIntersectionBlock
//
//
//-----------------------------------------------------------------------------

bool GridIntersectionBlock::parse()
{
    map<string, list<string>> name_value_map;

    for (auto & i : *this)
    {
        string & name = i.name;
        list<string> & value = i.value;

        if (name.empty() || value.empty())
            throw GridgenError("a record in " + label() + " has no name or no value");

        name_value_map[name] = value;
    }

    //check the names
    list<string> & value = name_value_map["grid"];
    if (!value.empty())
        grid_name = value.front();
    else
    {
        throw missing_record("grid");
    }

    value = name_value_map["shapefile"];
    if (!value.empty())
    {
        shapefile = value.front();
        //replace the '\\'
        string fullpath = get_folder_path();
        fullpath += shapefile;
        shapefile = fullpath;

        int pathsize = shapefile.size();
        for (int i = 0; i < pathsize; i++)
        {
            if (shapefile[i] == '\\') shapefile[i] = '/';
        }
    }
    else
    {
        throw missing_record("shapefile");
    }

    value = name_value_map["feature_type"];
    if (!value.empty())
        feature_type = tolower(value.front());
    else
    {
        throw missing_record("feature_type");
    }
    if (feature_type != "line" && feature_type != "polygon" && feature_type != "point")
        throw GridgenError("FEATURE_TYPE " + value.front() + " (GRID_INTERSECTION " + name +
                           ") is not line, polygon, or point");

    value = name_value_map["output_file"];
    if (!value.empty())
    {
        output_file = value.front();
        string fullpath = get_folder_path();
        fullpath += output_file;
        output_file = fullpath;

        int pathsize = output_file.size();
        for (int i = 0; i < pathsize; i++)
        {
            if (output_file[i] == '\\') output_file[i] = '/';
        }
    }
    else
    {
        throw missing_record("output_file");
    }

    value = name_value_map["output_shapefile"];
    if (!value.empty())
    {
        output_shapefile = value.front();
        string fullpath = get_folder_path();
        fullpath += output_shapefile;

        output_shapefile = fullpath;
        int pathsize = output_shapefile.size();
        for (int i = 0; i < pathsize; i++)
        {
            if (output_shapefile[i] == '\\') output_shapefile[i] = '/';
        }
    }

    value = name_value_map["layer"];
    if (!value.empty())
        layer = atoi(value.front().c_str()) - 1;
    else
    {
        throw missing_record("layer");
    }

    value = name_value_map["attributes"];
    if (!value.empty())
    {
        attributes.insert(attributes.end(), value.begin(), value.end());
    }
    else
    {
        attributes.push_back("all");
    }

    return true;
}

bool GridIntersectionBlock::initialize(map<string, DefinitionBlock *> & blocks)
{
    if (m_initialized) return true;

    this->grid = build_referenced_grid(blocks, "grid", grid_name);

    m_initialized = true;

    return true;
}

//-----------------------------------------------------------------------------
//
//
// ActiveDomainBlock
//
//
//-----------------------------------------------------------------------------

bool ActiveDomainBlock::parse()
{
    map<string, list<string>> name_value_map;

    for (auto & i : *this)
    {
        string & name = i.name;
        list<string> & value = i.value;

        if (name.empty() || value.empty())
            throw GridgenError("a record in " + label() + " has no name or no value");

        name_value_map[name] = value;
    }

    list<string> & value = name_value_map["shapefile"];
    if (!value.empty())
        shapefile = value.front();
    else
    {
        throw missing_record("shapefile");
    }

    value = name_value_map["feature_type"];
    if (!value.empty())
        feature_type = tolower(value.front());
    else
    {
        throw missing_record("feature_type");
    }
    if (feature_type != "line" && feature_type != "polygon" && feature_type != "point")
        throw GridgenError("FEATURE_TYPE " + value.front() + " (ACTIVE_DOMAIN " + name +
                           ") is not line, polygon, or point");

    value = name_value_map["include_boundary"];
    if (!value.empty()) include_boundary = (::tolower(value.front()[0]) == 'f') ? false : true;

    return true;
}

bool ActiveDomainBlock::initialize(map<string, DefinitionBlock *> & blocks)
{
    if (m_initialized) return true;

    m_initialized = true;

    return true;
}

//-----------------------------------------------------------------------------
//
//
// RefinementFeaturesBlock
//
//
//-----------------------------------------------------------------------------

bool RefinementFeaturesBlock::parse()
{
    map<string, list<string>> name_value_map;

    for (auto & i : *this)
    {
        string & name = i.name;
        list<string> & value = i.value;

        if (name.empty() || value.empty())
            throw GridgenError("a record in " + label() + " has no name or no value");

        name_value_map[name] = value;
    }

    //check the names
    list<string> & value = name_value_map["shapefile"];
    if (!value.empty())
        shapefile = value.front();
    else
    {
        throw missing_record("shapefile");
    }

    value = name_value_map["feature_type"];
    if (!value.empty())
        feature_type = value.front();
    else
    {
        throw missing_record("feature_type");
    }

    //a level from a shapefile attribute overrides refinement_level
    value = name_value_map["refinement_level"];
    if (!value.empty()) refinement_level = atoi(value.front().c_str());
    list<string> attribute =
        record(name_value_map, "refine_level_by_attribute", "refinement_level_by_attribute");
    if (!attribute.empty()) refinement_level_by_attribute = attribute.front();
    if (value.empty() && attribute.empty())
    {
        throw GridgenError(label() + " has neither a REFINEMENT_LEVEL nor a "
                                     "REFINEMENT_LEVEL_BY_ATTRIBUTE record");
    }

    return true;
}

//-----------------------------------------------------------------------------
//
//
// GridToUsgdataBlock
//
//
//-----------------------------------------------------------------------------

bool GridToUsgdataBlock::parse()
{
    map<string, list<string>> name_value_map;

    for (auto & i : *this)
    {
        string & name = i.name;
        list<string> & value = i.value;

        if (name.empty() || value.empty())
            throw GridgenError("a record in " + label() + " has no name or no value");

        name_value_map[name] = value;
    }

    //check the names
    list<string> & value = name_value_map["grid"];
    if (!value.empty())
        grid_name = value.front();
    else
    {
        throw missing_record("grid");
    }

    value = name_value_map["usg_data_prefix"];
    if (!value.empty())
        usg_data_prefix = value.front();
    else
    {
        throw missing_record("usg_data_prefix");
    }

    value = name_value_map["vertical_pass_through"];
    if (!value.empty())
    {
        string vs = tolower(value.front());
        if (vs == "true") vertical_pass_through = true;
    }

    value = name_value_map["gnc_vertical"];
    if (!value.empty()) gnc_vertical = (tolower(value.front()) == "true");

    return true;
}

bool GridToUsgdataBlock::initialize(map<string, DefinitionBlock *> & blocks)
{
    if (m_initialized) return true;

    this->grid = build_referenced_grid(blocks, "grid", grid_name);

    //vertical ghost nodes apply only to a quadtree grid
    if (gnc_vertical && dynamic_cast<QuadTree3D *>(this->grid) == nullptr)
    {
        warn_ignored("gnc_vertical", "a MODFLOW_GRID has no cells of different size, so no "
                                     "connection needs a ghost node");
        gnc_vertical = false;
    }

    m_initialized = true;

    return true;
}

//-----------------------------------------------------------------------------
//
//
// GridToShapefileBlock
//
//
//-----------------------------------------------------------------------------

bool GridToShapefileBlock::parse()
{
    map<string, list<string>> name_value_map;

    for (auto & i : *this)
    {
        string & name = i.name;
        list<string> & value = i.value;

        if (name.empty() || value.empty())
            throw GridgenError("a record in " + label() + " has no name or no value");

        name_value_map[name] = value;
    }

    //check the names
    list<string> & value = name_value_map["grid"];
    if (!value.empty())
        grid_name = value.front();
    else
    {
        throw missing_record("grid");
    }

    value = name_value_map["shapefile"];
    if (!value.empty())
    {
        shapefile = value.back(); //value.front();
        string fullpath = get_folder_path();
        fullpath += shapefile;
        shapefile = fullpath;

        int pathsize = shapefile.size();
        for (int i = 0; i < pathsize; i++)
        {
            if (shapefile[i] == '\\') shapefile[i] = '/';
        }
    }
    else
    {
        throw missing_record("shapefile");
    }

    value = name_value_map["node_numbering"];
    if (!value.empty())
    {
        if (tolower(value.front()) == "zero-based") one_based_node_numbering = false;
    }

    value = name_value_map["without_inactive"];
    if (!value.empty())
    {
        if (tolower(value.front()) == "false") without_inactive = false;
    }

    //CONCIDE is the misspelled name of CLOSE_POLYGONS
    value = name_value_map["concide"];
    if (!value.empty())
    {
        warn("CONCIDE (GRID_TO_SHAPEFILE " + name + ") is deprecated; use CLOSE_POLYGONS.");
        close_polygons = tolower(value.front()) != "false";
    }
    value = name_value_map["close_polygons"];
    if (!value.empty()) close_polygons = tolower(value.front()) != "false";

    value = name_value_map["feature_type"];
    if (!value.empty())
        feature_type = tolower(value.front());
    else
    {
        throw missing_record("feature_type");
    }
    if (feature_type != "line" && feature_type != "polygon" && feature_type != "point")
        throw GridgenError("FEATURE_TYPE " + value.front() + " (GRID_TO_SHAPEFILE " + name +
                           ") is not line, polygon, or point");

    return true;
}
bool GridToShapefileBlock::initialize(map<string, DefinitionBlock *> & blocks)
{
    if (m_initialized) return true;

    this->grid = build_referenced_grid(blocks, "grid", grid_name);

    warn_quadtree_export_records(this->grid);

    m_initialized = true;
    return true;
}

void DefinitionBlock::warn_if_not_action()
{
    if (is_action()) return;
    warn(label() + " is not an action, so running it writes nothing; run a QUADTREE_BUILDER, "
                   "GRID_TO_SHAPEFILE, GRID_TO_USGDATA, GRID_TO_VTKFILE, GRID_TO_GSF, or "
                   "GRID_INTERSECTION block to write output");
}

DefinitionBlock * DefinitionBlock::find_block(map<string, DefinitionBlock *> & blocks,
                                              const string & rname, const string & block_name)
{
    auto found = blocks.find(block_name);
    if (found == blocks.end() || found->second == nullptr) throw undefined_block(rname, block_name);
    return found->second;
}

ModflowGridBlock * DefinitionBlock::find_modflow_grid(map<string, DefinitionBlock *> & blocks,
                                                      const string & grid_name)
{
    DefinitionBlock * block = find_block(blocks, "modflow_grid", grid_name);
    auto * modflow_data = dynamic_cast<ModflowGridBlock *>(block);
    if (modflow_data == nullptr)
        throw wrong_block_type("modflow_grid", block, "a MODFLOW_GRID block");
    modflow_data->initialize(blocks);
    return modflow_data;
}

Grid * DefinitionBlock::build_referenced_grid(map<string, DefinitionBlock *> & blocks,
                                              const string & rname, const string & grid_name)
{
    DefinitionBlock * block = find_block(blocks, rname, grid_name);
    GridBlock * grid = dynamic_cast<GridBlock *>(block);
    if (grid == nullptr)
        throw wrong_block_type(rname, block, "a MODFLOW_GRID, QUADTREE, or QUADTREE_BUILDER block");

    grid->initialize(blocks);
    return grid->build();
}

//records that do not change any output of a quadtree grid
void DefinitionBlock::warn_quadtree_export_records(Grid * grid)
{
    if (dynamic_cast<QuadTree3D *>(grid) == nullptr) return;
    if (has_record("node_numbering"))
        warn_ignored("node_numbering",
                     "it has no effect on a QUADTREE grid, which keeps the node numbering "
                     "it was built with; set NODE_NUMBERING in the QUADTREE_BUILDER block");
    if (has_record("without_inactive"))
        warn_ignored("without_inactive", "it has no effect on a QUADTREE grid");
}

//-----------------------------------------------------------------------------
//
//
// GridToVtkfileBlock
//
//
//-----------------------------------------------------------------------------

bool GridToVtkfileBlock::parse()
{
    map<string, list<string>> name_value_map;

    for (auto & i : *this)
    {
        string & name = i.name;
        list<string> & value = i.value;

        if (name.empty() || value.empty())
            throw GridgenError("a record in " + label() + " has no name or no value");

        name_value_map[name] = value;
    }

    //check the names
    list<string> & value = name_value_map["grid"];
    if (!value.empty())
        grid_name = value.front();
    else
    {
        throw missing_record("grid");
    }

    value = name_value_map["vtkfile"];
    if (!value.empty())
    {
        vtkfile = value.back(); //value.front();
        string fullpath = get_folder_path();
        fullpath += vtkfile;
        vtkfile = fullpath;

        int pathsize = vtkfile.size();
        for (int i = 0; i < pathsize; i++)
        {
            if (vtkfile[i] == '\\') vtkfile[i] = '/';
        }
    }
    else
    {
        throw missing_record("vtkfile");
    }

    value = name_value_map["without_inactive"];
    if (!value.empty())
    {
        if (tolower(value.front()) == "false") without_inactive = false;
    }

    value = name_value_map["node_numbering"];
    if (!value.empty())
    {
        if (tolower(value.front()) == "zero-based") one_based_node_numbering = false;
    }

    value = name_value_map["share_vertex"];
    if (!value.empty()) b_share_vertex = (tolower(value.front()) == "true") ? true : false;

    return true;
}
bool GridToVtkfileBlock::initialize(map<string, DefinitionBlock *> & blocks)
{
    if (m_initialized) return true;

    this->grid = build_referenced_grid(blocks, "grid", grid_name);

    warn_quadtree_export_records(this->grid);

    m_initialized = true;
    return true;
}

//-----------------------------------------------------------------------------
//
//
// GridToGsfBlock
//
//
//-----------------------------------------------------------------------------

//a file name relative to the folder of the definition file, with forward slashes
static string output_path(const string & folder, const string & name)
{
    string path = folder + name;
    for (char & c : path)
        if (c == '\\') c = '/';
    return path;
}

bool GridToGsfBlock::parse()
{
    map<string, list<string>> name_value_map;

    for (auto & i : *this)
    {
        string & name = i.name;
        list<string> & value = i.value;

        if (name.empty() || value.empty())
            throw GridgenError("a record in " + label() + " has no name or no value");

        name_value_map[name] = value;
    }

    list<string> & value = name_value_map["grid"];
    if (value.empty()) throw missing_record("grid");
    grid_name = value.front();

    value = name_value_map["gsf_file"];
    if (value.empty()) throw missing_record("gsf_file");
    gsf_file = output_path(get_folder_path(), value.back());

    value = name_value_map["node_data_table_file"];
    if (!value.empty()) node_data_table_file = output_path(get_folder_path(), value.back());

    value = name_value_map["hanging_vertices"];
    if (!value.empty()) hanging_vertices = (tolower(value.front()) == "true");

    return true;
}

bool GridToGsfBlock::initialize(map<string, DefinitionBlock *> & blocks)
{
    if (m_initialized) return true;

    this->grid = build_referenced_grid(blocks, "grid", grid_name);

    if (hanging_vertices && dynamic_cast<QuadTree3D *>(this->grid) == nullptr)
        warn_ignored("hanging_vertices",
                     "the cells of a MODFLOW_GRID are all the same size, so no cell has "
                     "hanging vertices");

    m_initialized = true;
    return true;
}

} //end of namespace
