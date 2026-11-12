
#include "def_exporter.h"
#include "ModflowGrid.h"
#include "QuadTree3D.h"
#include "QuadTreeBuilder.h"
#include "GridVertices.h"
#include "VertexExport.h"
#include "Warnings.h"
#include <sstream>
#include <algorithm>
#include <filesystem>
#include <vector>

namespace cusg
{

inline ostream & operator<<(ostream & out, list<string> & strlist)
{
    for (auto & i : strlist)
    {
        out << " " << i;
    }
    return out;
}

string DefinitionExporter::output_file_name(const string & file_name,
                                            const string & read_folder) const
{
    if (read_folder == m_output_folder) return file_name;

    namespace fs = std::filesystem;
    string path = file_name;
    if (!read_folder.empty()) path = read_folder + "/" + file_name;
    replace(path.begin(), path.end(), '\\', '/');
    if (m_output_folder.empty()) return path;

    fs::path target = fs::absolute(fs::path(path)).lexically_normal();
    fs::path base = fs::absolute(fs::path(m_output_folder)).lexically_normal();
    fs::path relative = target.lexically_relative(base);
    //a file on another drive than the folder has no relative name
    return relative.empty() ? target.generic_string() : relative.generic_string();
}

list<string> DefinitionExporter::output_array(const list<string> & spec,
                                              const string & read_folder) const
{
    list<string> record = spec;
    string keyword = record.empty() ? "" : record.front();
    for (char & c : keyword) c = static_cast<char>(::tolower(static_cast<unsigned char>(c)));
    if (keyword == "open/close") record.back() = output_file_name(record.back(), read_folder);
    return record;
}

void DefinitionExporter::set_output_folder(const string & folder)
{
    m_output_folder = folder;
    replace(m_output_folder.begin(), m_output_folder.end(), '\\', '/');
}

bool DefinitionExporter::save_grid_definition(const string & filename, QuadtreeBuilderBlock & data)
{
    size_t slash = filename.find_last_of("/\\");
    set_output_folder((slash == string::npos) ? "" : filename.substr(0, slash));
    bool saved = save(filename, data);
    set_output_folder("");
    return saved;
}

bool DefinitionExporter::save(ostream & out, ModflowGridBlock & data)
{
    assert(data.type == "modflow_grid");
    out << "BEGIN MODFLOW_GRID " << data.name << "\n";
    out << "\tROTATION_ANGLE = " << RoundTrip(data.rotation_degrees) << "\n";
    out << "\tX_OFFSET = " << RoundTrip(data.x_offset) << "\n";
    out << "\tY_OFFSET = " << RoundTrip(data.y_offset) << "\n";
    out << "\tLENGTH_UNIT = " << data.length_unit << "\n";
    out << "\tNLAY = " << data.nlay << "\n";
    out << "\tNROW = " << data.nrow << "\n";
    out << "\tNCOL = " << data.ncol << "\n";
    //the arrays are read relative to the folder of the definition file of the block
    const string & read_folder = data.folder_path;
    list<string> delr = output_array(data.delr_str, read_folder);
    list<string> delc = output_array(data.delc_str, read_folder);
    list<string> top = output_array(data.top_str, read_folder);
    out << "\tDELR = " << delr << "\n";
    out << "\tDELC = " << delc << "\n";
    out << "\tTOP = " << top << "\n";

    for (auto & i : data.bot_str)
    {
        list<string> bottom = output_array(i.second, read_folder);
        out << "\tBOTTOM LAYER " << i.first + 1 << " = " << bottom << "\n";
    }

    out << "END MODFLOW_GRID\n";

    return true;
}

bool DefinitionExporter::save(ostream & out, QuadTree3D * grid)
{
    ModflowGrid * modflow_grid = grid->get_modflow_grid();
    //filenames, relative to the folder gridgen is run from; the name of the grid is that of
    //the file of its definition, which may include a folder
    const string & name = grid->name;
    size_t slash = name.find_last_of("/\\");
    string block_name = (slash == string::npos) ? name : name.substr(slash + 1);

    string tree_structure_file = name + ".tsf";

    //start to generate the block
    out << "BEGIN QUADTREE " << block_name << "\n";
    out << "\tMODFLOW_GRID = " << modflow_grid->name << "\n";
    out << "\tSTRUCTURE_File  = OPEN/CLOSE " << output_file_name(tree_structure_file, "") << "\n";

    //save quadtree structure...
    {
        ofstream fout(tree_structure_file.c_str());
        if (fout.good() == false) throw GridgenError("cannot write file " + tree_structure_file);
        fout << (*grid);
        fout.close();
        cout << "- Save Quadtree to file: " << tree_structure_file << endl;
    }

    //
    // Note that the following code works because nodes in quadtree are numbered layer by layer
    //

    //output tops
    double * top_ptr = grid->top;
    for (int i = 0; i < modflow_grid->nlay; i++)
    {
        int nsize_per_layer = grid->node_group.layer_cell_counts[i];
        stringstream ss;
        ss << name << ".top" << i + 1 << ".dat";
        out << "\tTOP LAYER " << i + 1 << " = OPEN/CLOSE " << output_file_name(ss.str(), "")
            << "\n";
        save(ss.str(), top_ptr, nsize_per_layer);
        top_ptr += nsize_per_layer;
    }

    //output bottoms
    double * bot_ptr = grid->bot;
    for (int i = 0; i < modflow_grid->nlay; i++)
    {
        int nsize_per_layer = grid->node_group.layer_cell_counts[i];
        stringstream ss;
        ss << name << ".bot" << i + 1 << ".dat";
        out << "\tBOTTOM LAYER " << i + 1 << " = OPEN/CLOSE " << output_file_name(ss.str(), "")
            << "\n";
        save(ss.str(), bot_ptr, nsize_per_layer);
        bot_ptr += nsize_per_layer;
    }

    out << "END QUADTREE\n";

    return true;
}

//save csr data into a set of files with "prefix"
bool DefinitionExporter::save(const string & prefix, CSRData & csr, int * nperlay,
                              double rotation_degrees)
{
    //filenames...
    string ia_data_name = prefix + ".ia.dat";
    string ja_data_name = prefix + ".ja.dat";
    string connection_direction_array = prefix + ".fldr.dat";
    string cell_area_array = prefix + ".area.dat";
    string connection_area_array = prefix + ".fahl.dat";
    string connection_length_array = prefix + ".c1.dat";
    string connection_length_array_reversed = prefix + ".c2.dat";
    string connection_count_array = prefix + ".iac.dat";
    string horizontal_connection_array = prefix + ".ihc.dat";
    string vertical_connection_array = prefix + ".ivc.dat";
    string face_width_area_array = prefix + ".hwva.dat";
    string face_normal_angle_array = prefix + ".angldegx.dat";

    //needed information about number of nodes per layer here

    //save csr data
    this->set_record_per_line(200);
    save(ia_data_name, csr.ia.data(), csr.nodes + 1);
    save(ja_data_name, csr.ja.data(), csr.nja, csr.ia.data());
    save(connection_direction_array, csr.fldr.data(), csr.nja, csr.ia.data());
    save(cell_area_array, csr.area.data(), csr.nodes, nperlay);
    save(connection_area_array, csr.fahl.data(), csr.nja, csr.ia.data());
    save(connection_length_array, csr.cl1.data(), csr.nja, csr.ia.data());
    save(connection_length_array_reversed, csr.cl2.data(), csr.nja, csr.ia.data());
    if (!csr.iac.empty()) save(connection_count_array, csr.iac.data(), csr.nodes);
    save(horizontal_connection_array, csr.ihc().data(), csr.nja, csr.ia.data());
    save(vertical_connection_array, csr.ivc().data(), csr.nja, csr.ia.data());
    save(face_width_area_array, csr.hwva.data(), csr.nja, csr.ia.data());
    save(face_normal_angle_array, csr.angldegx(rotation_degrees).data(), csr.nja, csr.ia.data());

    return true;
}

bool DefinitionExporter::save(ostream & out, QuadtreeBuilderBlock & data)
{
    bool r = save(out, *data.modflow_data);
    if (r == false) return false;
    out << "\n";
    return save(out, data.qtree_grid);
}

bool DefinitionExporter::save(GridToUsgdataBlock & data)
{
    QuadTree3D * quad_grid = dynamic_cast<QuadTree3D *>(data.grid);

    //needed information about number of nodes per layer here
    int * nperlay = nullptr;
    vector<int> mf_nperlay; //nodes per layer of a MODFLOW grid
    double rotation_degrees = 0.0;
    if (quad_grid != nullptr)
    {
        nperlay = quad_grid->get_layer_cell_counts();
        rotation_degrees = quad_grid->get_modflow_grid()->rotation_degrees;
    }
    else //get number of nodes per layer on the fly
    {
        ModflowGrid * modflow_grid = dynamic_cast<ModflowGrid *>(data.grid);
        if (modflow_grid != nullptr)
        {
            mf_nperlay.assign(modflow_grid->nlay, modflow_grid->ncol * modflow_grid->nrow);
            nperlay = mf_nperlay.data();
            rotation_degrees = modflow_grid->rotation_degrees;
        }
    }

    //vertical pass-through applies only to a quadtree grid
    if (data.vertical_pass_through)
    {
        if (quad_grid != nullptr)
            quad_grid->vertical_pass_through = true;
        else
            warn("Record VERTICAL_PASS_THROUGH in " + data.label() +
                 " is ignored; it applies only to a QUADTREE grid.");
    }

    //save the rest
    assert(data.grid);
    CSRData csr = data.grid->get_usg_csr_data();
    save(data.usg_data_prefix, csr, nperlay, rotation_degrees);
    write_vertices_cell2d(data.grid, data.usg_data_prefix);

    //save extra information for irregular map
    if (quad_grid != nullptr)
    {
        // Save ghost nodes to file
        string ghostnode_name = data.usg_data_prefix + ".gnc.dat";
        quad_grid->save_ghost_node(ghostnode_name, data.gnc_vertical);

        // save nodes per layer
        string nodesperlay_data_name = data.usg_data_prefix + ".nodesperlay.dat";
        ModflowGrid * modflow_grid = quad_grid->get_modflow_grid();
        save(nodesperlay_data_name, nperlay, modflow_grid->nlay);
    }

    return true;
}

bool DefinitionExporter::save_node_coordinate(const string & prefix, Grid * grid)
{
    int id_offset = 1;
    string filename = prefix + ".nod";
    ofstream fout(filename.c_str());
    if (fout.good() == false) throw GridgenError("cannot write file " + filename);

    ///////////////////////////////////////////////////////
    QuadTree3D * qtree_grid = dynamic_cast<QuadTree3D *>(grid);
    ModflowGrid * modflow_grid = dynamic_cast<ModflowGrid *>(grid);

    //count the total number of cells needed to be saved
    int total = 0, connection_count = 0;
    for (auto b : grid->node_group)
    {
        if (b->active && b->is_leaf)
        {
            total++;
            if (qtree_grid != nullptr)
            {
                if (qtree_grid->vertical_pass_through)
                {
                    vector<double> con_areas;
                    connection_count +=
                        qtree_grid->get_node_connections_with_pass_through(b, &con_areas).size();
                }
                else
                    connection_count += qtree_grid->get_node_connections(b).size();
            }
            else if (modflow_grid != nullptr)
            {
                connection_count += modflow_grid->get_node_connections(b->id).size();
            }
        }
    }

    //save the NACTIVECELLS and NCONNECTIONS  to the file
    fout << total << "\t" << connection_count << endl;

    vector<pair<long, Box *>> sorted_boxes;
    sorted_boxes.reserve(total);

    for (auto box : grid->node_group)
    {
        if (box->active == false || box->is_leaf == false) continue;
        sorted_boxes.push_back(make_pair(box->number + id_offset, box));
    }

    std::sort(sorted_boxes.begin(), sorted_boxes.end());

    // cell centers, rotated and offset
    GridVertices gv(grid);

    for (auto & entry : sorted_boxes)
    {
        Box * box = entry.second;
        const Point2d & center = gv.center(box);
        fout << entry.first << " " << box->layer + id_offset << " " << RoundTrip(center[0]) << " "
             << RoundTrip(center[1]) << " " << RoundTrip(box->z) << " " << RoundTrip(box->dx) << " "
             << RoundTrip(box->dy) << " " << RoundTrip(box->dz) << "\n";
    }

    fout.close();
    cout << "- Save data to file: " << filename << endl;
    return true;
}

} //end namespace
