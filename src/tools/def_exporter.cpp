
#include "def_exporter.h"
#include "ModflowGrid.h"
#include "QuadTree3D.h"
#include "QuadTreeBuilder.h"
#include "GridVertices.h"
#include "VertexExport.h"
#include "Warnings.h"
#include <sstream>
#include <algorithm>
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
    out << "\tDELR = " << data.delr_str << "\n";
    out << "\tDELC = " << data.delc_str << "\n";
    out << "\tTOP = " << data.top_str << "\n";

    for (auto & i : data.bot_str)
    {
        out << "\tBOTTOM LAYER " << i.first + 1 << " = " << i.second << "\n";
    }

    out << "END MODFLOW_GRID\n";

    return true;
}

bool DefinitionExporter::save(ostream & out, QuadTree3D * grid)
{
    ModflowGrid * modflow_grid = grid->get_modflow_grid();
    //filenames...
    const string & name = grid->name;

    string ghostnode_name = name + ".gnc.dat"; //ghost nodes are saved in get_usg_csr_data()

    string tree_structure_file = name + ".tsf";
    string nodesperlay_data_name = name + ".nodesperlay.dat";

    //start to generate the block
    out << "BEGIN QUADTREE " << name << "\n";
    out << "\tMODFLOW_GRID = " << modflow_grid->name << "\n";
    out << "\tSTRUCTURE_File  = OPEN/CLOSE " << tree_structure_file << "\n";

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
        out << "\tTOP LAYER " << i + 1 << " = OPEN/CLOSE " << ss.str() << "\n";
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
        out << "\tBOTTOM LAYER " << i + 1 << " = OPEN/CLOSE " << ss.str() << "\n";
        save(ss.str(), bot_ptr, nsize_per_layer);
        bot_ptr += nsize_per_layer;
    }

    out << "END QUADTREE\n";

    return true;
}

//save csr data into a set of files with "prefix"
bool DefinitionExporter::save(const string & prefix, CSRData & csr, int * nperlay)
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
    if (quad_grid != nullptr)
    {
        nperlay = quad_grid->get_layer_cell_counts();
    }
    else //get number of nodes per layer on the fly
    {
        ModflowGrid * modflow_grid = dynamic_cast<ModflowGrid *>(data.grid);
        if (modflow_grid != nullptr)
        {
            mf_nperlay.assign(modflow_grid->nlay, modflow_grid->ncol * modflow_grid->nrow);
            nperlay = mf_nperlay.data();
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
    save(data.usg_data_prefix, csr, nperlay);
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
