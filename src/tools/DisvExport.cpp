#include "DisvExport.h"

#include <cmath>
#include <fstream>
#include <iostream>
#include <map>

#include "GhostNode.h"
#include "GridVertices.h"
#include "GridgenError.h"
#include "ModflowGrid.h"
#include "QuadTree3D.h"
#include "RoundTrip.h"

namespace cusg
{

namespace
{

//a leaf cell with the row and column of the base grid cell that contains it
struct LayerCell
{
    Box * box;
    int row;
    int column;
};

//the leaf cells of box in the order they are numbered
void collect_leaves(Box * box, int row, int column, vector<LayerCell> & cells)
{
    if (box->is_leaf)
    {
        cells.push_back({box, row, column});
        return;
    }
    for (int child : {0, 1, 3, 2}) collect_leaves(box->p_children[child], row, column, cells);
}

//true if other is divided into the same cells as first
bool same_cells(const Box * first, const Box * other)
{
    if (first->is_leaf != other->is_leaf) return false;
    if (first->is_leaf) return true;
    for (int child = 0; child < 4; child++)
        if (!same_cells(first->p_children[child], other->p_children[child])) return false;
    return true;
}

ofstream open_output(const string & filename)
{
    ofstream fout(filename.c_str());
    if (fout.good() == false) throw GridgenError("cannot write file " + filename);
    return fout;
}

//write the values of each layer, starting each layer on a new line
template <typename T> void write_layers(ofstream & fout, const vector<T> & values, int ncpl)
{
    const int values_per_line = 10;
    for (size_t i = 0; i < values.size(); i++)
    {
        fout << round_trip_value(values[i]);
        bool line_end = ((i % ncpl) + 1) % values_per_line == 0 || (i + 1) % ncpl == 0;
        fout << (line_end ? "\n" : " ");
    }
}

} // namespace

DisvData disv_data(Grid * grid, bool vertical_pass_through)
{
    QuadTree3D * qtree = dynamic_cast<QuadTree3D *>(grid);
    ModflowGrid * mfgrid =
        (qtree != nullptr) ? qtree->get_modflow_grid() : dynamic_cast<ModflowGrid *>(grid);
    if (mfgrid == nullptr)
        throw GridgenError("the grid to export is neither a MODFLOW grid nor a quadtree grid");
    if (qtree != nullptr) qtree->number_nodes();

    DisvData data;
    data.nlay = mfgrid->nlay;

    //the cells of each layer, in the order of the first layer
    vector<vector<LayerCell>> layers(data.nlay);
    for (int k = 0; k < data.nlay; k++)
    {
        for (int r = 0; r < mfgrid->nrow; r++)
        {
            for (int c = 0; c < mfgrid->ncol; c++)
            {
                Box * root = grid->node_group[mfgrid->get_box_index(k, r, c)];
                Box * first = grid->node_group[mfgrid->get_box_index(0, r, c)];
                if (!same_cells(first, root))
                    throw GridgenError(
                        "layer " + to_string(k + 1) + " of grid " + grid->name +
                        " is refined differently from layer 1, so the grid cannot be written as "
                        "MODFLOW 6 DISV input, which has the same cells in every layer; refine "
                        "every layer with the same features, or write DISU input with a "
                        "GRID_TO_USGDATA block");
                collect_leaves(root, r, c, layers[k]);
            }
        }
    }
    data.ncpl = (int)layers[0].size();

    //vertices numbered in the order the cells use them, with the corners of
    //smaller cells on the edges of each cell, active or not, since a cell
    //inactive in one layer may be active in another
    GridVertices grid_vertices(grid);
    map<int, int> vertex_numbers;
    for (const LayerCell & cell : layers[0])
    {
        Cell2d cell2d;
        const Point2d & center = grid_vertices.center(cell.box);
        cell2d.x = center[0];
        cell2d.y = center[1];
        for (int vertex : grid_vertices.outline(cell.box, true))
        {
            auto found = vertex_numbers.find(vertex);
            if (found == vertex_numbers.end())
            {
                found = vertex_numbers.insert(make_pair(vertex, (int)data.vertices.size())).first;
                data.vertices.push_back(grid_vertices.vertices()[vertex]);
            }
            cell2d.vertices.push_back(found->second);
        }
        data.cell2d.push_back(cell2d);
    }

    //elevations of the cells, and of the surfaces of the base grid for an
    //inactive cell
    auto surface = [&](int index, const LayerCell & cell)
    { return mfgrid->botm(index, cell.row, cell.column); };
    auto cell_top = [&](int k, const LayerCell & cell)
    { return (qtree != nullptr) ? qtree->top[cell.box->number] : surface(k, cell); };
    auto cell_bottom = [&](int k, const LayerCell & cell)
    { return (qtree != nullptr) ? qtree->bot[cell.box->number] : surface(k + 1, cell); };

    for (int j = 0; j < data.ncpl; j++)
    {
        const LayerCell & cell = layers[0][j];
        data.top.push_back(cell.box->active ? cell_top(0, cell) : surface(0, cell));
    }
    for (int k = 0; k < data.nlay; k++)
    {
        for (int j = 0; j < data.ncpl; j++)
        {
            const LayerCell & cell = layers[k][j];
            double bottom;
            if (cell.box->active) bottom = cell_bottom(k, cell);
            //the top of an active cell below an inactive one is the bottom of the
            //inactive cell
            else if (k + 1 < data.nlay && layers[k + 1][j].box->active)
                bottom = cell_top(k + 1, layers[k + 1][j]);
            else
                bottom = surface(k + 1, cell);
            data.botm.push_back(bottom);
            data.idomain.push_back(cell.box->active ? 1 : (vertical_pass_through ? -1 : 0));

            if (k == 0 || !cell.box->active || !layers[k - 1][j].box->active) continue;
            double top = cell_top(k, cell);
            double above = cell_bottom(k - 1, layers[k - 1][j]);
            if (fabs(top - above) > 1.0e-9 * (fabs(top) + fabs(above) + 1.0))
                data.tops_not_aligned++;
        }
    }

    if (qtree == nullptr) return data;

    //horizontal ghost nodes by layer and cell; every layer has the same cells,
    //so no vertical connection needs one
    data.has_ghost_nodes = true;
    map<long, pair<int, int>> layer_cell;
    for (int k = 0; k < data.nlay; k++)
        for (int j = 0; j < data.ncpl; j++)
            if (layers[k][j].box->active) layer_cell[layers[k][j].box->number] = make_pair(k, j);
    vector<GhostNode> ghost_nodes;
    data.inexact_ghost_nodes = compute_ghost_nodes(qtree, ghost_nodes, false);
    for (const GhostNode & ghost_node : ghost_nodes)
    {
        DisvGhostNode disv_ghost_node;
        disv_ghost_node.layer = layer_cell.at(ghost_node.large_node).first;
        disv_ghost_node.large_cell = layer_cell.at(ghost_node.large_node).second;
        disv_ghost_node.small_cell = layer_cell.at(ghost_node.small_node).second;
        disv_ghost_node.contributor1 = layer_cell.at(ghost_node.j1).second;
        disv_ghost_node.contributor2 = layer_cell.at(ghost_node.j2).second;
        disv_ghost_node.alpha1 = ghost_node.alpha1;
        disv_ghost_node.alpha2 = ghost_node.alpha2;
        data.ghost_nodes.push_back(disv_ghost_node);
    }
    return data;
}

void write_disv_data(const DisvData & data, const string & prefix)
{
    vector<string> written;

    //the DIMENSIONS block
    string name = prefix + ".dimensions.dat";
    ofstream fout = open_output(name);
    fout << "NLAY " << data.nlay << "\n";
    fout << "NCPL " << data.ncpl << "\n";
    fout << "NVERT " << data.vertices.size() << "\n";
    fout.close();
    written.push_back(name);

    name = prefix + ".vertices.dat";
    fout = open_output(name);
    for (size_t i = 0; i < data.vertices.size(); i++)
        fout << i + 1 << " " << RoundTrip(data.vertices[i][0]) << " "
             << RoundTrip(data.vertices[i][1]) << "\n";
    fout.close();
    written.push_back(name);

    name = prefix + ".cell2d.dat";
    fout = open_output(name);
    for (size_t j = 0; j < data.cell2d.size(); j++)
    {
        const Cell2d & cell = data.cell2d[j];
        fout << j + 1 << " " << RoundTrip(cell.x) << " " << RoundTrip(cell.y) << " "
             << cell.vertices.size();
        for (int vertex : cell.vertices) fout << " " << vertex + 1;
        fout << "\n";
    }
    fout.close();
    written.push_back(name);

    name = prefix + ".top.dat";
    fout = open_output(name);
    write_layers(fout, data.top, data.ncpl);
    fout.close();
    written.push_back(name);

    name = prefix + ".botm.dat";
    fout = open_output(name);
    write_layers(fout, data.botm, data.ncpl);
    fout.close();
    written.push_back(name);

    name = prefix + ".idomain.dat";
    fout = open_output(name);
    write_layers(fout, data.idomain, data.ncpl);
    fout.close();
    written.push_back(name);

    //cellids of n, m, and the two contributing cells, then the two factors
    if (data.has_ghost_nodes)
    {
        name = prefix + ".gnc.dat";
        fout = open_output(name);
        for (const DisvGhostNode & ghost_node : data.ghost_nodes)
        {
            int layer = ghost_node.layer + 1;
            fout << layer << " " << ghost_node.large_cell + 1 << " " << layer << " "
                 << ghost_node.small_cell + 1 << " " << layer << " " << ghost_node.contributor1 + 1
                 << " " << layer << " " << ghost_node.contributor2 + 1 << " "
                 << RoundTrip(ghost_node.alpha1) << " " << RoundTrip(ghost_node.alpha2) << "\n";
        }
        fout.close();
        written.push_back(name);
    }

    for (const string & file : written) cout << "- Save data to file: " << file << endl;
}

} // namespace cusg
