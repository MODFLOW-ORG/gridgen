#include "VertexExport.h"

#include <algorithm>
#include <fstream>
#include <iostream>
#include <map>
#include <vector>

#include "GridToVTK.h"
#include "GridVertices.h"
#include "GridgenError.h"
#include "ModflowGrid.h"
#include "QuadTree3D.h"
#include "RoundTrip.h"

namespace cusg
{

namespace
{

//active cells in node number order
vector<Box *> active_cells(Grid * grid)
{
    vector<pair<long, Box *>> numbered;
    for (auto box : grid->node_group)
        if (box->active && box->is_leaf) numbered.push_back(make_pair(box->number, box));
    sort(numbered.begin(), numbered.end());
    vector<Box *> cells;
    for (auto & entry : numbered) cells.push_back(entry.second);
    return cells;
}

ModflowGrid * base_grid(Grid * grid)
{
    QuadTree3D * qtree = dynamic_cast<QuadTree3D *>(grid);
    if (qtree != nullptr) return qtree->get_modflow_grid();
    ModflowGrid * mfgrid = dynamic_cast<ModflowGrid *>(grid);
    if (mfgrid == nullptr)
        throw GridgenError("the grid to export is neither a MODFLOW grid nor a quadtree grid");
    return mfgrid;
}

ofstream open_output(const string & filename)
{
    ofstream fout(filename.c_str());
    if (fout.good() == false) throw GridgenError("cannot write file " + filename);
    return fout;
}

void collect_base_cells(Box * box, int layer, int row, int column, map<Box *, array<int, 3>> & base)
{
    if (box->is_leaf)
    {
        base[box] = {layer, row, column};
        return;
    }
    for (auto child : box->p_children) collect_base_cells(child, layer, row, column, base);
}

} // namespace

void write_gsf(Grid * grid, const string & filename, bool hanging_vertices)
{
    ModflowGrid * mfgrid = base_grid(grid);
    QuadTree3D * qtree = dynamic_cast<QuadTree3D *>(grid);
    if (qtree != nullptr) qtree->number_nodes();

    GridVertices gv(grid);
    map<pair<int, int>, double> elevations = vertex_elevations(grid, gv);
    vector<Box *> cells = active_cells(grid);

    //vertices at the top and bottom of each cell, numbered in the order the
    //cells use them
    map<pair<int, int>, int> vertex_numbers;
    vector<pair<int, int>> vertices;
    vector<vector<int>> cell_vertices(cells.size());
    for (size_t n = 0; n < cells.size(); n++)
    {
        Box * box = cells[n];
        vector<int> ring;
        if (hanging_vertices)
            ring = gv.outline(box);
        else
            ring.assign(gv.corners(box), gv.corners(box) + 4);
        for (int surface = box->layer; surface <= box->layer + 1; surface++)
        {
            for (int vertex : ring)
            {
                pair<int, int> key(vertex, surface);
                auto found = vertex_numbers.find(key);
                if (found == vertex_numbers.end())
                {
                    found = vertex_numbers.insert(make_pair(key, (int)vertices.size() + 1)).first;
                    vertices.push_back(key);
                }
                cell_vertices[n].push_back(found->second);
            }
        }
    }

    ofstream fout = open_output(filename);
    fout << "UNSTRUCTURED GWF\n";
    fout << cells.size() << " " << mfgrid->nlay << " 1 1\n";
    fout << vertices.size() << "\n";
    for (auto & key : vertices)
    {
        const Point2d & xy = gv.vertices()[key.first];
        fout << RoundTrip(xy[0]) << " " << RoundTrip(xy[1]) << " " << RoundTrip(elevations.at(key))
             << "\n";
    }
    for (size_t n = 0; n < cells.size(); n++)
    {
        Box * box = cells[n];
        const Point2d & center = gv.center(box);
        fout << box->number + 1 << " " << RoundTrip(center[0]) << " " << RoundTrip(center[1]) << " "
             << RoundTrip(box->z) << " " << box->layer + 1 << " " << cell_vertices[n].size();
        for (int vertex : cell_vertices[n]) fout << " " << vertex;
        fout << "\n";
    }
    fout.close();
    cout << "- Save grid specification file: " << filename << endl;
}

void write_node_data_table(Grid * grid, const string & filename)
{
    ModflowGrid * mfgrid = base_grid(grid);
    map<Box *, array<int, 3>> base;
    for (int k = 0; k < mfgrid->nlay; k++)
        for (int r = 0; r < mfgrid->nrow; r++)
            for (int c = 0; c < mfgrid->ncol; c++)
                collect_base_cells(grid->node_group[mfgrid->get_box_index(k, r, c)], k, r, c, base);

    ofstream fout = open_output(filename);
    for (auto box : active_cells(grid))
    {
        const array<int, 3> & cell = base.at(box);
        fout << box->number + 1 << " " << cell[0] + 1 << " " << cell[1] + 1 << " " << cell[2] + 1
             << "\n";
    }
    fout.close();
    cout << "- Save node data table: " << filename << endl;
}

void write_vertices_cell2d(Grid * grid, const string & prefix)
{
    base_grid(grid);
    GridVertices gv(grid);
    vector<Box *> cells = active_cells(grid);

    //vertices numbered in the order the cells use them
    map<int, int> vertex_numbers;
    vector<int> vertices;
    vector<vector<int>> cell_vertices(cells.size());
    for (size_t n = 0; n < cells.size(); n++)
    {
        for (int vertex : gv.outline(cells[n]))
        {
            auto found = vertex_numbers.find(vertex);
            if (found == vertex_numbers.end())
            {
                found = vertex_numbers.insert(make_pair(vertex, (int)vertices.size() + 1)).first;
                vertices.push_back(vertex);
            }
            cell_vertices[n].push_back(found->second);
        }
    }

    string vertices_name = prefix + ".vertices.dat";
    ofstream fout = open_output(vertices_name);
    for (size_t i = 0; i < vertices.size(); i++)
    {
        const Point2d & xy = gv.vertices()[vertices[i]];
        fout << i + 1 << " " << RoundTrip(xy[0]) << " " << RoundTrip(xy[1]) << "\n";
    }
    fout.close();
    cout << "- Save data to file: " << vertices_name << endl;

    string cell2d_name = prefix + ".cell2d.dat";
    fout = open_output(cell2d_name);
    for (size_t n = 0; n < cells.size(); n++)
    {
        const Point2d & center = gv.center(cells[n]);
        fout << cells[n]->number + 1 << " " << RoundTrip(center[0]) << " " << RoundTrip(center[1])
             << " " << cell_vertices[n].size();
        for (int vertex : cell_vertices[n]) fout << " " << vertex;
        fout << "\n";
    }
    fout.close();
    cout << "- Save data to file: " << cell2d_name << endl;
}

} // namespace cusg
