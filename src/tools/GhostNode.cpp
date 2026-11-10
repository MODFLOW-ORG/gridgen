#include "GhostNode.h"
#include "GridgenError.h"
#include "RoundTrip.h"
#include <fstream>
using namespace std;

namespace cusg
{
//check if c is between [a, b]
bool between(double a, double b, double c)
{
    if ((a < c && c < b) || (a > c && c > b))
        return true;
    else
        return false;
}
//the barycentric weights of the point offset by (dx, dy) from the center of
//cell in the triangle of the centers of cell, first, and second, so a head
//that varies linearly is reproduced; false if the triangle has no area
bool barycentric_weights(Box * cell, Box * first, Box * second, double dx, double dy,
                         double & alpha1, double & alpha2)
{
    double x1 = first->x - cell->x;
    double y1 = first->y - cell->y;
    double x2 = second->x - cell->x;
    double y2 = second->y - cell->y;
    double determinant = x1 * y2 - x2 * y1;
    if (fabs(determinant) <= 1.0e-9 * fabs(cell->dx * cell->dy)) return false;
    alpha1 = (dx * y2 - dy * x2) / determinant;
    alpha2 = (x1 * dy - y1 * dx) / determinant;
    return true;
}

//interpolate the ghost node, offset by (dx, dy) from the center of
//large_cell, from contributor and a second neighbor of large_cell when
//contributor is offset along the connection, as a coarser neighbor is; of
//the neighbors that give weights that are not negative and sum to less than
//one, the one whose weights sum least. Unchanged if contributor is aligned or
//there is none
void add_second_contributor(Box * small_cell, Box * large_cell, Box * contributor, double dx,
                            double dy, GhostNode & ghost_node)
{
    double offset_x = contributor->x - large_cell->x;
    double offset_y = contributor->y - large_cell->y;
    double cross = dx * offset_y - dy * offset_x;
    if (fabs(cross) <= 1.0e-9 * fabs(large_cell->dx * large_cell->dy)) return;

    Box * best = nullptr;
    double best_alpha1 = 0.0, best_alpha2 = 0.0;
    for (int side = 0; side < 4; side++)
    {
        list<Box *> neighbors;
        large_cell->get_neighbors(neighbors, side);
        for (Box * candidate : neighbors)
        {
            if (candidate->number == -1 || candidate == small_cell || candidate == contributor)
                continue;
            double alpha1, alpha2;
            if (!barycentric_weights(large_cell, contributor, candidate, dx, dy, alpha1, alpha2))
                continue;
            if (alpha1 < 0 || alpha2 < 0 || alpha1 + alpha2 >= 1) continue;
            if (best == nullptr || alpha1 + alpha2 < best_alpha1 + best_alpha2)
            {
                best = candidate;
                best_alpha1 = alpha1;
                best_alpha2 = alpha2;
            }
        }
    }
    if (best == nullptr) return;

    ghost_node.j1 = contributor->number;
    ghost_node.j2 = best->number;
    ghost_node.alpha1 = best_alpha1;
    ghost_node.alpha2 = best_alpha2;
}

void compute_gn(Box * small_cell, Box * large_cell, int direction, vector<GhostNode> & ghost_nodes)
{ //decide the node j first

    GhostNode ghost_node;
    ghost_node.small_node = small_cell->number;
    ghost_node.large_node = large_cell->number;
    //the cell the ghost node is interpolated from when there is one
    Box * contributor = nullptr;

    //if it is in north or south
    if (direction == 0 || direction == 2)
    {
        list<Box *> neighbors;
        if (between(large_cell->x - large_cell->dx, large_cell->x, small_cell->x))
        { //search the 3 direction(west)
            large_cell->get_neighbors(neighbors, 3);
        }
        else if (between(large_cell->x, large_cell->x + large_cell->dx, small_cell->x))
        { //search for the 1 direction
            large_cell->get_neighbors(neighbors, 1);
        }

        if (neighbors.empty()) return;

        if (neighbors.size() > 2)
        {
            throw GridgenError("a cell in layer " + to_string(large_cell->layer + 1) +
                               " has more than two neighbors on one side, so its ghost nodes "
                               "cannot be computed; smooth the grid so that neighboring cells "
                               "differ by at most one level of refinement");
        }

        if (neighbors.size() == 1)
        {
            contributor = neighbors.front();
            ghost_node.j1 = ghost_node.j2 = neighbors.front()->number;
            ghost_node.alpha1 =
                ((large_cell->x - small_cell->x) / (large_cell->x - neighbors.front()->x)) /
                2; // fabs((large_cell->x - small_cell->x)/large_cell->dx);
            ghost_node.alpha2 = ghost_node.alpha1;
        }
        else
        {
            assert(neighbors.size() == 2);
            ghost_node.j1 = neighbors.front()->number;
            ghost_node.j2 = neighbors.back()->number;
            //they have the same X coordinates
            ghost_node.alpha1 =
                ((large_cell->x - small_cell->x) / (large_cell->x - neighbors.front()->x)) /
                2; //fabs( ( (large_cell->x - small_cell->x)/large_cell->dx ) / 2 );
            ghost_node.alpha2 = ghost_node.alpha1;
        }
    }
    else if (direction == 1 || direction == 3)
    {
        list<Box *> neighbors;
        if (between(large_cell->y + large_cell->dy, large_cell->y, small_cell->y))
            large_cell->get_neighbors(neighbors, 0);
        else if (between(large_cell->y, large_cell->y - large_cell->dy, small_cell->y))
            large_cell->get_neighbors(neighbors, 2);

        if (neighbors.empty()) return;

        if (neighbors.size() > 2)
        {
            throw GridgenError("a cell in layer " + to_string(large_cell->layer + 1) +
                               " has more than two neighbors on one side, so its ghost nodes "
                               "cannot be computed; smooth the grid so that neighboring cells "
                               "differ by at most one level of refinement");
        }

        if (neighbors.size() == 1)
        {
            contributor = neighbors.front();
            ghost_node.j1 = ghost_node.j2 = neighbors.front()->number;
            ghost_node.alpha1 = ghost_node.alpha2 =
                ((large_cell->y - small_cell->y) / (large_cell->y - neighbors.front()->y)) /
                2; //fabs( (large_cell->y - small_cell->y) / large_cell->dy);
        }
        else
        {
            assert(neighbors.size() == 2);
            ghost_node.j1 = neighbors.front()->number;
            ghost_node.j2 = neighbors.back()->number;
            ghost_node.alpha1 = ghost_node.alpha2 =
                ((large_cell->y - small_cell->y) / (large_cell->y - neighbors.front()->y)) /
                2; //fabs( ( (large_cell->y - small_cell->y) / large_cell->dy) / 2 );
        }
    }
    else
        assert(false);

    if (ghost_node.j1 == -1 || ghost_node.j2 == -1)
        return; //don't add it to ghost_nodes if j1 or j2 is inactive...

    if (contributor != nullptr)
    {
        //the ghost node is level with the center of large_cell across the
        //connection and with the center of small_cell along it
        bool north_south = (direction == 0 || direction == 2);
        double dx = north_south ? small_cell->x - large_cell->x : 0.0;
        double dy = north_south ? 0.0 : small_cell->y - large_cell->y;
        add_second_contributor(small_cell, large_cell, contributor, dx, dy, ghost_node);
    }

    ghost_nodes.push_back(ghost_node);
}
void compute_gn_neigh(Box * small_cell, list<Box *> & neighbors, int direction,
                      vector<GhostNode> & ghost_nodes)
{
    //neighbors
    for (auto neighbor : neighbors)
    {
        if (neighbor->number == -1) continue;

        if (neighbor->depth >= small_cell->depth) continue;

        //the node should be considered
        compute_gn(small_cell, neighbor, direction, ghost_nodes);
    }
}

//the active neighbor of cell on side direction that a ghost node offset from
//the center of cell by transverse_offset along that side interpolates from; of
//two finer neighbors, the one on the side of the offset; nullptr if none
Box * side_neighbor(Box * cell, int direction, double transverse_offset)
{
    list<Box *> neighbors;
    cell->get_neighbors(neighbors, direction);
    if (neighbors.size() > 2)
    {
        throw GridgenError("a cell in layer " + to_string(cell->layer + 1) +
                           " has more than two neighbors on one side, so its ghost nodes "
                           "cannot be computed; smooth the grid so that neighboring cells "
                           "differ by at most one level of refinement");
    }

    Box * neighbor = nullptr;
    if (neighbors.size() == 1)
        neighbor = neighbors.front();
    else
    {
        //north and south sides run along x, east and west sides along y
        bool along_x = (direction == 0 || direction == 2);
        for (Box * candidate : neighbors)
        {
            double offset = along_x ? candidate->x - cell->x : candidate->y - cell->y;
            if ((offset > 0) == (transverse_offset > 0)) neighbor = candidate;
        }
    }

    if (neighbor == nullptr || neighbor->number == -1) return nullptr;
    return neighbor;
}

//add the ghost node of the vertical connection between large_cell and the
//smaller small_cell, placed in large_cell below or above the center of
//small_cell and interpolated from large_cell and its neighbors on the side of
//small_cell in x and y
void compute_vertical_gn(Box * small_cell, Box * large_cell, vector<GhostNode> & ghost_nodes)
{
    double dx = small_cell->x - large_cell->x;
    double dy = small_cell->y - large_cell->y;

    //directions 0 to 3 are north, east, south, and west
    Box * x_neighbor = side_neighbor(large_cell, (dx > 0) ? 1 : 3, dy);
    Box * y_neighbor = side_neighbor(large_cell, (dy > 0) ? 0 : 2, dx);
    if (x_neighbor == nullptr && y_neighbor == nullptr) return;

    GhostNode ghost_node;
    ghost_node.large_node = large_cell->number;
    ghost_node.small_node = small_cell->number;

    if (x_neighbor != nullptr && y_neighbor != nullptr)
    {
        //neighbors in x and y are never in line with large_cell
        ghost_node.j1 = x_neighbor->number;
        ghost_node.j2 = y_neighbor->number;
        barycentric_weights(large_cell, x_neighbor, y_neighbor, dx, dy, ghost_node.alpha1,
                            ghost_node.alpha2);
    }
    else
    {
        //one direction only, written as one cell twice with half the factor
        //each, as for a horizontal connection
        Box * neighbor = (x_neighbor != nullptr) ? x_neighbor : y_neighbor;
        double alpha = (x_neighbor != nullptr) ? dx / (neighbor->x - large_cell->x)
                                               : dy / (neighbor->y - large_cell->y);
        ghost_node.j1 = ghost_node.j2 = neighbor->number;
        ghost_node.alpha1 = ghost_node.alpha2 = alpha / 2;
    }

    ghost_nodes.push_back(ghost_node);
}

void compute_ghost_nodes(QuadTree3D * qtree, vector<GhostNode> & ghost_nodes, bool vertical)
{
    //loop over all the leaf node
    for (auto b : qtree->node_group)
    {
        if (!b->is_leaf || b->number == -1) continue;

        list<Box *> nei[4];
        for (int i = 0; i < 4; i++)
        {
            b->get_neighbors(nei[i], i);
        }

        for (int j = 0; j < 4; j++)
        {
            compute_gn_neigh(b, nei[j], j, ghost_nodes);
        }
    }

    if (!vertical) return;

    //after the horizontal ghost nodes, so they are written as without; the
    //vertical connections are those of the CSR data, including the vertical
    //pass-through connections across inactive cells
    for (auto b : qtree->node_group)
    {
        if (!b->is_leaf || b->number == -1) continue;

        list<pair<int, int>> connections =
            qtree->vertical_pass_through ? qtree->get_node_connections_with_pass_through(b, nullptr)
                                         : qtree->get_node_connections(b);
        for (const pair<int, int> & connection : connections)
        {
            //fldir is +3 for up and -3 for down
            if (abs(connection.second) != 3) continue;
            Box * neighbor = qtree->node_group[connection.first];
            if (neighbor->number == -1 || neighbor->depth <= b->depth) continue;
            compute_vertical_gn(neighbor, b, ghost_nodes);
        }
    }
}

void write_ghost_node(string file_path, vector<GhostNode> & ghost_nodes, int offset = 1)
{
    ofstream ofile;
    ofile.open(file_path.c_str(), std::ofstream::out | std::ofstream::trunc);
    for (auto & ghost_node : ghost_nodes)
    {
        ofile << ghost_node.large_node + offset << "\t" << ghost_node.small_node + offset << "\t"
              << ghost_node.j1 + offset << "\t" << ghost_node.j2 + offset << "\t"
              << RoundTrip(ghost_node.alpha1) << "\t" << RoundTrip(ghost_node.alpha2) << endl;
    }
    ofile.close();
}
}