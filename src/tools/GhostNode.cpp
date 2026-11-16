#include "GhostNode.h"
#include "GridgenError.h"
#include "RoundTrip.h"
#include <algorithm>
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

//roundoff allowed in a contributing factor, so that a factor that is zero in
//exact arithmetic is not rejected as negative and equal sums tie
constexpr double FACTOR_TOLERANCE = 1.0e-9;

//true if the contributing factors are not negative and sum to less than one,
//within roundoff; a factor that is negative by roundoff only is set to zero
bool interpolates(double & alpha1, double & alpha2)
{
    if (alpha1 < -FACTOR_TOLERANCE || alpha2 < -FACTOR_TOLERANCE) return false;
    if (alpha1 + alpha2 >= 1.0 - FACTOR_TOLERANCE) return false;
    alpha1 = max(alpha1, 0.0);
    alpha2 = max(alpha2, 0.0);
    return true;
}

//true if contributor alone reproduces a linear head at the ghost node offset
//by (dx, dy) from the center of large_cell, as when it is in line with them
bool in_line(Box * large_cell, Box * contributor, double dx, double dy)
{
    double cross = dx * (contributor->y - large_cell->y) - dy * (contributor->x - large_cell->x);
    return fabs(cross) <= 1.0e-9 * fabs(large_cell->dx * large_cell->dy);
}

//interpolate the head at the ghost node, offset by (dx, dy) from the center
//of large_cell, from contributor and a second neighbor of large_cell when
//contributor is offset along the connection, as a coarser neighbor is; of
//the neighbors that give contributing factors that are not negative and sum
//to less than one, the one whose factors sum least. False, and unchanged, if
//no neighbor does; true if contributor is in line
bool add_second_contributor(Box * small_cell, Box * large_cell, Box * contributor, double dx,
                            double dy, GhostNode & ghost_node)
{
    if (in_line(large_cell, contributor, dx, dy)) return true;

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
            if (!interpolates(alpha1, alpha2)) continue;
            if (best == nullptr || alpha1 + alpha2 < best_alpha1 + best_alpha2 - FACTOR_TOLERANCE)
            {
                best = candidate;
                best_alpha1 = alpha1;
                best_alpha2 = alpha2;
            }
        }
    }
    if (best == nullptr) return false;

    ghost_node.j1 = contributor->number;
    ghost_node.j2 = best->number;
    ghost_node.alpha1 = best_alpha1;
    ghost_node.alpha2 = best_alpha2;
    return true;
}

//add the active neighbors of cell, other than skipped and those already in
//candidates, to candidates
void add_active_neighbors(Box * cell, Box * skipped, vector<Box *> & candidates)
{
    for (int side = 0; side < 4; side++)
    {
        list<Box *> neighbors;
        cell->get_neighbors(neighbors, side);
        for (Box * neighbor : neighbors)
        {
            if (neighbor->number == -1 || neighbor == skipped) continue;
            if (find(candidates.begin(), candidates.end(), neighbor) != candidates.end()) continue;
            candidates.push_back(neighbor);
        }
    }
}

//the pair of candidates whose contributing factors give the head of a linear
//head field at the ghost node offset by (dx, dy) from the center of
//large_cell, are not negative, sum to less than one, and sum least; false if
//no pair does
bool best_pair(Box * large_cell, const vector<Box *> & candidates, double dx, double dy,
               GhostNode & ghost_node)
{
    bool found = false;
    double best_sum = 0.0;
    for (size_t i = 0; i < candidates.size(); i++)
    {
        for (size_t k = i + 1; k < candidates.size(); k++)
        {
            double alpha1, alpha2;
            if (!barycentric_weights(large_cell, candidates[i], candidates[k], dx, dy, alpha1,
                                     alpha2))
                continue;
            if (!interpolates(alpha1, alpha2)) continue;
            if (found && alpha1 + alpha2 >= best_sum - FACTOR_TOLERANCE) continue;
            found = true;
            best_sum = alpha1 + alpha2;
            ghost_node.j1 = candidates[i]->number;
            ghost_node.j2 = candidates[k]->number;
            ghost_node.alpha1 = alpha1;
            ghost_node.alpha2 = alpha2;
        }
    }
    return found;
}

//interpolate the head at the ghost node, offset by (dx, dy) from the center of
//large_cell, where the cells beside large_cell on the side of the offset are
//inactive, missing, or coarser, so that a linear head field is reproduced: from
//two active neighbors of large_cell, or small_cell when use_small_cell is true,
//and otherwise from the neighbors of those neighbors too. The head is not
//projected beyond the active cells, which would take a negative contributing
//factor; false if no pair interpolates it
bool interpolated_contributors(Box * small_cell, Box * large_cell, double dx, double dy,
                               bool use_small_cell, GhostNode & ghost_node)
{
    //small_cell contributes to a horizontal connection only; directly above or
    //below a vertical ghost node, its factor would be one and the connection
    //would carry no flow
    Box * skipped = use_small_cell ? large_cell : small_cell;
    vector<Box *> candidates;
    add_active_neighbors(large_cell, skipped, candidates);
    if (best_pair(large_cell, candidates, dx, dy, ghost_node)) return true;

    //a cell with one active neighbor, or neighbors in line with it
    size_t first_ring = candidates.size();
    for (size_t i = 0; i < first_ring; i++)
        add_active_neighbors(candidates[i], skipped, candidates);
    candidates.erase(remove(candidates.begin(), candidates.end(), large_cell), candidates.end());
    return best_pair(large_cell, candidates, dx, dy, ghost_node);
}

//true if cell has an inactive neighbor in its layer, at the edge of an active domain
bool next_to_inactive_cell(Box * cell)
{
    for (int side = 0; side < 4; side++)
    {
        list<Box *> neighbors;
        cell->get_neighbors(neighbors, side);
        for (Box * neighbor : neighbors)
        {
            if (neighbor->number == -1) return true;
        }
    }
    return false;
}

//add the ghost node of the horizontal connection between large_cell and the
//smaller small_cell; inexact_at_domain_edge counts the connections at the edge
//of an active domain whose ghost node does not reproduce a linear head field
void compute_gn(Box * small_cell, Box * large_cell, int direction, vector<GhostNode> & ghost_nodes,
                int & inexact_at_domain_edge)
{ //decide the node j first

    GhostNode ghost_node;
    ghost_node.small_node = small_cell->number;
    ghost_node.large_node = large_cell->number;
    //the cell the ghost node is interpolated from when there is one
    Box * contributor = nullptr;

    //the ghost node is level with the center of large_cell across the
    //connection and with the center of small_cell along it
    bool north_south = (direction == 0 || direction == 2);
    double dx = north_south ? small_cell->x - large_cell->x : 0.0;
    double dy = north_south ? 0.0 : small_cell->y - large_cell->y;

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

        if (neighbors.empty())
        {
            if (interpolated_contributors(small_cell, large_cell, dx, dy, true, ghost_node))
                ghost_nodes.push_back(ghost_node);
            else if (next_to_inactive_cell(large_cell))
                inexact_at_domain_edge++;
            return;
        }

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

        if (neighbors.empty())
        {
            if (interpolated_contributors(small_cell, large_cell, dx, dy, true, ghost_node))
                ghost_nodes.push_back(ghost_node);
            else if (next_to_inactive_cell(large_cell))
                inexact_at_domain_edge++;
            return;
        }

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

    //an inactive neighbor on the side of the offset
    if (ghost_node.j1 == -1 || ghost_node.j2 == -1)
    {
        if (interpolated_contributors(small_cell, large_cell, dx, dy, true, ghost_node))
            ghost_nodes.push_back(ghost_node);
        else
            inexact_at_domain_edge++;
        return;
    }

    //a coarser neighbor that, with no second neighbor, does not reproduce a
    //linear head field; the ghost node keeps the coarser neighbor alone if no
    //pair interpolates it
    if (contributor != nullptr &&
        !add_second_contributor(small_cell, large_cell, contributor, dx, dy, ghost_node))
    {
        GhostNode interpolated = ghost_node;
        if (interpolated_contributors(small_cell, large_cell, dx, dy, true, interpolated))
            ghost_node = interpolated;
        else if (next_to_inactive_cell(large_cell))
            inexact_at_domain_edge++;
    }

    ghost_nodes.push_back(ghost_node);
}
void compute_gn_neigh(Box * small_cell, list<Box *> & neighbors, int direction,
                      vector<GhostNode> & ghost_nodes, int & inexact_at_domain_edge)
{
    //neighbors
    for (auto neighbor : neighbors)
    {
        if (neighbor->number == -1) continue;

        if (neighbor->depth >= small_cell->depth) continue;

        //the node should be considered
        compute_gn(small_cell, neighbor, direction, ghost_nodes, inexact_at_domain_edge);
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
//small_cell in x and y; inexact_at_domain_edge counts the connections at the
//edge of an active domain whose ghost node does not reproduce a linear head field
void compute_vertical_gn(Box * small_cell, Box * large_cell, vector<GhostNode> & ghost_nodes,
                         int & inexact_at_domain_edge)
{
    double dx = small_cell->x - large_cell->x;
    double dy = small_cell->y - large_cell->y;

    //directions 0 to 3 are north, east, south, and west
    Box * x_neighbor = side_neighbor(large_cell, (dx > 0) ? 1 : 3, dy);
    Box * y_neighbor = side_neighbor(large_cell, (dy > 0) ? 0 : 2, dx);
    GhostNode ghost_node;
    ghost_node.large_node = large_cell->number;
    ghost_node.small_node = small_cell->number;

    //a missing neighbor in x or y, as at the edge of an active domain or of the
    //grid; if no pair interpolates the ghost node, it is interpolated in one
    //direction only, or not at all
    if (x_neighbor == nullptr || y_neighbor == nullptr)
    {
        if (interpolated_contributors(small_cell, large_cell, dx, dy, false, ghost_node))
        {
            ghost_nodes.push_back(ghost_node);
            return;
        }
        if (next_to_inactive_cell(large_cell)) inexact_at_domain_edge++;
    }
    if (x_neighbor == nullptr && y_neighbor == nullptr) return;

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

int compute_ghost_nodes(QuadTree3D * qtree, vector<GhostNode> & ghost_nodes, bool vertical)
{
    int inexact_at_domain_edge = 0;
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
            compute_gn_neigh(b, nei[j], j, ghost_nodes, inexact_at_domain_edge);
        }
    }

    if (!vertical) return inexact_at_domain_edge;

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
            compute_vertical_gn(neighbor, b, ghost_nodes, inexact_at_domain_edge);
        }
    }
    return inexact_at_domain_edge;
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