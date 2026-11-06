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
void computeGN(Box * small_cell, Box * large_cell, int direction, vector<GhostNode> & ghost_nodes)
{ //decide the node j first

    GhostNode ghost_node;
    ghost_node.small_node = small_cell->number;
    ghost_node.large_node = large_cell->number;

    //if it is in north or south
    if (direction == 0 || direction == 2)
    {
        list<Box *> neighbors;
        if (between(large_cell->x - large_cell->dx, large_cell->x, small_cell->x))
        { //search the 3 direction(west)
            large_cell->getNeighbors(neighbors, 3);
        }
        else if (between(large_cell->x, large_cell->x + large_cell->dx, small_cell->x))
        { //search for the 1 direction
            large_cell->getNeighbors(neighbors, 1);
        }

        if (neighbors.empty()) return;

        if (neighbors.size() > 2)
        {
            throw GridgenError("the node " + to_string(large_cell->id) + " has more than two (" +
                               to_string(neighbors.size()) +
                               ") neighbor box nodes in direction=" + to_string(direction));
        }

        if (neighbors.size() == 1)
        {
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
            large_cell->getNeighbors(neighbors, 0);
        else if (between(large_cell->y, large_cell->y - large_cell->dy, small_cell->y))
            large_cell->getNeighbors(neighbors, 2);

        if (neighbors.empty()) return;

        if (neighbors.size() > 2)
        {
            throw GridgenError("the node " + to_string(large_cell->id) + " has more than two (" +
                               to_string(neighbors.size()) +
                               ") neighbor box nodes in direction=" + to_string(direction));
        }

        if (neighbors.size() == 1)
        {
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

    ghost_nodes.push_back(ghost_node);
}
void computeGNNeigh(Box * small_cell, list<Box *> & neighbors, int direction,
                    vector<GhostNode> & ghost_nodes)
{
    //neighbors
    for (auto neighbor : neighbors)
    {
        if (neighbor->number == -1) continue;

        if (neighbor->depth >= small_cell->depth) continue;

        //the node should be considered
        computeGN(small_cell, neighbor, direction, ghost_nodes);
    }
}

void computeGNC(QuadTree3D * qtree, vector<GhostNode> & ghost_nodes)
{
    //loop over all the leaf node
    for (auto b : qtree->nodegroup)
    {
        if (!b->isLeaf || b->number == -1) continue;

        list<Box *> nei[4];
        for (int i = 0; i < 4; i++)
        {
            b->getNeighbors(nei[i], i);
        }

        for (int j = 0; j < 4; j++)
        {
            computeGNNeigh(b, nei[j], j, ghost_nodes);
        }
    }
}

void writeGhostNode(string filePath, vector<GhostNode> & ghost_nodes, int offset = 1)
{
    ofstream ofile;
    ofile.open(filePath.c_str(), std::ofstream::out | std::ofstream::trunc);
    for (auto & ghost_node : ghost_nodes)
    {
        ofile << ghost_node.large_node + offset << "\t" << ghost_node.small_node + offset << "\t"
              << ghost_node.j1 + offset << "\t" << ghost_node.j2 + offset << "\t"
              << round_trip(ghost_node.alpha1) << "\t" << round_trip(ghost_node.alpha2) << endl;
    }
    ofile.close();
}
}