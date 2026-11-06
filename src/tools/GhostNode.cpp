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
void computeGN(Box * small_cell, Box * large_cell, int dir, vector<GhostNode> & gns)
{ //decide the node j first

    GhostNode gn;
    gn.small_node = small_cell->number;
    gn.large_node = large_cell->number;

    //if it is in north or south
    if (dir == 0 || dir == 2)
    {
        list<Box *> nl;
        if (between(large_cell->x - large_cell->dx, large_cell->x, small_cell->x))
        { //search the 3 direction(west)
            large_cell->getNeighbors(nl, 3);
        }
        else if (between(large_cell->x, large_cell->x + large_cell->dx, small_cell->x))
        { //search for the 1 direction
            large_cell->getNeighbors(nl, 1);
        }

        if (nl.empty()) return;

        if (nl.size() > 2)
        {
            throw GridgenError("the node " + to_string(large_cell->id) + " has more than two (" +
                               to_string(nl.size()) +
                               ") neighbor box nodes in dir=" + to_string(dir));
        }

        if (nl.size() == 1)
        {
            gn.j1 = gn.j2 = nl.front()->number;
            gn.alpha1 = ((large_cell->x - small_cell->x) / (large_cell->x - nl.front()->x)) /
                        2; // fabs((large_cell->x - small_cell->x)/large_cell->dx);
            gn.alpha2 = gn.alpha1;
        }
        else
        {
            assert(nl.size() == 2);
            gn.j1 = nl.front()->number;
            gn.j2 = nl.back()->number;
            //they have the same X coordinates
            gn.alpha1 = ((large_cell->x - small_cell->x) / (large_cell->x - nl.front()->x)) /
                        2; //fabs( ( (large_cell->x - small_cell->x)/large_cell->dx ) / 2 );
            gn.alpha2 = gn.alpha1;
        }
    }
    else if (dir == 1 || dir == 3)
    {
        list<Box *> nl;
        if (between(large_cell->y + large_cell->dy, large_cell->y, small_cell->y))
            large_cell->getNeighbors(nl, 0);
        else if (between(large_cell->y, large_cell->y - large_cell->dy, small_cell->y))
            large_cell->getNeighbors(nl, 2);

        if (nl.empty()) return;

        if (nl.size() > 2)
        {
            throw GridgenError("the node " + to_string(large_cell->id) + " has more than two (" +
                               to_string(nl.size()) +
                               ") neighbor box nodes in dir=" + to_string(dir));
        }

        if (nl.size() == 1)
        {
            gn.j1 = gn.j2 = nl.front()->number;
            gn.alpha1 = gn.alpha2 =
                ((large_cell->y - small_cell->y) / (large_cell->y - nl.front()->y)) /
                2; //fabs( (large_cell->y - small_cell->y) / large_cell->dy);
        }
        else
        {
            assert(nl.size() == 2);
            gn.j1 = nl.front()->number;
            gn.j2 = nl.back()->number;
            gn.alpha1 = gn.alpha2 =
                ((large_cell->y - small_cell->y) / (large_cell->y - nl.front()->y)) /
                2; //fabs( ( (large_cell->y - small_cell->y) / large_cell->dy) / 2 );
        }
    }
    else
        assert(false);

    if (gn.j1 == -1 || gn.j2 == -1) return; //don't add it to gns if j1 or j2 is inactive...

    gns.push_back(gn);
}
void computeGNNeigh(Box * small_cell, list<Box *> & nlist, int dir, vector<GhostNode> & gns)
{
    //neighbors
    for (auto tb : nlist)
    {
        if (tb->number == -1) continue;

        if (tb->depth >= small_cell->depth) continue;

        //the node should be considered
        computeGN(small_cell, tb, dir, gns);
    }
}

void computeGNC(QuadTree3D * qtree, vector<GhostNode> & gns)
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
            computeGNNeigh(b, nei[j], j, gns);
        }
    }
}

void writeGhostNode(string filePath, vector<GhostNode> & gns, int offset = 1)
{
    ofstream ofile;
    ofile.open(filePath.c_str(), std::ofstream::out | std::ofstream::trunc);
    for (auto & ghn : gns)
    {
        ofile << ghn.large_node + offset << "\t" << ghn.small_node + offset << "\t"
              << ghn.j1 + offset << "\t" << ghn.j2 + offset << "\t" << round_trip(ghn.alpha1)
              << "\t" << round_trip(ghn.alpha2) << endl;
    }
    ofile.close();
}
}