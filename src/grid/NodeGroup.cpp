/*! \file NodeGroup.cpp
 
\brief Container for storing references to `Box` nodes
\author <a href="http://masc.cs.gmu.edu">MASC Group</a>,
        George Mason University
\author <a href="http://profile.usgs.gov/langevin">Christian Langevin</a>,
        USGS
\bug    No known bugs.

*/

//   $Id: $

#include "NodeGroup.h"
#include "Grid.h"

#include <cassert>
using namespace std;

namespace cusg
{

/*!

    This object emulates a container by storing references to all node
    objects.  This object manages all of the node objects and stores
    reference to them in the nodearray, which is a numpy array.

    Added a removed nodes capability so that when a node is removed,
    the space in the nodegroup array can be reused.  The add_node
    method first checks to see if there are any removed node spaces
    that can be used and it uses them if possible.

    Whenever nodes are added or removed, the rebuild_nodeid_array flag
    is set to true indicating that the nodes need to be renumbered.

*/

//constructor
NodeGroup::NodeGroup(Grid * g)
{
    grid = g;
    nlay = 1; //number of layers (? why is this here)
    rebuild_nodeid_array = true;
    nextnodenumber = 0;
    chunksize = 10;
    total_leaf_number = 0;
    active_leaf_number = 0;
}

int NodeGroup::add(Box * box)
{
    int nodeid = 0;

    if (removed_nodes.empty() == false)
    {
        nodeid = removed_nodes.front();
        removed_nodes.pop_front();
        (*this)[nodeid] = box;
    }
    else
    {
        nodeid = size();
        push_back(box);
    }

    owned_boxes.emplace_back(box);
    box->id = nodeid; //nodeid is the position of the box in the array (group)

    rebuild_nodeid_array = true;

    if (box->isLeaf)
    {
        total_leaf_number++;
        if (box->active) active_leaf_number += 1;
    }

    if (box->layer > nlay) nlay = box->layer;

    return nodeid;
}

void NodeGroup::share(const NodeGroup & other)
{
    vector<Box *>::operator=(other);
    owned_boxes.clear();
    grid = other.grid;
    nodeid_array = other.nodeid_array;
    nodelay_array = other.nodelay_array;
    nlay = other.nlay;
    rebuild_nodeid_array = other.rebuild_nodeid_array;
    nextnodenumber = other.nextnodenumber;
    chunksize = other.chunksize;
    removed_nodes = other.removed_nodes;
    total_leaf_number = other.total_leaf_number;
    active_leaf_number = other.active_leaf_number;
}

/*!
    Return an array of size of the number of active leaf nodes
    and include in the array the nodeid.

    nodeid_array[nodenumber] = nodeid
*/
int * NodeGroup::get_nodeid_array()
{
    if (rebuild_nodeid_array) build_node_info();
    return nodeid_array.data();
}

/*!
 * return a nodelay_array, which contains the number of
    nodes per layer:
    nodelay_array[layer] = nodes in layer
    This method numbers the nodes by layer, starting with the
    lowest layer number.
 *
 */
int * NodeGroup::get_nodelay_array()
{
    if (rebuild_nodeid_array) build_node_info();
    return nodelay_array.data();
}

///number of (active) leaf nodes in this group
NodeGroup::size_type NodeGroup::node_size()
{
    if (rebuild_nodeid_array) build_node_info();
    return active_leaf_number;
}

///get the total number of leaf size including inactive leaves
NodeGroup::size_type NodeGroup::get_total_leaf_size()
{
    if (rebuild_nodeid_array) build_node_info();
    return total_leaf_number;
}

void NodeGroup::disable_nodes_in_layer(int layer)
{
    for (auto box : *this)
    {
        if (box == nullptr) continue;
        if (box->layer != layer) continue;
        box->active = false;
    }
}

/*!
 * build information of nodelay_array and nodeid_array
 *
 */
void NodeGroup::build_node_info()
{

    //count the number of leaf nodes
    this->total_leaf_number = 0;
    this->active_leaf_number = 0;
    for (auto & i : *this)
    {
        if (i == nullptr) continue;
        if (i->isLeaf)
        {
            total_leaf_number++;
            if (i->active) active_leaf_number += 1;
        }
    }

    //reset things...
    nodeid_array.assign(this->active_leaf_number, 0);
    nodelay_array.assign(nlay, 0); // number of nodes per layer

    //number the nodes, starting with the lowest layer number
    //this is the id of the node in CSR/COO
    //only leaves have valid node number
    grid->number_nodes();

    //node number is ordered by layer
    //node id is the position in this class
    for (auto box : *this)
    {
        if (box == nullptr) continue;

        if (box->isLeaf)
        {
            if (box->active)
            {
                nodeid_array[box->number] = box->id;
                nodelay_array[box->layer] += 1;
            }
        }
    } //end i

    rebuild_nodeid_array = false; //set the flag since the array was just built
}

} // end namespace cusg
