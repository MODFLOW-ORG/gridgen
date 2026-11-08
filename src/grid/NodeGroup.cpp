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

    Whenever nodes are added or removed, the rebuild_box_index_array flag
    is set to true indicating that the nodes need to be renumbered.

*/

//constructor
NodeGroup::NodeGroup(Grid * g)
{
    grid = g;
    nlay = 1; //number of layers (? why is this here)
    rebuild_box_index_array = true;
    nextnodenumber = 0;
    chunksize = 10;
    total_leaf_number = 0;
    active_leaf_number = 0;
}

int NodeGroup::add(Box * box)
{
    int box_index = 0;

    if (removed_nodes.empty() == false)
    {
        box_index = removed_nodes.front();
        removed_nodes.pop_front();
        (*this)[box_index] = box;
    }
    else
    {
        box_index = size();
        push_back(box);
    }

    owned_boxes.emplace_back(box);
    box->id = box_index; //box_index is the position of the box in the array (group)

    rebuild_box_index_array = true;

    if (box->isLeaf)
    {
        total_leaf_number++;
        if (box->active) active_leaf_number += 1;
    }

    if (box->layer > nlay) nlay = box->layer;

    return box_index;
}

void NodeGroup::share(const NodeGroup & other)
{
    vector<Box *>::operator=(other);
    owned_boxes.clear();
    grid = other.grid;
    box_index_array = other.box_index_array;
    layer_cell_counts = other.layer_cell_counts;
    nlay = other.nlay;
    rebuild_box_index_array = other.rebuild_box_index_array;
    nextnodenumber = other.nextnodenumber;
    chunksize = other.chunksize;
    removed_nodes = other.removed_nodes;
    total_leaf_number = other.total_leaf_number;
    active_leaf_number = other.active_leaf_number;
}

/*!
    Return an array of size of the number of active leaf nodes
    and include in the array the box_index.

    box_index_array[cell] = box_index
*/
int * NodeGroup::get_box_index_array()
{
    if (rebuild_box_index_array) build_node_info();
    return box_index_array.data();
}

/*!
 * return a layer_cell_counts, which contains the number of
    nodes per layer:
    layer_cell_counts[layer] = nodes in layer
    This method numbers the nodes by layer, starting with the
    lowest layer number.
 *
 */
int * NodeGroup::get_layer_cell_counts()
{
    if (rebuild_box_index_array) build_node_info();
    return layer_cell_counts.data();
}

///number of (active) leaf nodes in this group
NodeGroup::size_type NodeGroup::node_size()
{
    if (rebuild_box_index_array) build_node_info();
    return active_leaf_number;
}

///get the total number of leaf size including inactive leaves
NodeGroup::size_type NodeGroup::get_total_leaf_size()
{
    if (rebuild_box_index_array) build_node_info();
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
 * build information of layer_cell_counts and box_index_array
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
    box_index_array.assign(this->active_leaf_number, 0);
    layer_cell_counts.assign(nlay, 0); // number of nodes per layer

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
                box_index_array[box->number] = box->id;
                layer_cell_counts[box->layer] += 1;
            }
        }
    } //end i

    rebuild_box_index_array = false; //set the flag since the array was just built
}

} // end namespace cusg
