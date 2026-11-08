/*! \file NodeGroup.h
 
\brief Container for storing references to `Box` nodes
\author <a href="http://masc.cs.gmu.edu">MASC Group</a>,
        George Mason University
\author <a href="http://profile.usgs.gov/langevin">Christian Langevin</a>,
        USGS
\bug    No known bugs.

*/

//   $Id: $

#pragma once

#include <memory>
#include <vector>
using namespace std;

#include "Box.h"

namespace cusg
{

class Grid; //defined in Grid.h

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
class NodeGroup : public vector<Box *>
{
public:
    //constructor
    NodeGroup(Grid * g = nullptr);

    NodeGroup(const NodeGroup &) = delete;
    NodeGroup & operator=(const NodeGroup &) = delete;
    NodeGroup(NodeGroup &&) = default;
    NodeGroup & operator=(NodeGroup &&) = default;

    ///add a box, which this group then owns, and return its node id
    int add(Box * box);

    ///refer to the boxes of another group, which keeps ownership of them
    void share(const NodeGroup & other);

    /*!
            Return an array of size of the number of active leaf nodes
            and include in the array the box_index.

            box_index_array[cell] = box_index
        */
    int * get_box_index_array();

    /*!
         * return layer_cell_counts, which contains the number of
            nodes per layer:
            layer_cell_counts[layer] = nodes in layer
            This method numbers the nodes by layer, starting with the
            lowest layer number.
         *
         */
    int * get_layer_cell_counts();

    ///make all nodes in the given layer as disabled
    void disable_nodes_in_layer(int layer);

    ///number of active leaf nodes in this group
    size_type node_size();

    ///get the total number of leaf size including inactive leaves
    size_type get_total_leaf_size();

    //
    //
    // Data (TODO: move these to protected section)
    //
    //
    Grid * grid;                   //the hosting Grid
    vector<int> box_index_array;   //node ids (mapping node number to node id)
    vector<int> layer_cell_counts; //number of nodes per layer
    int nlay;                      //total number of layers

    bool rebuild_box_index_array; //rebuild node id

    int nextnodenumber;

    int chunksize; ///what is a chunk size??
    list<int> removed_nodes;

protected:
    vector<unique_ptr<Box>> owned_boxes; //the boxes added to this group

    int total_leaf_number;  //number of leaves store in this group (total number of leaves)
    int active_leaf_number; //number of leaves store in this group (active nodes only)

    /*!
         * build information of layer_cell_counts and box_index_array
         *
         */
    void build_node_info();
};

} // end namespace cusg
