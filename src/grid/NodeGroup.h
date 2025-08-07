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
    
        Whenever nodes are added or removed, the rebuild_nodeid_array flag
        is set to true indicating that the nodes need to be renumbered.

    */
    class NodeGroup : public vector<Box*>
    {
    public:

        //constructor
        NodeGroup(Grid * g=NULL);

        int add(Box * box);

        int remove(Box * box);

        /*!
            Return an array of size of the number of active leaf nodes
            and include in the array the nodeid.

            nodeid_array[nodenumber] = nodeid
        */
        int * get_nodeid_array();

        /*!
         * return a nodelay_array, which contains the number of
            nodes per layer:
            nodelay_array[layer] = nodes in layer
            This method numbers the nodes by layer, starting with the
            lowest layer number.
         *
         */
        int * get_nodelay_array();


        ///make all nodes in the given layer as diabled
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
        Grid * grid;          //the hosting Grid
        int * nodeid_array;   //node ids (mapping node number to node id)
        int * nodelay_array;  //number of nodes per layer
        int nlay;             //total number of layers

        bool rebuild_nodeid_array; //rebuid node id

        int nextnodenumber;


        int chunksize; ///what is a chunk size??
        list<int> removed_nodes;

    protected:

        int total_leaf_number; //number of leaves store in this group (totoal number of leaves)
        int active_leaf_number; //number of leaves store in this group (active nodes only)

        /*!
         * build information of nodelay_array and nodeid_array
         *
         */
        void build_node_info();
    };

} // end namespace cusg

