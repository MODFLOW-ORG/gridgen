
/*! \file ModflowGrid.h

\brief Modflow grid in 3D

\author <a href="http://masc.cs.gmu.edu">MASC group</a>, George Mason University 
\author <a href="http://profile.usgs.gov/langevin/">Christian Langevin</a>, USGS
\bug    No known bugs.
*/

#pragma once

#include "ModflowGrid2D.h"
#include "def_parser.h"
#include <limits.h>

namespace cusg
{

/*! \brief MODFLOW grid 3D (nlay, nrow, ncol)

    ModflowGrid is a MODFLOW grid class.  This is a structured grid of shape 
    (nlay, nrow, ncol).  The structure assumes that the first layer is the top 
    layer, and that the first row corresponds to the north side of the grid, 
    so that index [1, 1, 1] corresponds to the top layer, and the northwest 
    corner of the grid.
    
    ModflowGrid is a subclass of ModflowGrid2D.
    
    ModflowGrid is also a subclass of the Grid, which contains 
    generic methods for all grid classes.
    
*/

class ModflowGrid : public ModflowGrid2D
{
public:
    ModflowGrid(const modflow_grid_raw_data & rawdata);

    ModflowGrid(int lay, int row, int col, double dr, double dc, vector<double *> & bot,
                double xoffset, double yoffset, double angle);

    virtual ~ModflowGrid() {}

    /// Return a numpy array of size (nodes) that contains the area for all cells.
    double * get_cell_areas();

    ///Return the nodeid for the modflow grid.
    virtual int get_nodeid(int k, int i, int j);

    virtual int getFirstNLayer();

    ///Return the 3D indices for the specified nodeid.
    Index3d get_indices(int nodeid);

    ///get box from id
    Box * get_nodeobj(int nodeid);

    virtual Box * get_nodeobj(int i, int j, int k);

    ///get box connections from id
    ///Return a sorted list of connections, where a connection is
    ///a list [tonodeid, fldir].  fldir is -1, +1, -2, +2, -3, +3 for -x, +x, -y, +y, -z, +z.
    list<pair<int, int>> get_node_connections(int nodeid);

    /*! \brief Obtain CSRData

        Return a UsgCsrData object, which contains all of the necessary
        information for creating an unstructured grid model.

        not implemented yet
     */
    CSRData get_usg_csr_data();

    ///number the nodes (i.e. leaves) in the grid
    virtual void number_nodes();

    ////rotate all the boxes
    void getRotatePara(double & cx, double & cy, double & angle)
    {
        cx = this->X[0] - 0.5 * delr[0];
        cy = this->Y[nrow - 1] - 0.5 * delc[nrow - 1];
        //cy = this->Y[nrow - 1]  + 0.5 * delc[nrow - 1];//- 0.5 * delc[nrow - 1];
        angle = rotation;
    }

protected:
    void init(int nlay, int nrow, int ncol, const vector<double *> & botm);

    ///create nodes and add them to nodegroup
    void create_nodes();

    ///add connection between boxes
    ///fldir is -1, +1, -2, +2, -3, +3 for -x, +x, -y, +y, -z, +z
    void add_connections_3d();
};

} //namespace cusg
