
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
    ModflowGrid(const ModflowGridBlock & rawdata);

    ModflowGrid(int lay, int row, int col, double delr, double delc, vector<double *> & bot,
                double xoffset, double yoffset, double angle);

    ~ModflowGrid() override {}

    /// Return a numpy array of size (nodes) that contains the area for all cells.
    vector<double> get_cell_areas() override;

    ///Return the box_index for the modflow grid.
    virtual int get_box_index(int k, int i, int j);

    ///Return the 3D indices for the specified box_index.
    Index3d get_indices(int box_index);

    ///get box from id
    //a MODFLOW grid numbers its cells as its boxes
    Box * get_cell_box(int cell) override;

    //the box in layer k, row i, column j
    virtual Box * get_cell_box(int k, int i, int j);

    ///get box connections from id
    ///Return a sorted list of connections, where a connection is
    ///a list [to_box_index, fldir].  fldir is -1, +1, -2, +2, -3, +3 for -x, +x, -y, +y, -z, +z.
    list<pair<int, int>> get_node_connections(int cell) override;

    /*! \brief Obtain CSRData

        Return a UsgCsrData object, which contains all of the necessary
        information for creating an unstructured grid model.

        not implemented yet
     */
    CSRData get_usg_csr_data() override;

    ///number the nodes (i.e. leaves) in the grid
    void number_nodes() override;

    ////rotate all the boxes
    void get_rotation(double & cx, double & cy, double & angle)
    {
        cx = this->x_centers[0] - 0.5 * delr[0];
        cy = this->y_centers[nrow - 1] - 0.5 * delc[nrow - 1];
        angle = rotation;
    }

protected:
    void init(int nlay, int nrow, int ncol, const vector<double *> & botm);

    ///create nodes and add them to node_group
    void create_nodes();

    ///add connection between boxes
    ///fldir is -1, +1, -2, +2, -3, +3 for -x, +x, -y, +y, -z, +z
    void add_connections_3d();
};

} //namespace cusg
