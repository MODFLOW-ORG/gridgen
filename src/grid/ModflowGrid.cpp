
/*! \file ModflowGrid.cpp

\brief Modflow grid in 3D

\author <a href="http://masc.cs.gmu.edu">MASC group</a>, George Mason University 
\author <a href="http://profile.usgs.gov/langevin/">Christian Langevin</a>, USGS
\bug    No known bugs.
*/

#include "ModflowGrid.h"

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

ModflowGrid::ModflowGrid(const modflow_grid_raw_data & rawdata)
    : ModflowGrid2D(rawdata.nrow, rawdata.ncol, rawdata.delr, rawdata.delc, rawdata.x_offset,
                    rawdata.y_offset, rawdata.rotation_angle)
{
    init(rawdata.nlay, rawdata.nrow, rawdata.ncol, rawdata.bot);
}

ModflowGrid::ModflowGrid(int lay, int row, int col, double dr, double dc, vector<double *> & bot,
                         double xoffset, double yoffset, double angle)
    : ModflowGrid2D(row, col, dr, dc, xoffset, yoffset, angle)
{
    init(lay, row, col, bot);
}

/// Return a numpy array of size (nodes) that contains the area for all cells.
double * ModflowGrid::get_cell_areas()
{
    double * area = new double[nodes];
    assert(area);

    for (int k = 0; k < nlay; k++)
    {
        for (int i = 0; i < nrow; i++)
        {
            for (int j = 0; j < ncol; j++)
            {
                int nodeid = get_nodeid(k, i, j);
                area[nodeid] = delr[j] * delc[i];
            } //end j
        } //end i
    } //end k

    return area;
}

///Return the nodeid for the modflow grid.
int ModflowGrid::get_nodeid(int k, int i, int j)
{
    return k * nrow * ncol + i * ncol + j;
}

int ModflowGrid::getFirstNLayer()
{
    return 0;
}

///Return the 3D indices for the specified nodeid.
Index3d ModflowGrid::get_indices(int nodeid)
{

    if (nodeid > (nlay * nrow * ncol))
    {
        cerr << "Nodeid = " << nodeid << " is not valid..." << endl;
        assert(false);
    }

    int k = int(nodeid / nrow / ncol);
    int i = int((nodeid - k * nrow * ncol) / ncol);
    int j = nodeid - k * nrow * ncol - i * ncol;

    return Index3d(k, i, j);
}

///get box from id
Box * ModflowGrid::get_nodeobj(int nodeid)
{
    return nodegroup[nodeid];
}

Box * ModflowGrid::get_nodeobj(int i, int j, int k)
{
    int id = get_nodeid(i, j, k);
    return get_nodeobj(id);
}

///get box connections from id
///Return a sorted list of connections, where a connection is
///a list [tonodeid, fldir].  fldir is -1, +1, -2, +2, -3, +3 for -x, +x, -y, +y, -z, +z.
list<pair<int, int>> ModflowGrid::get_node_connections(int nodeid)
{
    list<pair<int, int>> connections;
    Index3d kij = get_indices(nodeid);
    int k = kij[0];
    int i = kij[1];
    int j = kij[2];

    //up
    if (k > 0)
    {
        int nn = get_nodeid(k - 1, i, j);
        connections.push_back(make_pair(nn, 3));
    }

    //back
    if (i > 0)
    {
        int nn = get_nodeid(k, i - 1, j);
        connections.push_back(make_pair(nn, 2));
    }
    //left
    if (j > 0)
    {
        int nn = get_nodeid(k, i, j - 1);
        connections.push_back(make_pair(nn, -1));
    }
    //right
    if (j < ncol - 1)
    {
        int nn = get_nodeid(k, i, j + 1);
        connections.push_back(make_pair(nn, 1));
    }
    //front
    if (i < nrow - 1)
    {
        int nn = get_nodeid(k, i + 1, j);
        connections.push_back(make_pair(nn, -2));
    }
    //down
    if (k < nlay - 1)
    {
        int nn = get_nodeid(k + 1, i, j);
        connections.push_back(make_pair(nn, -3));
    }

    return connections;
}

/*! \brief Obtain CSRData

        Return a UsgCsrData object, which contains all of the necessary
        information for creating an unstructured grid model.

        TODO: We also need to check if the neighboring cells are active or not
     */
CSRData ModflowGrid::get_usg_csr_data()
{
    cout << "- Creating the UsgCsrData object." << endl;

    int box_size = nodegroup.size();

    //go through each cell and count diagonal position and connections
    int nja = (2 * (ncol - 1) * nrow * nlay + //x connections
               2 * ncol * (nrow - 1) * nlay + //y connections
               2 * ncol * nrow * (nlay - 1) + //z connections
               ncol * nrow * nlay);           //diagonal

    int * ia = new int[box_size + 1];
    int * ja = new int[nja];
    double * area = get_cell_areas();
    double * cl1 = new double[nja];
    double * cl2 = new double[nja];
    double * fahl = new double[nja];
    int * fldr = new int[nja];
    int * iac = new int[box_size]; // number of connection +1 for each cell, array of size  nodes

    assert(ia && ja && area && cl1 && cl2 && fahl && fldr && iac);

    //Go through each connection and fill the arrays.
    int japos = 0;
    int nodenumber = 0;
    double dz = 1.0;

    for (int k = 0; k < nlay; k++)
    {
        for (int i = 0; i < nrow; i++)
        {
            for (int j = 0; j < ncol; j++)
            {
                //diagonal position
                ia[nodenumber] = japos + 1;
                ja[japos] = nodenumber + 1;
                cl1[japos] = 0;
                cl2[japos] = 0;
                fahl[japos] = 0;
                fldr[japos] = 0;
                iac[nodenumber] = 1;
                japos += 1;

                double dzkij = 0;

                if (is3d)
                {
                    dzkij = botm(k, i, j) - botm(k + 1, i, j);
                }

                if (k > 0) //up
                {
                    int tonodenumber = get_nodeid(k - 1, i, j);
                    ja[japos] = tonodenumber + 1;
                    fldr[japos] = +3;
                    cl1[japos] = 0.5 * (botm(k, i, j) - botm(k + 1, i, j));
                    cl2[japos] = 0.5 * (botm(k - 1, i, j) - botm(k, i, j));
                    fahl[japos] = delr[j] * delc[i];
                    iac[nodenumber]++;
                    japos += 1;
                }

                if (i > 0) //back
                {
                    int tonodenumber = get_nodeid(k, i - 1, j);
                    ja[japos] = tonodenumber + 1;
                    fldr[japos] = +2;
                    cl1[japos] = 0.5 * (delc[i]);
                    cl2[japos] = 0.5 * (delc[i - 1]);
                    if (is3d) dz = 0.5 * (botm(k, i - 1, j) - botm(k + 1, i - 1, j) + dzkij);

                    fahl[japos] = delr[j] * dz;
                    iac[nodenumber]++;
                    japos += 1;
                }

                if (j > 0) //left
                {
                    int tonodenumber = get_nodeid(k, i, j - 1);
                    ja[japos] = tonodenumber + 1;
                    fldr[japos] = -1;
                    cl1[japos] = 0.5 * (delr[j]);
                    cl2[japos] = 0.5 * (delr[j - 1]);
                    if (is3d) dz = 0.5 * (botm(k, i, j - 1) - botm(k + 1, i, j - 1) + dzkij);

                    fahl[japos] = delc[i] * dz;
                    iac[nodenumber]++;
                    japos += 1;
                }

                if (j < ncol - 1) //right
                {
                    int tonodenumber = get_nodeid(k, i, j + 1);
                    ja[japos] = tonodenumber + 1;
                    fldr[japos] = +1;
                    cl1[japos] = 0.5 * (delr[j]);
                    cl2[japos] = 0.5 * (delr[j + 1]);
                    if (is3d) dz = 0.5 * (botm(k, i, j + 1) - botm(k + 1, i, j + 1) + dzkij);

                    fahl[japos] = delc[i] * dz;
                    iac[nodenumber]++;
                    japos += 1;
                }

                if (i < nrow - 1) //front
                {
                    int tonodenumber = get_nodeid(k, i + 1, j);
                    ja[japos] = tonodenumber + 1;
                    fldr[japos] = -2;
                    cl1[japos] = 0.5 * (delc[i]);
                    cl2[japos] = 0.5 * (delc[i + 1]);
                    if (is3d) dz = 0.5 * (botm(k, i + 1, j) - botm(k + 1, i + 1, j) + dzkij);

                    fahl[japos] = delr[j] * dz;
                    iac[nodenumber]++;
                    japos += 1;
                }

                if (k < nlay - 1) //down
                {
                    int tonodenumber = get_nodeid(k + 1, i, j);
                    ja[japos] = tonodenumber + 1;
                    fldr[japos] = -3;
                    cl1[japos] = 0.5 * (botm(k, i, j) - botm(k + 1, i, j));
                    cl2[japos] = 0.5 * (botm(k + 1, i, j) - botm(k + 2, i, j));
                    fahl[japos] = delr[j] * delc[i];
                    iac[nodenumber]++;
                    japos += 1;
                }

                nodenumber += 1;

            } //end for j
        } //end for i
    } //end for k

    ia[nodenumber] = japos + 1;

    return CSRData(nodes, nja, area, ia, ja, cl1, cl2, fahl, fldr, iac);
}

///number the nodes (i.e. leaves) in the grid
void ModflowGrid::number_nodes()
{
    unsigned int number = 0;
    for (int l = 0; l < nlay; l++)
    {
        for (int r = 0; r < nrow; r++)
        {
            for (int c = 0; c < ncol; c++)
            {
                int bid = get_nodeid(l, r, c);
                Box * box = nodegroup[bid];
                if (box->active) box->number = number++;
            } //end c
        } //end r
    } //end l
}

void ModflowGrid::init(int nlay, int nrow, int ncol, const vector<double *> & botm)
{
    //log the creation
    cout << "- A ModflowGrid object was instantiated." << endl;

    //assign the arguments to the object.
    this->nlay = nlay;
    this->nodes = nlay * nrow * ncol;
    if (nlay > 1 || botm.size() > 1) is3d = true;

    //copy bottom layers
    if (botm.empty() == false)
    {
        int csize = nrow * ncol;
        this->bot = new double[this->nodes + csize];
        this->top = new double[csize];
        ;
        assert(bot && top);

        memcpy(top, botm[0], csize * sizeof(double));

        for (int l = 0; l < nlay + 1; l++)
        {
            int id = l * csize;
            memcpy(&bot[id], botm[l], csize * sizeof(double));
        } //end for
    } //

    create_nodes();
}

///create nodes and add them to nodegroup
void ModflowGrid::create_nodes()
{
    // Create the nodegroup.
    nodegroup = NodeGroup(this);
    nodegroup.nlay = nlay;

    for (int k = 0; k < nlay; k++)
    {
        for (int i = 0; i < nrow; i++)
        {
            double y = Y[i];

            for (int j = 0; j < ncol; j++)
            {
                double x = X[j];

                Point3d position;
                Vector3d dxdydz;

                if (bot != NULL)
                {
                    double z = 0.5 * (botm(k, i, j) + botm(k + 1, i, j));
                    double dz = (botm(k, i, j) - botm(k + 1, i, j));
                    position = Point3d(x, y, z);
                    dxdydz = Vector3d(delr[j], delc[i], dz);
                }
                else
                {
                    position.set(x, y, 0);
                    dxdydz.set(delr[j], delc[i], 0);
                }

                Box * nodeobj = new Box(position, dxdydz);
                assert(nodeobj);
                nodeobj->layer = k;

                //TODO: add nodeobj to something?
                nodegroup.add(nodeobj);

            } //end i
        } //end j
    } //end k

    //Create the connections for the root nodes.
    if (is3d)
        add_connections_3d();
    else
        add_connections();

    nodegroup.get_nodeid_array(); //this forces to build node number/id
}

///add connection between boxes
///fldir is -1, +1, -2, +2, -3, +3 for -x, +x, -y, +y, -z, +z.
void ModflowGrid::add_connections_3d()
{
    int last_lay = nlay - 1;
    int last_row = nrow - 1;
    int last_col = ncol - 1;

    for (int k = 0; k < nlay; k++)
    {
        for (int i = 0; i < nrow; i++)
        {
            for (int j = 0; j < ncol; j++)
            {
                int nodeid = get_nodeid(k, i, j);

                Box * nodeobj = nodegroup[nodeid];

                if (i < last_row)
                {
                    //make connection to front
                    int tonode = get_nodeid(k, i + 1, j);
                    int fldir = -2;
                    nodeobj->addConnection(nodegroup[tonode], fldir);
                }

                if (i > 0)
                {
                    //make connection to back
                    int tonode = get_nodeid(k, i - 1, j);
                    int fldir = 2;
                    nodeobj->addConnection(nodegroup[tonode], fldir);
                }

                if (j < last_col)
                {
                    //make connection to right
                    int tonode = get_nodeid(k, i, j + 1);
                    int fldir = 1;
                    nodeobj->addConnection(nodegroup[tonode], fldir);
                }

                if (j > 0)
                {
                    //make connection to left
                    int tonode = get_nodeid(k, i, j - 1);
                    int fldir = -1;
                    nodeobj->addConnection(nodegroup[tonode], fldir);
                }

                if (k < last_lay)
                {
                    //make connection down
                    int tonode = get_nodeid(k + 1, i, j);
                    int fldir = -3;
                    nodeobj->addConnection(nodegroup[tonode], fldir);
                }

                if (k > 0)
                {
                    //make connection up
                    int tonode = get_nodeid(k - 1, i, j);
                    int fldir = 3;
                    nodeobj->addConnection(nodegroup[tonode], fldir);
                }

            } //end j
        } //end i
    } //end k
}

} //namespace cusg
