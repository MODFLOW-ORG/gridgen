
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

	
	/*! \brief Instantiate a ModflowGrid object.
        
        This method is used to instantiate a ModflowGrid object.  The 
        following arguments are required:
           \param nlay: number of layers
           \param nrow: number of rows
           \param ncol: number of columns
           \param delr: a ModflowArray object of size [ncol]
           \param delc: a ModflowArray object of size [nrow]
           \param botm: a ModflowArray object of size [nlay + 1, nrow, ncol]
            
        Optional arguments that can be specified using "...":

           \param xoffset: the offset of the grid in the x direction
           \param yoffset: the offset of the grid in the y direction
           \param rotation: the grid rotation angle in degrees relative to the lower  
                      left corner.
                      
        The local coordinate system has a (0, 0) x,y location that corresponds 
        to the lower left corner of the grid.
        
        \note This is a variadic function.
    */
	ModflowGrid::ModflowGrid(int nlay, int nrow, int ncol, double * delr, double * delc, vector< vector< vector<int> > >& botm,  float xoffset, float yoffset, float rotation)
	:ModflowGrid2D(nrow,ncol,delr,delc,xoffset,yoffset,rotation)
	{
		init(nlay, nrow, ncol, botm);
	}
   
	ModflowGrid::ModflowGrid(int nlay, int nrow, int ncol, double delr, double delc, 
		vector< vector< vector<int> > >& botm,  float xoffset, float yoffset, float rotation)
	:ModflowGrid2D(nrow,ncol,delr,delc,xoffset,yoffset,rotation)
	{
		init(nlay, nrow, ncol, botm);
	}
 //  
    ModflowGrid::ModflowGrid(const ModflowGrid &other) : ModflowGrid2D(other.nrow,
        other.ncol, other.delr[0], other.delc[0], other.xoffset,
        other.yoffset, other.rotation)
    {
        bot=other.bot;
        init(other.nlay, other.nrow, other.ncol, other.bot);
    }

    ModflowGrid::ModflowGrid(const modflow_grid_raw_data &rawdata) :
        ModflowGrid2D(rawdata.nrow, rawdata.ncol, rawdata.delr, rawdata.delc, rawdata.x_offset, rawdata.y_offset, rawdata.rotation_angle)
    {
        init(rawdata.nlay, rawdata.nrow, rawdata.ncol, rawdata.bot);
    }

	ModflowGrid::ModflowGrid(int lay,int row,int col,double dr,double dc, vector<double*>& bot,double xoffset,double yoffset,double angle)
		:ModflowGrid2D(row, col, dr, dc, xoffset, yoffset, angle)
	{
		init(lay, row, col, bot);
	}
    
    ModflowGrid::ModflowGrid() { /* be sure to call init() once you set up the grid! */
		
	};
    
	/// Return a numpy array of size (nodes) that contains the area for all cells.
	double *  ModflowGrid::get_cell_areas()
	{   
        double * area=new double[nodes];
        assert(area);    
        
        for( int  k=0; k< nlay; k++)
        {
            for ( int i=0; i<nrow; i++)
            {
                for( int j=0; j< ncol; j++)
                {
                    int nodeid = get_nodeid(k, i, j);
                    area[nodeid] = delr[j] * delc[i];
                }//end j
            }//end i
        }//end k
        
        return area;
	}
        
	/*!
        \brief Return the four vertices for the specified nodeid.
        
        \verbatim
        3-------2   7-------6
        |  top  |   |  bot  |
        |       |   |       |
        0-------1   4-------5
        \endverbatim
    */    
	vector<Point3d> ModflowGrid::get_vertices(int nodeid)
	{
        Index3d kij = get_indices(nodeid);
        int k=kij[0];
        int i=kij[1];
        int j=kij[2];
        
        vector<Point3d> vertices(8, Point3d() );
        
        double top=botm(k,i,j);
        vertices[0].set(Xe[j], Ye[i + 1], top);
        vertices[1].set(Xe[j + 1], Ye[i + 1], top);
        vertices[2].set(Xe[j + 1], Ye[i], top);
        vertices[3].set(Xe[j], Ye[i], top);

        double bot=botm(k+1,i,j);
        vertices[4].set(Xe[j], Ye[i + 1], bot);
        vertices[5].set(Xe[j + 1], Ye[i + 1], bot);
        vertices[6].set(Xe[j + 1], Ye[i], bot);
        vertices[7].set(Xe[j], Ye[i], bot);

        return vertices;
	}
	
	int ModflowGrid::get_nodeid(int i, int j){ return ModflowGrid2D::get_nodeid(i,j);}

	///Return the nodeid for the modflow grid.
	int ModflowGrid::get_nodeid(int k, int i, int j)
	{
        return k * nrow * ncol + i * ncol + j;
    }

	int ModflowGrid::getFirstNLayer(){return 0;}
    
    ///Return the 3D indices for the specified nodeid.
    Index3d ModflowGrid::get_indices(int nodeid)
    {
		
        if (nodeid >(nlay * nrow * ncol))
        { 
            cerr<<"Nodeid = "<<nodeid<<" is not valid..."<<endl;
            assert(false);
        }
            
        int k = int( nodeid / nrow / ncol );
        int i = int( (nodeid - k * nrow * ncol) / ncol );
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
        int id=get_nodeid(i,j,k);
        return get_nodeobj(id);
    }
    
    ///get box connections from id
    ///Return a sorted list of connections, where a connection is
    ///a list [tonodeid, fldir].  fldir is -1, +1, -2, +2, -3, +3 for -x, +x, -y, +y, -z, +z.
    list< pair<int,int> >  ModflowGrid::get_node_connections(int nodeid)
    {
    	list< pair<int,int> > connections;
    	Index3d ijk=get_indices(nodeid);
    	int i=ijk[0];
    	int j=ijk[1];
    	int k=ijk[2];
    	    	
    	//up
    	if( k > 0 )
    	{
            int nn = get_nodeid(k - 1, i, j);
            connections.push_back(make_pair(nn, 3));
        }   
        
    	//back
    	if(i>0)
    	{
            int nn = get_nodeid(k, i-1, j);
            connections.push_back(make_pair(nn, 2));
    	}
 		//left
        if(j > 0)
        {
            int nn = get_nodeid(k, i, j - 1);
            connections.push_back(make_pair(nn, -1));
        }
    	//right
        if(j < ncol - 1)
        {
            int nn = get_nodeid(k, i, j + 1);
            connections.push_back(make_pair(nn, 1));
    	}
        //front
        if( i < nrow - 1 )
        {
            int  nn = get_nodeid(k, i + 1, j);
            connections.push_back( make_pair(nn, -2) );
        }   
        //down
        if( k < nlay - 1 )
        {
            int nn = get_nodeid(k + 1, i, j);
            connections.push_back( make_pair(nn, -3) );
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
        cout<<"- Creating the UsgCsrData object."<<endl;

        int box_size=nodegroup.size();

        //go through each cell and count diagonal position and connections
        int nja = ( 2 * (ncol - 1) * nrow * nlay +     //x connections
                    2 * ncol * (nrow - 1) * nlay +     //y connections
                    2 * ncol * nrow * (nlay - 1) +     //z connections
                    ncol * nrow * nlay );              //diagonal


        int * ia = new int[box_size+1];
        int * ja = new int[nja];
        double * area = get_cell_areas();
        float * cl1 = new float[nja];
        float * cl2 = new float[nja];
        float * fahl = new float[nja];
        int * fldr = new int[nja];
        int * iac = new int[box_size]; // number of connection +1 for each cell, array of size  nodes

        assert(ia && ja && area && cl1 && cl2 && fahl && fldr && iac);

        //Go through each connection and fill the arrays.
        int japos = 0;
        int nodenumber = 0;
        int dz = 1.0;

        for( int  k=0; k< nlay; k++)
        {
            for ( int i=0; i<nrow; i++)
            {
                for( int j=0; j< ncol; j++)
                {
                    //diagonal position
                    ia[nodenumber] = japos+1;
                    ja[japos] = nodenumber+1;
                    cl1[japos] = 0;
                    cl2[japos] = 0;
                    fahl[japos] = 0;
                    fldr[japos] = 0;
                    iac[nodenumber] = 1;
                    japos += 1;

                    int dzkij=0;

                    if(is3d)
                    {
                        dzkij = botm(k,i,j) - botm(k+1,i,j);
                    }

                    if( k > 0 ) //up
                    {
                        int tonodenumber = get_nodeid(k - 1, i, j);
                        ja[japos] = tonodenumber+1;
                        fldr[japos] = -3;
                        cl1[japos] = 0.5 * (botm(k,i,j) - botm(k+1,i,j));
                        cl2[japos] = 0.5 * (botm(k-1,i,j) - botm(k,i,j));
                        fahl[japos] = delr[j] * delc[i];
                        iac[nodenumber]++;
                        japos += 1;
                    }

                    if( i > 0 )//back
                    {
                        int tonodenumber = get_nodeid(k, i - 1, j);
                        ja[japos] = tonodenumber+1;
                        fldr[japos] = +2;
                        cl1[japos] = 0.5 * (delc[i]);
                        cl2[japos] = 0.5 * (delc[i - 1]);
                        if(is3d)
                            dz = 0.5 * (botm(k,i-1,j) -  botm(k+1,i-1,j) + dzkij);

                        fahl[japos] = delr[j] * dz;
                        iac[nodenumber]++;
                        japos += 1;
                    }

                    if( j > 0 ) //left
                    {
                        int tonodenumber = get_nodeid(k, i, j - 1);
                        ja[japos] = tonodenumber+1;
                        fldr[japos] = -1;
                        cl1[japos] = 0.5 * (delr[j]);
                        cl2[japos] = 0.5 * (delr[j - 1]);
                        if( is3d )
                            dz = 0.5 * (botm(k,i,j - 1) - botm(k + 1,i,j - 1) + dzkij);

                        fahl[japos] = delc[i] * dz;
                        iac[nodenumber]++;
                        japos += 1;
                    }

                    if( j < ncol - 1 ) //right
                    {
                        int tonodenumber = get_nodeid(k, i, j + 1);
                        ja[japos] = tonodenumber+1;
                        fldr[japos] = +1;
                        cl1[japos] = 0.5 * (delr[j]);
                        cl2[japos] = 0.5 * (delr[j + 1]);
                        if(is3d)
                            dz = 0.5 * (botm(k,i,j - 1) - botm(k + 1,i,j - 1) + dzkij);

                        fahl[japos] = delc[i] * dz;
                        iac[nodenumber]++;
                        japos += 1;
                    }

                    if( i < nrow - 1 ) //front
                    {
                        int tonodenumber = get_nodeid(k, i + 1, j);
                        ja[japos] = tonodenumber+1;
                        fldr[japos] = -2;
                        cl1[japos] = 0.5 * (delc[i]);
                        cl2[japos] = 0.5 * (delc[i + 1]);
                        if(is3d)
                            dz = 0.5 * (botm(k,i,j - 1) - botm(k + 1,i,j - 1) + dzkij);

                        fahl[japos] = delr[j] * dz;
                        iac[nodenumber]++;
                        japos += 1;
                    }

                    if( k < nlay - 1 ) //down
                    {
                        int tonodenumber = get_nodeid(k + 1, i, j);
                        ja[japos] = tonodenumber+1;
                        fldr[japos] = -3;
                        cl1[japos] = 0.5 * ( botm(k,i,j) -  botm(k + 1,i,j) );
                        cl2[japos] = 0.5 * ( botm(k + 1,i,j) - botm(k + 2,i,j));
                        fahl[japos] = delr[j] * delc[i];
                        iac[nodenumber]++;
                        japos += 1;
                    }

                    nodenumber += 1;

                }  //end for j
            }//end for i
        }//end for k

        ia[nodenumber] = japos;

        return CSRData(nodes, nja, area, ia, ja, cl1, cl2, fahl, fldr, iac);
    }

	///this is a temporary SWIG/Python interface. Will find a better way to handle this
	void ModflowGrid::get_usg_csr_data(void * data)
	{
		CSRData * csr=(CSRData*)data;
		*csr=get_usg_csr_data();
	}

	///number the nodes (i.e. leaves) in the grid
    void ModflowGrid::number_nodes()
    {
        unsigned int number=0;
        for(int l=0;l<nlay;l++)
        {
            for(int r=0;r<nrow;r++)
              {
                  for(int c=0;c<ncol;c++)
                  {
                      int bid=get_nodeid(l,r,c);
                      Box * box=nodegroup[bid];
                      if(box->active)
                      box->number=number++;
                  }//end c
            }//end r
        }//end l
    }

	void ModflowGrid::updateNeighbors()
	{
		for(int l=0;l<nlay;l++)
		{
			for(int r=0;r<nrow;r++)
			{
				for(int c=0;c<ncol;c++)
				{
					int bid=get_nodeid(l,r,c);
					Box * box=nodegroup[bid];
					//pDown and pUp
					if(l == 0)
						box->pDown = NULL;
					else
						box->pDown = nodegroup[get_nodeid(l-1,r,c)];
					if(l == nlay-1)
						box->pUp = NULL;
					else
						box->pUp = nodegroup[get_nodeid(l+1, r, c)];

		
					if(r == 0)
						box->pChildren[0] = NULL;
					else
						box->pChildren[0] = nodegroup[get_nodeid(l,r-1,c)];
					if(r == nrow-1)
						box->pChildren[2] = NULL;
					else
						box->pChildren[2] = nodegroup[get_nodeid(l, r+1, c)];

					if(c == 0)
						box->pChildren[3] = NULL;
					else
						box->pChildren[3] = nodegroup[get_nodeid(l,r,c-1)];

					if(c == ncol - 1)
						box->pChildren[1] = NULL;
					else
						box->pChildren[1] = nodegroup[get_nodeid(l, r, c+1)];
				}
			}
		}
	}

    void ModflowGrid::init(int nlay, int nrow, int ncol, double * botm)
    {
        //log the creation
        cout<<"- A ModflowGrid object was instantiated."<<endl;

        //assign the arguments to the object.
        this->nlay = nlay;
        this->nodes = nlay * nrow * ncol;

        //copy bottom layers
        if(botm!=NULL)
        {	
        	long botsize= (nlay+1) * nrow * ncol;
            this->bot=new double[botsize];
            this->top=bot;
            assert(bot);
            memcpy(bot,botm,botsize*sizeof(double));
        }

        create_nodes();
    }

	void ModflowGrid::init(int nlay, int nrow, int ncol, const vector< vector< vector<int> > > &botm)
	{
		//log the creation
        cout<<"- A ModflowGrid object was instantiated."<<endl;

        //assign the arguments to the object.
        this->nlay = nlay;
        this->nodes = nlay * nrow * ncol;
        if(nlay>1 || botm.size()>1) is3d=true;

        //copy bottom layers
        if(botm.empty()==false)
        {
	        long botsize= (nlay+1) * nrow * ncol;
            this->bot=new double[botsize];
            this->top=bot;
            assert(bot);

            for(int l=0, id=0;l<nlay+1;l++) for(int r=0;r<nrow;r++) for(int c=0;c<ncol;c++) bot[id++]=botm[l][r][c];
        }

        create_nodes();
	}

    void ModflowGrid::init(int nlay, int nrow, int ncol, const vector<double*> &botm)
    {
        //log the creation
        cout<<"- A ModflowGrid object was instantiated."<<endl;

        //assign the arguments to the object.
        this->nlay = nlay;
        this->nodes = nlay * nrow * ncol;
        if(nlay>1 || botm.size()>1 ) is3d=true;

        //copy bottom layers
        if(botm.empty()==false)
        {
            int csize=nrow * ncol;
            this->bot=new double[this->nodes+csize];
            this->top=new double[csize];;
            assert(bot && top);

            memcpy(top,botm[0],csize*sizeof(double));

            for(int l=0;l<nlay+1;l++)
            {
                int id=l*csize;
                memcpy(&bot[id],botm[l],csize*sizeof(double));
            }//end for
        }//

        create_nodes();
    }

    ///create nodes and add them to nodegroup
    void ModflowGrid::create_nodes()
    {
        // Create the nodegroup.
        nodegroup = NodeGroup(this);
        nodegroup.nlay=nlay;

        for( int k=0; k< nlay; k++)
        {
            for( int i=0; i< nrow; i++)
            {
                double y = Y[i];

                for( int j=0; j< ncol; j++)
                {
                    double x = X[j];

                    Point3d position;
                    Vector3d dxdydz;

                    if(bot!=NULL)
                    {
                        double z = 0.5 * (botm(k,i,j)+botm(k+1,i,j));
                        double dz = (botm(k,i,j) - botm(k+1,i,j));
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
                    nodeobj->layer=k;

                    //TODO: add nodeobj to something?
                    nodegroup.add(nodeobj);

                }//end i
            }//end j
        }//end k

        //Create the connections for the root nodes.
        if(is3d)
            add_connections_3d();
        else
            add_connections();

        nodegroup.get_nodeid_array(); //this forces to build node number/id
    }

    ///add connection between boxes
    ///fldir is -1, +1, -2, +2, -3, +3 for -x, +x, -y, +y, -z, +z.
    void ModflowGrid::add_connections_3d()
    {
        int last_lay=nlay - 1;
        int last_row=nrow - 1;
        int last_col=ncol - 1;

        for( int k=0;k<nlay;k++)
        {
            for( int i=0;i<nrow;i++)
            {
                for( int j=0;j<ncol;j++)
                {
                    int nodeid = get_nodeid(k, i, j);

                    //cout<<"Adding connections for nodeid: "<<nodeid<<endl;
                    Box * nodeobj = nodegroup[nodeid];

                    if( i < last_row )
                    {
                        //make connection to front
                        int tonode = get_nodeid(k, i + 1, j);
                        int fldir = -2;
                        nodeobj->addConnection(nodegroup[tonode], fldir);
                    }

                    if( i > 0 )
                    {
                        //make connection to back
                        int tonode = get_nodeid(k, i - 1, j);
                        int fldir = 2;
                        nodeobj->addConnection(nodegroup[tonode], fldir);
                    }

                    if( j < last_col )
                    {
                        //make connection to right
                        int tonode = get_nodeid(k, i, j + 1);
                        int fldir = 1;
                        nodeobj->addConnection(nodegroup[tonode], fldir);
                    }

                    if( j > 0 )
                    {
                        //make connection to left
                        int tonode = get_nodeid(k, i, j - 1);
                        int fldir = -1;
                        nodeobj->addConnection(nodegroup[tonode], fldir);
                    }

                    if( k < last_lay )
                    {
                        //make connection down
                        int tonode = get_nodeid(k + 1, i, j);
                        int fldir = -3;
                        nodeobj->addConnection(nodegroup[tonode], fldir);
                    }

                    if( k > 0 )
                    {
                        //make connection up
                        int tonode = get_nodeid(k - 1, i, j);
                        int fldir = 3;
                        nodeobj->addConnection(nodegroup[tonode], fldir);
                    }

                }//end j
            }//end i
        }//end k
    }
	

	//void ModflowGrid::rotate()
	//{
	//	//if(rotation==0)
	//	//	return;
	//	double cx = this->X[0] - 0.5 * delr[0]; 
	//	double cy = this->Y[nrow - 1] - 0.5 * delc[nrow - 1];
	//	NodeGroup::iterator nit = nodegroup.begin();
	//	for(; nit != nodegroup.end(); ++nit)
	//	{
	//		Box* b = *nit;
	//		b->ldX = cx + (b->x - b->dx / 2  - cx) * cos(rotation) - (b->y - b->dy / 2 - cy) * sin(rotation);
	//		b->ldY = cy + (b->y - b->dy / 2 - cy) * cos(rotation) + (b->x - b->dx / 2 - cx) * sin(rotation);
	//		
	//		b->luX = cx + (b->x - b->dx / 2  - cx) * cos(rotation) - (b->y + b->dy / 2 - cy) * sin(rotation);
	//		b->luY = cy + (b->y + b->dy / 2 - cy) * cos(rotation) + (b->x - b->dx / 2 - cx) * sin(rotation);
	//		
	//		b->rdX = cx + (b->x + b->dx / 2  - cx) * cos(rotation) - (b->y - b->dy / 2 - cy) * sin(rotation);
	//		b->rdY = cy + (b->y - b->dy / 2 - cy) * cos(rotation) + (b->x + b->dx / 2 - cx) * sin(rotation);
	//		
	//		b->ruX = cx + (b->x + b->dx / 2  - cx) * cos(rotation) - (b->y + b->dy / 2 - cy) * sin(rotation);
	//		b->ruY = cy + (b->y + b->dy / 2 - cy) * cos(rotation) + (b->x + b->dx / 2 - cx) * sin(rotation);
	//	}
	//	this->bIsRotated = true;
	//}


} //namespace cusg

