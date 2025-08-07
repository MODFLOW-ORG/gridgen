
/*!

\brief ...

\author <a href="masc.cs.gmu.edu/">MASC group</a>, George Mason University,
        <a href="profile.usgs.gov/langevin/">Christian Langevin</a>, USGS
\bug    No known bugs.

*/

#include "ModflowGrid2D.h"
#include <climits>
#define IndexEPSION 1e-10

namespace cusg
{

	
	/*! \brief Instantiate a ModflowGrid2D object.
        
		\param nrow The number of rows
        \param *ncol* The number of columns.
		\param *delr*: The column spacings along a row.
		\param *delc*: The row spacings along a column.
            
        Optional arguments that can be specified using "...":

        \param *xoffset*: The offset of the grid in the x direction.

        \param *yoffset*: The offset of the grid in the y direction.

        \param *rotation*: The grid rotation angle in degrees relative to the lower  
                      left corner.  DOES NOT WORK YET
                      
        The local coordinate system has a (0, 0) x,y location that corresponds 
        to the lower left corner of the grid.
        
        \note This is a variadic function.
    */
	ModflowGrid2D::ModflowGrid2D(int nrow, int ncol, double delr, double delc, double xoffset, double yoffset, double rotation)
	{
        //log the creation
        cout<<"- A ModflowGrid2D object was instantiated."<<endl;
		
        //assign the arguments to the object.
        this->nlay = 1;
        this->nrow = nrow;
        this->ncol = ncol;
        this->nodes = nrow * ncol;
        this->delr = new double[ncol]; //ModflowArray(delr, arrayshape=(ncol,)).as_ndarray();
        this->delc = new double[nrow]; //ModflowArray(delc, arrayshape=(nrow,)).as_ndarray();
        
    	assert(this->delr && this->delc);
    	for(int i=0;i<ncol;i++) this->delr[i]=delr;
    	for(int i=0;i<nrow;i++) this->delc[i]=delc;

    	height=delc*nrow;
    	width=delr*ncol;

        //set defaults
        this->xoffset = xoffset;
        this->yoffset = yoffset;
        this->rotation = rotation;

        init2D();
	}
	
    /*! \brief Instantiate a ModflowGrid2D object.

        \param nrow The number of rows
        \param *ncol* The number of columns.
        \param *delr*: The column spacings along a row.
        \param *delc*: The row spacings along a column.

        Optional arguments that can be specified using "...":

        \param *xoffset*: The offset of the grid in the x direction.

        \param *yoffset*: The offset of the grid in the y direction.

        \param *rotation*: The grid rotation angle in degrees relative to the lower
                      left corner.  DOES NOT WORK YET

        The local coordinate system has a (0, 0) x,y location that corresponds
        to the lower left corner of the grid.

        \note This is a variadic function.
    */
    ModflowGrid2D::ModflowGrid2D(int nrow, int ncol, double * delr, double * delc, double xoffset, double yoffset, double rotation)
    {
        //log the creation
        cout<<"- A ModflowGrid2D object was instantiated."<<endl;

        //assign the arguments to the object.
        this->nlay = 1;
        this->nrow = nrow;
        this->ncol = ncol;
        this->nodes = nrow * ncol;
        this->delr = new double[ncol]; //ModflowArray(delr, arrayshape=(ncol,)).as_ndarray();
        this->delc = new double[nrow]; //ModflowArray(delc, arrayshape=(nrow,)).as_ndarray();

        assert(this->delr && this->delc);

        height=width=0;
        for(int i=0;i<ncol;i++){ this->delr[i]=delr[i]; width+=delr[i]; }
        for(int i=0;i<nrow;i++){ this->delc[i]=delc[i]; height+=delc[i]; }

        //set defaults
        this->xoffset = xoffset;
        this->yoffset = yoffset;
        this->rotation = rotation;

        init2D();
    }


    ModflowGrid2D::ModflowGrid2D() : Grid()
    {
        /* be sure to call init2D() once you set up the grid! */
        nlay = 0;
        nlay=nrow=ncol=nodes=0;
        delr=delc=NULL;
        Xe=Ye=NULL; //x/y coord of all edges
        xoffset = 0;
        yoffset = 0;
        rotation= 0;
        width=height=0;
    }
	
	/// # of nodes in all layers (int *)
	int * ModflowGrid2D::get_nodelay()
	{
        int * nodelay = new int[nlay]; 
        assert(nodelay);
        
        for(int k=0; k<nlay; k++)
            nodelay[k] = this->nrow * this->ncol;

        return nodelay;
	}
	
	///x coords for the center of all cells
    double * ModflowGrid2D::get_local_x_array()
    {
    	//allocate space
    	double * x = new double[ncol];
    	assert(x);
    	//accumulate the values of delr
    	for(int i=0;i<ncol;i++) x[i]=((i==0)?(0):x[i-1])+delr[i]; 
        for(int i=0;i<ncol;i++) x[i]=x[i]+xoffset-0.5*delr[i]; 
        
        return x;  
    }
    
    
    /// #  y coord for the center of the cell (and for all cells)
    double * ModflowGrid2D::get_local_y_array()
    {
  //  	//allocate space
  //  	double * y = new double[nrow];
  //  	assert(y);
		//double Ly=0;
  //  	for(int i=0;i<nrow;i++) Ly+=delc[i];
  //  	for(int i=0;i<nrow;i++) y[i] = ((i==0)?(0):y[i-1]+delc[i]); 
		////modified by Guilin
  //  	//for(int i=0;i<nrow;i++) y[i] = Ly-y[i]-0.5*delc[i]+yoffset;
  //      for(int i=0;i<nrow;i++) y[i] = Ly-y[i]+0.5*delc[i]+yoffset;

	//alocate space
		double * y =  new double[nrow];
		assert(y);
		for (int i = nrow - 1; i >=0; i--) y[i] = ((i ==nrow-1)?(0):y[i + 1]) +delc[i];
		for (int i = 0; i < nrow; i++) y[i] = y[i] + yoffset - 0.5*delc[i];

        return y;  
    }
    
	/// x coord of all edges
    double * ModflowGrid2D::get_local_Xe_array()
    {
        //The Xe array contains the edge coordinate for each cell
        int size=ncol+1;
        double * Xe = new double[size];
        assert(Xe);
        Xe[0]=0;
        for(int i=1;i<size;i++) Xe[i] = Xe[i-1]+delr[i-1]; 
        for(int i=0;i<size;i++) Xe[i] += xoffset;

        return Xe;
	}
	
	/// y coord of all edges
    double * ModflowGrid2D::get_local_Ye_array()
    {
        //allocate space
        int size=nrow+1;
    	double * Ye = new double[size];
    	assert(Ye);
		Ye[size - 1] = 0;
		for (int i = size - 2; i >=0; i--) Ye[i] = Ye[i + 1] + delc[i];
		for (int i = 0; i < size; i++) Ye[i] +=yoffset;

		//double Ly=0;
		//for(int i=0;i<nrow;i++) Ly+=delc[i];
		//Ye[0]=0;
		//for(int i=1;i<size;i++) Ye[i] = Ye[i-1] + delc[i-1];
  //      for(int i=1;i<size;i++) Ye[i] = Ly - Ye[i];
  //      Ye[0] = Ly;
		//for(int i=0;i<size;i++) Ye[i]+=yoffset;
		
		return Ye;
    }
    
    
	/// Return a numpy array of size (nodes) that contains the area for all cells.
	double * ModflowGrid2D::get_cell_areas()
	{
		//log the creation
        cout<<"Getting cell areas."<<endl;
        
        double * area=new double(nrow*ncol);
        assert(area);
        
        for( int r=0;r<nrow;r++)
        {
            for(int c=0;c<ncol;c++)
            {
                int nodeid = get_nodeid(r, c);
                area[nodeid] = delr[c] * delc[r];
            }
        }        
        
        return area;
	}


        
	/*!
        \brief Return the four vertices for the specified nodeid.

        3-------2
        |       |
        |       |
        0-------1
    */

	vector<Point2d> ModflowGrid2D::get_vertices(int nodeid)
	{
        Index2d ij = get_indices(nodeid);
        int i=ij[0];
        int j=ij[1];
                
        Point2d p0(Xe[j], Ye[i + 1]);
        Point2d p1(Xe[j + 1], Ye[i + 1]);
        Point2d p2(Xe[j + 1], Ye[i]);
        Point2d p3(Xe[j], Ye[i]);
        
        vector<Point2d> vertices;
        vertices.push_back(p0);
        vertices.push_back(p1);
        vertices.push_back(p2);
        vertices.push_back(p3);
                         
        return vertices;
	}

	/*! 
        \brief Return a tuple containing two points (the lower left and the upper 
        right) in global grid coordinates.
    */    	
	pair<Point2d,Point2d> ModflowGrid2D::get_extent()
	{
		Point2d min, max;
		
        min[0] = Xe[0];
        min[1] = Ye[nrow];
        
        max[0] = Xe[ncol];
        max[1] = Ye[0];        
        
        return make_pair(min,max);
	}
	
	///Return the nodeid for the modflow grid.
	int ModflowGrid2D::get_nodeid(int i, int j)
	{
        return i * ncol + j;
    }
    
    ///Return the indices for the specified nodeid.
    Index2d ModflowGrid2D::get_indices(int nodeid)
    {
		
        if( nodeid > nrow * ncol )
        {
            cerr<<"Nodeid ("<<nodeid<<" is not valid..."<<endl;
            assert(false);
        }

        int i = int( nodeid / this->ncol );
        int j = nodeid - i * this->ncol;
        return Index2d(i, j);
    }
        
	///get box from id
    Box * ModflowGrid2D::create_nodeobj(int nodeid)
    {
        Index2d ij = get_indices(nodeid);
        int i=ij[0];
        int j=ij[1];
        return create_nodeobj(i,j);
    }

    Box * ModflowGrid2D::create_nodeobj(int i, int j)
    {
        double x = this->X[j];
        double y = this->Y[i];
        Point3d position(x, y,0);
        
        double dx = this->delr[j];
        double dy = this->delc[i];
        Vector3d dxdydz(dx, dy, 0);
        
        Box * nodeobj = new Box(position, dxdydz); //(position, dxdydz, save=False);
        
        assert(nodeobj);
        
        return nodeobj;
    }
    
    ///
    Box * ModflowGrid2D::find_nodeobj(const Point2d& pt)
    {
        Index2d id2d=getIndex(pt);

        if(id2d[0]==INT_MAX || id2d[1]==INT_MAX|| id2d[0] ==INT_MIN || id2d[1]==INT_MIN) return NULL;

        //get the box
        int id=get_nodeid(id2d[0],id2d[1]);
        return nodegroup[id];
    }

	/// return a box containtain the given id
    Box * ModflowGrid2D::find_nodeobj(int i, int j)
	{
        //get the box
        int id=get_nodeid(i,j);

		if(nodegroup.empty() || id<0 || ((unsigned)id)>nodegroup.size())
		{
			cerr<<"! Error: ModflowGrid2D::find_nodeobj index out of bound"<<endl;
			//return NULL;
			exit(1);
		}

        return nodegroup[id];
	}

    ///get box connections from id
    ///Return a sorted list of connections, where a connection is
    ///a list [tonodeid, fldir].  fldir is -1, +1, -2, +2 for -x, +x, -y, +y.  
    list< pair<int,int> > ModflowGrid2D::get_node_connections(int nodeid)
    {
        list< pair<int,int> > connections;
    	Index2d ij=get_indices(nodeid);
    	int i=ij[0];
    	int j=ij[1];
    	
    	//back
    	if(i>0)
    	{
            int nn = get_nodeid(i-1, j);
            connections.push_back(make_pair(nn, 2));
    	}
 		//left
        if(j > 0)
        {
            int nn = get_nodeid(i, j - 1);
            connections.push_back(make_pair(nn, -1));
        }
    	//right
        if(j < ncol - 1)
        {
            int nn = get_nodeid(i, j + 1);
            connections.push_back(make_pair(nn, 1));
    	}
        //front
        if( i < nrow - 1 )
        {
            int  nn = get_nodeid(i + 1, j);
            connections.push_back( make_pair(nn, -2) );
        }   
        
        return connections;
        
    }
    
    //	/*!
    //	    Method for intersecting the model grid with a point, line, rectangle,
    //        or polygon.  A rectangle intersection is faster than a polygon
    //        intersection.
    //     */
    //    void intersection()
    //    {
    //		//nothing yet
    //    }

    ///number the nodes (i.e. leaves) in the grid
    void ModflowGrid2D::number_nodes()
    {
        unsigned int number=0;
        for(int r=0;r<nrow;r++)
          {
              for(int c=0;c<ncol;c++)
              {
                  int bid=get_nodeid(r,c);
                  Box * box=nodegroup[bid];
                  if(box->active)
                      box->number=number++;
                  else
                      box->number=-1;
              }//end c
        }//end r
    }

    //access the value of the values in bottom
    double ModflowGrid2D::botm(int layer, int row, int col)
    {
        int id=layer*(ncol*nrow)+row*ncol+col;
        return bot[id];
    }

    void ModflowGrid2D::init2D() {
        X = get_local_x_array();
        Y = get_local_y_array();
        Xe = get_local_Xe_array();
        Ye = get_local_Ye_array();

        //create and connect actual nodes (boxes)
        //create_nodes();
    }

    ///create nodes and add them to nodegroup
    void ModflowGrid2D::create_nodes()
    {
        // Create the nodegroup.
        nodegroup = NodeGroup(this);
        nodegroup.nlay=1;

        for( int i=0; i< nrow; i++)
        {
            double y = Y[i];
            for( int j=0; j< ncol; j++)
            {
                double x = X[j];

                Point3d position;
                Vector3d dxdydz;

                double z = 0.5 * (botm(0,i,j)+botm(1,i,j));
                double dz = (botm(0,i,j) - botm(1,i,j));

                position = Point3d(x, y, z);
                dxdydz = Vector3d(delr[j], delc[i], dz);
                Box * nodeobj = new Box(position, dxdydz);
                assert(nodeobj);
                nodeobj->layer=0;
                nodegroup.add(nodeobj);
            }//end i
        }//end j

        add_connections();
    }

    ///add connection between boxes
    ///fldir is -1, +1, -2, +2, -3, +3 for -x, +x, -y, +y, -z, +z.
    void ModflowGrid2D::add_connections()
    {

        //grid.logger.info('Adding the 2d level zero connections.')
        //cout<<"Adding the 2d level zero connections"<<endl;

        for( int i=0; i<nrow; i++)
        {
            for( int j=0; j<ncol; j++)
            {
                int nodeid = get_nodeid(i,j);

                //cout<<"Adding connections for nodeid: "<< nodeid<<endl;

                Box * nodeobj = nodegroup[nodeid];

                if( i < nrow - 1)
                {
                    //make connection to front
                    int tonode = get_nodeid(i + 1, j);
                    int fldir = -2;
                    nodeobj->addConnection(nodegroup[tonode], fldir);
                }

                if( i > 0 )
                {
                    //make connection to back
                    int tonode = get_nodeid(i - 1, j);
                    int fldir = 2;
                    nodeobj->addConnection(nodegroup[tonode], fldir);
                }
                if( j < ncol - 1 )
                {
                    //make connection to right
                    int tonode = get_nodeid(i, j + 1);
                    int fldir = 1;
                    nodeobj->addConnection(nodegroup[tonode], fldir);
                }
                if( j > 0 )
                {
                    //make connection to left
                    int tonode = get_nodeid(i, j - 1);
                    int fldir = -1;
                    nodeobj->addConnection(nodegroup[tonode], fldir);
                }
            }//end j
        }//end i
    }


    ///find the index of the box containing pt
	///row major (row_id, col_id)
    Index2d ModflowGrid2D::getIndex(const Point2d& pt)
    {
        int i_x=INT_MAX;
        int i_y=INT_MAX;
		bool bContain = true;

		double * Xe=this->Xe;
		double * Ye=this->Ye;

		bool xAsc = delr[0] >0;
		bool yAsc = delc[0] >0;

		if(xAsc &&(pt[0] <Xe[0] - IndexEPSION || pt[0] >Xe[ncol] + IndexEPSION))
		{
			//assert(Xe[0] < Xe[ncol]);
			if(pt[0] <Xe[0] - IndexEPSION)
				i_x = INT_MIN;
			else if(pt[0] > Xe[ncol] + IndexEPSION)
				i_x = INT_MAX;
			bContain = false;
		}
		if( (!xAsc)&& (pt[0] <Xe[ncol] - IndexEPSION || pt[0] >Xe[0] + IndexEPSION))
		{
			//assert(Xe[0] > Xe[ncol] );
			if(pt[0] >Xe[0]+ IndexEPSION)
				i_x = INT_MIN;
			else if(pt[0] <Xe[ncol] - IndexEPSION)
				i_x = INT_MAX;
			bContain = false;
		}
		if(yAsc &&(pt[1]>Ye[0] + IndexEPSION|| pt[1] <Ye[nrow] - IndexEPSION))
		{
			//assert(Ye[0] > Ye[nrow] );
			if(pt[1] > Ye[0]+ IndexEPSION)
				i_y = INT_MIN;
			else if(pt[1] <Ye[nrow] - IndexEPSION)
				i_y = INT_MAX;
			bContain = false;
		}
		if( (!yAsc) && (pt[1] <Ye[0]  - IndexEPSION || pt[1] >Ye[nrow] + IndexEPSION))
		{
			//assert(Ye[0] < Ye[nrow]);
			if(pt[1] <Ye[0] - IndexEPSION)
				i_y = INT_MIN;
			else if(pt[1] >Ye[nrow] + IndexEPSION)
				i_y = INT_MAX;
			bContain = false;
		}

		if(!bContain)
		{
			return  Index2d(i_y,i_x);
		}
		if(yAsc)
		{
			//cerr<<"The coordinate order is not correct!\n";
			assert(Ye[0] > Ye[nrow]);
		}

		//use the binary search
		i_x = getXIndexByBinarySearch(pt, Xe, 0, this->ncol, xAsc);
		i_y = getYIndexByBinarySearch(pt, Ye, 0, this->nrow, yAsc);


		////find i_x
		//for (int i = 0; i <this->ncol; i++)
		//{
		//	int n = i+1;
		//	if(xAsc && (Xe[i] <=pt[0] && Xe[n] >=pt[0]))
		//	{
		//		i_x = i; break;
		//	}
		//	if((!xAsc) && (Xe[n]<=pt[0] && Xe[i] >=pt[0]))
		//	{
		//		i_x = i; break;
		//	}
		//}
		////find i_y
		//for (int i = 0; i <this->nrow; i++)
		//{
		//	int n = i+1;
		//	if(yAsc &&(Ye[i]>=pt[1] && Ye[n]<=pt[1]))
		//	{
		//		i_y = i; break;
		//	}
		//	if((!yAsc) && (Ye[i] <=pt[1] && Ye[n] >=pt[1]))
		//	{
		//		i_y = i; break;
		//	}
		//}
		//cout<<i_y<<" - "<<i_x<<"\t";

		return Index2d(i_y,i_x);
		/******************************************************************/
		//the below code is by Jyh-Ming Lien
        //find i_x
        //double * Xe=this->Xe;
        //for(int i=0;i<this->ncol;i++)
        //{
        //    int n=i+1;
        //    if(Xe[i]<=pt[0] && Xe[n]>=pt[0]){
        //        i_x=i;
        //        break;
        //    }
        //}

        ////find i_y
        //double * Ye=this->Ye;
        //for(int i=0;i<this->nrow;i++)
        //{
        //    int n=i+1;
        //    if(Ye[i]>=pt[1] && Ye[n]<=pt[1]){//this is not correct since the order is from n-1 to 0 in y, instead of from 0 to n-1 in x.
        //        i_y=i;
        //        break;
        //    }
        //}
    }

	int ModflowGrid2D::getXIndexByBinarySearch(const Point2d& pt, double * Xe, int left, int right, bool ascending)
	{
		if(left == right)
			return left;
		if(left + 1 == right)
			return left;

		int mid = (left + right) / 2;
		if(Xe[mid] == pt[0])
			return mid;
		else if(Xe[mid] > pt[0])
		{
			if(ascending) 
				return getXIndexByBinarySearch(pt, Xe, left, mid, ascending);
			else
				return getXIndexByBinarySearch(pt, Xe, mid, right, ascending);
		}
		else // if(Xe[mid] <pt[0])
		{
			if(!(Xe[mid] < pt[0]))
				cerr<<"! Error: ModflowGrid2D::getXIndexByBinarySearch may suffer precision problem!\n";
			if(ascending)
				return getXIndexByBinarySearch(pt, Xe, mid, right, ascending);
			else
				return getXIndexByBinarySearch(pt, Xe, left, mid, ascending);
		}
		//else
			//return mid;

		//return -1;
	}


	int  ModflowGrid2D::getYIndexByBinarySearch(const Point2d& pt, double*  Ye, int left, int right, bool ascending)
	{
		if(left == right)
			return left;
		if(left + 1 == right)
			return left;

		int mid = (left + right) / 2;
		if(Ye[mid] == pt[1])
			return mid;
		else if(Ye[mid] > pt[1])
		{
			if(ascending)
				return getYIndexByBinarySearch(pt, Ye, mid, right, ascending);
			else 
				return getYIndexByBinarySearch(pt, Ye, left, mid, ascending);
		}
		else //if(Ye[mid] <pt[1])
		{
			if(!(Ye[mid] <pt[1]))
				cerr<<"Error: ModflowGrid2D::getYIndexByBinarySearch may have precision problem!\n";
			if(ascending)
				return getYIndexByBinarySearch(pt, Ye, left, mid, ascending);
			else 
				return getYIndexByBinarySearch(pt, Ye, mid, right, ascending);
		}

		//return -1;
	}

} //namespace cusg

