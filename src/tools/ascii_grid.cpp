#include "ascii_grid.h"

#include <cctype>

namespace cusg
{

/*

Read file in the following format:
(see http://docs.codehaus.org/display/GEOTOOLS/ArcInfo+ASCII+Grid+format)

ncols 157
nrows 171
xllcorner -156.08749650000
yllcorner 18.870890200000
cellsize 0.00833300
0 0 1 1 1 2 3 3 5 6 8 9 12 14 18 21 25 30 35 41 47 53
59 66 73 79 86 92 97 102 106 109 112 113 113 113 111 109 106
103 98 94 89 83 78 72 67 61 56 51 46 41 37 32 29 25 22 19
etc...

*/


bool arcinfo_ascii_grid::build(istream& ins)
{
	//read header
	string tmp;
	for(int i=0;i<6;i++) //loop 6 times at most
	{
		if( isdigit(ins.peek()) ) break;

		ins>>tmp;
		
		tmp=tolower(tmp);

		if(tmp=="ncols") ins>>ncols;
		if(tmp=="nrows") ins>>nrows;
		if(tmp=="xllcorner") ins>>xllcorner;
		if(tmp=="yllcorner") ins>>yllcorner;
		if(tmp=="cellsize") ins>>cellsize;
		if(tmp=="nodata_value") ins>>nodata_value;
	}

	if(cellsize<=0){ cerr<<"! Error: cellsize<=0"<<endl; return false; }
	if(ncols<=0){ cerr<<"! Error: ncols<=0"<<endl; return false; }
	if(nrows<=0){ cerr<<"! Error: nrows<=0"<<endl; return false; }
	
	//allocate mem
	unsigned int size=nrows*ncols; 
	values=new double[size];
	assert(values);

	//read values
	for(unsigned int i=0;i<size;i++)
	{
		//int id=(nrows-i/ncols-1)*(ncols)+(i%ncols);
		ins>>values[i];
	}

	//determine yulcorner
	//yulcorner=yllcorner-nrows*cellsize;

	return true;
}

//perform bilinear interpolation between (i,j), (i+1,j), (i,j+1), (i+1,j+1)
//x and y must be between 0 and 1
double bilinear_interpolation(ascii_grid * grid, unsigned int i, unsigned int j, double x, double y)
{
	//get values
	double ul=grid->getvalue(i,j);
	double ur=grid->getvalue(i,j+1);
	double lr=grid->getvalue(i+1,j+1);
	double ll=grid->getvalue(i+1,j);

	//check the validness of the values...
	if( grid->isvalid(ul)==false && grid->isvalid(ur)==false && grid->isvalid(ll)==false && grid->isvalid(lr)==false)
		return ul; //All values are invalid

	if( grid->isvalid(ul) && grid->isvalid(ur)==false ){ ur=ul; }
	if( grid->isvalid(ul)==false && grid->isvalid(ur) ){ ul=ur; }
	if( grid->isvalid(lr) && grid->isvalid(ll)==false ){ ll=lr; }
	if( grid->isvalid(lr)==false && grid->isvalid(ll) ){ lr=ll; }

	//interpolate in x dir
	double x1, x2;

	if( grid->isvalid(ul)==false && grid->isvalid(ur)==false) x1=ul; //avoid interpolating invalid values...
	else x1=(1-x)*ul+x*ur;

	if( grid->isvalid(ll)==false && grid->isvalid(lr)==false) x2=ll; //avoid interpolating invalid values...
	else x2=(1-x)*ll+x*lr;
	//else x2=(1-x)*ll+x*lr;

	if( grid->isvalid(x1) && grid->isvalid(x2)==false ){ x2=x1;  }
	if( grid->isvalid(x1)==false && grid->isvalid(x2) ){ x1=x2;  }

	//interpolate in y dir
	return (1-y)*x1+y*x2;
}


//perform trilinear interpolation between (i,j), (i+1,j), (i,j+1), (i+1,j+1) and grid1 and grid2
//x and y and z must be between 0 and 1
double trilinear_interpolation(ascii_grid * grid1, ascii_grid * grid2, unsigned int i, unsigned int j, double x, double y, double z)
{
	double z1=bilinear_interpolation(grid1,i,j,x,y);
	double z2=bilinear_interpolation(grid2,i,j,x,y);
	
	//interpolate in z dir
	if( grid1->isvalid(z1) && grid2->isvalid(z2)==false ) z2=z1;
	if( grid1->isvalid(z1)==false && grid2->isvalid(z2) ) z1=z2;

	return (1-z)*z1+z*z2;
}


} //end of namespace cusg