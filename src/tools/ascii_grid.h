#pragma once

/*
Read file in the ARC ASCIIGRID format:
(see http://docs.codehaus.org/display/GEOTOOLS/ArcInfo+ASCII+Grid+format)

Usually the file extension is .asc, but recent versions of ESRI software also recognize the extension .grd.
*/

#include <iostream>
#include <string>
#include <fstream>
#include <cassert>
using namespace std;

#include "ModflowGrid.h"

namespace cusg
{

struct ascii_grid_cell
{
	double x, y, dx, dy; 
};

///bare bone of ascci grid
class ascii_grid
{
	public:

	ascii_grid()
	{
		 ncols=nrows=0;
		 nodata_value=-9999; //the ESRI default is -9999
		 values=NULL;
	}

	//
	//
	// Access functions
	//
	//

	double * getvalues() { return values; }
	double getvalue(unsigned int i, unsigned int j)
	{
		if(i>=nrows || j>=ncols)
		{
			//cerr<<" ! Error: arcinfo_ascii_grid::getvalue index out of bound ("<<i<<"/"<<nrows<<","<<j<<"/"<<ncols<<")"<<endl;
			return nodata_value;
		}

		unsigned int index=i*ncols+j;
		return values[index];
	}

	//access other data members
	unsigned int getncols() { return ncols; }
	unsigned int getnrows() { return nrows; }

	//check if the given valule is valid or not
	bool isvalid(double v){ return v!=nodata_value; }

	//given a 2d point (x,y) returns the index of the point
	virtual Index2d getIndex(double x, double y)=0;

	//
	virtual ascii_grid_cell getCell(unsigned int i, unsigned int j)=0;

	double getNoDATA_value() const { return nodata_value; }

	// get X and get Y
	virtual double getX(unsigned int xid) = 0;

	virtual double getY(unsigned int yid) = 0;
	
protected:

	//The origin of the grid is the upper left and terminus at the lower right.
	unsigned int ncols; //ncols refers to the number of columns in the grid
	unsigned int nrows; //nrows refers to the number of rows in the grid 

	double nodata_value; //nodata_value refers to the value that represents missing data
	                     //the ESRI default is -9999

	double * values;
};

class mf_ascii_grid : public ascii_grid
{
	public:
		
	mf_ascii_grid(ModflowGrid * mfg, double * v, unsigned int col, unsigned int row, double nodata=-9999 ) : ascii_grid()
	{
		values=v;
		ncols=col;
		nrows=row;
		nodata_value=nodata;
		m_mfg=mfg;
	}

	//get the index of the given point
	virtual Index2d getIndex(double x, double y)
	{
		return m_mfg->getIndex(Point2d(x,y));
	}

	//create a cell and return it
	virtual ascii_grid_cell getCell(unsigned int i, unsigned int j)
	{
		Box * box=m_mfg->find_nodeobj(i,j);
		ascii_grid_cell cell;
		cell.dx=box->dx;
		cell.dy=box->dy;
		cell.x=box->x;
		cell.y=box->y;

		return cell;
	}

	virtual double getX(unsigned int xid)
	{
		Box * box = m_mfg->find_nodeobj(xid, 0);
		return box->x;
	}
	virtual double getY(unsigned int yid)
	{
		Box * box = m_mfg->find_nodeobj(0, yid);
		return box->x;
	}

private:

	ModflowGrid * m_mfg;
};

class arcinfo_ascii_grid : public ascii_grid
{
public:

	arcinfo_ascii_grid() : ascii_grid()
	{
		 cellsize=xllcorner=yllcorner=0; //=yulcorner=0;
	}

	~arcinfo_ascii_grid() { delete [] values; }

	bool build(const string& filename )
	{
		ifstream fin(filename.c_str());
		if(fin.good()==false)
		{
			cerr<<"! Error: arcinfo_ascii_grid::build cannot open file: ("<<filename<<")"<<endl;
			return false;
		}

		cout<<"    * Read surface elevation from file: "<<filename<<endl;

		bool r=build(fin);

		fin.close();

		return r;
	}

	bool build(istream& out);

	//
	//
	// Access functions
	//
	//

	//get the x coordinate of the CENTER of a cell
	virtual double getX(unsigned int i)
	{
		if(i>=ncols)
		{
			cerr<<" ! Error: arcinfo_ascii_grid::getX index out of bound"<<endl;
			i=ncols-1;
		}

		return xllcorner+cellsize*(i+0.5f);
	}

	//get the y coordinate of the CENTER of a cell
	virtual double getY(unsigned int i)
	{
		if(i>=nrows)
		{
			cerr<<" ! Error: arcinfo_ascii_grid:: getY index out of bound"<<endl;
			i=nrows-1;
		}

		return yllcorner+cellsize*( (nrows-i-1)+0.5f);
	}

	//get the index of the given point
	virtual Index2d getIndex(double x, double y)
	{
		int i=nrows-1-floor((y-yllcorner)/cellsize);
		int j=floor((x-xllcorner)/cellsize);
		
		if(i<0 || j<0 || ((unsigned int)i)>=nrows || ((unsigned int)j)>=ncols)
		{
			if(i<0) i=0;
			if(j<0) j=0;
			if(i>=(int)nrows) i=(int)(nrows-1);
			if(j>=(int)ncols) j=(int)(ncols-1);
		}

		return Index2d(i,j);
	}

	//create a cell and return it
	virtual ascii_grid_cell getCell(unsigned int i, unsigned int j)
	{
		ascii_grid_cell cell;
		cell.x=getX(j);
		cell.y=getY(i);
		cell.dx=cellsize;
		cell.dy=cellsize;

		return cell;
	}

	//access other data members

	double getxllcorner(){ return xllcorner; }
	double getyllcorner(){ return yllcorner; }
	double getcellsize(){ return cellsize; }

protected:

	string tolower(const string& s)
	{
		string lower;
		for(string::const_iterator i=s.begin();i!=s.end();i++)
		{
			char c=::tolower(*i);
			lower.push_back(c);
		}
		return lower;
	}

	//The origin of the grid is the upper left and terminus at the lower right.

	//xllcorner and yllcorner are given as the EDGES of the grid, NOT the centers of the edge cells.
	double xllcorner; //xllcorner refers to the western edge of the grid 
	double yllcorner;  //yllcorner refers to the southern edge of the grid 

	//determined from the input values
	//double yulcorner;  //yllcorner refers to the north edge of the grid 

	double cellsize;  //cellsize refers to the resolution of the grid 
};

//perform bilinear interpolation between (i,j), (i+1,j), (i,j+1), (i+1,j+1)
//x and y must be between 0 and 1
double bilinear_interpolation(ascii_grid * grid, unsigned int i, unsigned int j, double x, double y);


//perform trilinear interpolation between (i,j), (i+1,j), (i,j+1), (i+1,j+1) and grid1 and grid2
//x and y and z must be between 0 and 1
double trilinear_interpolation(ascii_grid * grid1, ascii_grid * grid2, unsigned int i, unsigned int j, double x, double y, double z);


}//end of namespace cusg