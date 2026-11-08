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
#include "GridgenError.h"

namespace cusg
{

struct ascii_grid_cell
{
    double x, y, dx, dy;
};

///bare bones of an ASCII grid
class ascii_grid
{
public:
    ascii_grid()
    {
        ncols = nrows = 0;
        nodata_value = -9999; //the ESRI default is -9999
        values = nullptr;
    }

    //
    //
    // Access functions
    //
    //

    double * getvalues() { return values; }
    double getvalue(unsigned int i, unsigned int j)
    {
        if (i >= nrows || j >= ncols)
        {
            return nodata_value;
        }

        unsigned int index = i * ncols + j;
        return values[index];
    }

    //access other data members
    unsigned int getncols() { return ncols; }
    unsigned int getnrows() { return nrows; }

    //check if the given value is valid or not
    bool isvalid(double value) { return value != nodata_value; }

    //given a 2d point (x,y) returns the index of the point
    virtual Index2d getIndex(double x, double y) = 0;

    //
    virtual ascii_grid_cell getCell(unsigned int i, unsigned int j) = 0;

    double getNoDATA_value() const { return nodata_value; }

    // get X and get Y
    virtual double getX(unsigned int col) = 0;

    virtual double getY(unsigned int row) = 0;

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
    mf_ascii_grid(ModflowGrid * mfgrid, double * grid_values, unsigned int col, unsigned int row,
                  double nodata = -9999)
        : ascii_grid()
    {
        values = grid_values;
        ncols = col;
        nrows = row;
        nodata_value = nodata;
        m_mfgrid = mfgrid;
    }

    //get the index of the given point
    Index2d getIndex(double x, double y) override { return m_mfgrid->getIndex(Point2d(x, y)); }

    //create a cell and return it
    ascii_grid_cell getCell(unsigned int i, unsigned int j) override
    {
        Box * box = m_mfgrid->find_box(i, j);
        ascii_grid_cell cell;
        cell.dx = box->dx;
        cell.dy = box->dy;
        cell.x = box->x;
        cell.y = box->y;

        return cell;
    }

    double getX(unsigned int col) override
    {
        Box * box = m_mfgrid->find_box(col, 0);
        return box->x;
    }
    double getY(unsigned int row) override
    {
        Box * box = m_mfgrid->find_box(0, row);
        return box->x;
    }

private:
    ModflowGrid * m_mfgrid;
};

class arcinfo_ascii_grid : public ascii_grid
{
public:
    arcinfo_ascii_grid() : ascii_grid()
    {
        cellsize = xllcorner = yllcorner = 0; //=yulcorner=0;
    }

    ~arcinfo_ascii_grid() { delete[] values; }

    bool build(const string & filename)
    {
        ifstream fin(filename.c_str());
        if (fin.good() == false) throw GridgenError("cannot open ASCII grid file " + filename);

        cout << "    * Read surface elevation from file: " << filename << endl;

        bool built = build(fin);

        fin.close();

        return built;
    }

    bool build(istream & out);

    //
    //
    // Access functions
    //
    //

    //get the x coordinate of the CENTER of a cell
    double getX(unsigned int i) override
    {
        //an index past the last column is taken as the last
        if (i >= ncols) i = ncols - 1;

        return xllcorner + cellsize * (i + 0.5f);
    }

    //get the y coordinate of the CENTER of a cell
    double getY(unsigned int i) override
    {
        //an index past the last row is taken as the last
        if (i >= nrows) i = nrows - 1;

        return yllcorner + cellsize * ((nrows - i - 1) + 0.5f);
    }

    //get the index of the given point
    Index2d getIndex(double x, double y) override
    {
        int i = nrows - 1 - floor((y - yllcorner) / cellsize);
        int j = floor((x - xllcorner) / cellsize);

        if (i < 0 || j < 0 || ((unsigned int)i) >= nrows || ((unsigned int)j) >= ncols)
        {
            if (i < 0) i = 0;
            if (j < 0) j = 0;
            if (i >= (int)nrows) i = (int)(nrows - 1);
            if (j >= (int)ncols) j = (int)(ncols - 1);
        }

        return Index2d(i, j);
    }

    //create a cell and return it
    ascii_grid_cell getCell(unsigned int i, unsigned int j) override
    {
        ascii_grid_cell cell;
        cell.x = getX(j);
        cell.y = getY(i);
        cell.dx = cellsize;
        cell.dy = cellsize;

        return cell;
    }

    //access other data members

    double getxllcorner() { return xllcorner; }
    double getyllcorner() { return yllcorner; }
    double getcellsize() { return cellsize; }

protected:
    string tolower(const string & s)
    {
        string lower;
        for (char i : s)
        {
            char c = ::tolower(i);
            lower.push_back(c);
        }
        return lower;
    }

    //The origin of the grid is the upper left and terminus at the lower right.

    //xllcorner and yllcorner are given as the EDGES of the grid, NOT the centers of the edge cells.
    double xllcorner; //xllcorner refers to the western edge of the grid
    double yllcorner; //yllcorner refers to the southern edge of the grid

    //determined from the input values
    //double yulcorner;  //yllcorner refers to the north edge of the grid

    double cellsize; //cellsize refers to the resolution of the grid
};

//perform bilinear interpolation between (i,j), (i+1,j), (i,j+1), (i+1,j+1)
//x and y must be between 0 and 1
double bilinear_interpolation(ascii_grid * grid, unsigned int i, unsigned int j, double x,
                              double y);

} //end of namespace cusg