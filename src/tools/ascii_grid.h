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

struct AsciiGridCell
{
    double x, y, dx, dy;
};

///bare bones of an ASCII grid
class AsciiGrid
{
public:
    AsciiGrid()
    {
        m_ncols = m_nrows = 0;
        m_nodata_value = -9999; //the ESRI default is -9999
        m_values = nullptr;
    }

    //
    //
    // Access functions
    //
    //

    double * get_values() { return m_values; }
    double get_value(unsigned int i, unsigned int j)
    {
        if (i >= m_nrows || j >= m_ncols)
        {
            return m_nodata_value;
        }

        unsigned int index = i * m_ncols + j;
        return m_values[index];
    }

    //access other data members
    unsigned int getncols() { return m_ncols; }
    unsigned int getnrows() { return m_nrows; }

    //check if the given value is valid or not
    bool is_valid(double value) { return value != m_nodata_value; }

    //given a 2d point (x,y) returns the index of the point
    virtual Index2d get_index(double x, double y) = 0;

    //
    virtual AsciiGridCell get_cell(unsigned int i, unsigned int j) = 0;

    double get_nodata_value() const { return m_nodata_value; }

    // get X and get Y
    virtual double get_x(unsigned int col) = 0;

    virtual double get_y(unsigned int row) = 0;

protected:
    //The origin of the grid is the upper left and terminus at the lower right.
    unsigned int m_ncols; //ncols refers to the number of columns in the grid
    unsigned int m_nrows; //nrows refers to the number of rows in the grid

    double m_nodata_value; //nodata_value refers to the value that represents missing data
                           //the ESRI default is -9999

    double * m_values;
};

class ModflowAsciiGrid : public AsciiGrid
{
public:
    ModflowAsciiGrid(ModflowGrid * mfgrid, double * grid_values, unsigned int col, unsigned int row,
                     double nodata = -9999)
        : AsciiGrid()
    {
        m_values = grid_values;
        m_ncols = col;
        m_nrows = row;
        m_nodata_value = nodata;
        m_mfgrid = mfgrid;
    }

    //get the index of the given point
    Index2d get_index(double x, double y) override { return m_mfgrid->get_index(Point2d(x, y)); }

    //create a cell and return it
    AsciiGridCell get_cell(unsigned int i, unsigned int j) override
    {
        Box * box = m_mfgrid->find_box(i, j);
        AsciiGridCell cell;
        cell.dx = box->dx;
        cell.dy = box->dy;
        cell.x = box->x;
        cell.y = box->y;

        return cell;
    }

    double get_x(unsigned int col) override
    {
        Box * box = m_mfgrid->find_box(col, 0);
        return box->x;
    }
    double get_y(unsigned int row) override
    {
        Box * box = m_mfgrid->find_box(0, row);
        return box->x;
    }

private:
    ModflowGrid * m_mfgrid;
};

class ArcinfoAsciiGrid : public AsciiGrid
{
public:
    ArcinfoAsciiGrid() : AsciiGrid()
    {
        m_cellsize = m_xllcorner = m_yllcorner = 0; //=yulcorner=0;
    }

    ~ArcinfoAsciiGrid() { delete[] m_values; }

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
    double get_x(unsigned int i) override
    {
        //an index past the last column is taken as the last
        if (i >= m_ncols) i = m_ncols - 1;

        return m_xllcorner + m_cellsize * (i + 0.5f);
    }

    //get the y coordinate of the CENTER of a cell
    double get_y(unsigned int i) override
    {
        //an index past the last row is taken as the last
        if (i >= m_nrows) i = m_nrows - 1;

        return m_yllcorner + m_cellsize * ((m_nrows - i - 1) + 0.5f);
    }

    //get the index of the given point
    Index2d get_index(double x, double y) override
    {
        int i = m_nrows - 1 - floor((y - m_yllcorner) / m_cellsize);
        int j = floor((x - m_xllcorner) / m_cellsize);

        if (i < 0 || j < 0 || ((unsigned int)i) >= m_nrows || ((unsigned int)j) >= m_ncols)
        {
            if (i < 0) i = 0;
            if (j < 0) j = 0;
            if (i >= (int)m_nrows) i = (int)(m_nrows - 1);
            if (j >= (int)m_ncols) j = (int)(m_ncols - 1);
        }

        return Index2d(i, j);
    }

    //create a cell and return it
    AsciiGridCell get_cell(unsigned int i, unsigned int j) override
    {
        AsciiGridCell cell;
        cell.x = get_x(j);
        cell.y = get_y(i);
        cell.dx = m_cellsize;
        cell.dy = m_cellsize;

        return cell;
    }

    //access other data members

    double get_xllcorner() { return m_xllcorner; }
    double get_yllcorner() { return m_yllcorner; }
    double get_cellsize() { return m_cellsize; }

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
    double m_xllcorner; //xllcorner refers to the western edge of the grid
    double m_yllcorner; //yllcorner refers to the southern edge of the grid

    //determined from the input values
    //double yulcorner;  //yllcorner refers to the north edge of the grid

    double m_cellsize; //cellsize refers to the resolution of the grid
};

//perform bilinear interpolation between (i,j), (i+1,j), (i,j+1), (i+1,j+1)
//x and y must be between 0 and 1
double bilinear_interpolation(AsciiGrid * grid, unsigned int i, unsigned int j, double x, double y);

} //end of namespace cusg