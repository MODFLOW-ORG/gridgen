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

bool ArcinfoAsciiGrid::build(istream & ins)
{
    //read header
    string keyword;
    for (int i = 0; i < 6; i++) //loop 6 times at most
    {
        if (isdigit(ins.peek())) break;

        ins >> keyword;

        keyword = tolower(keyword);

        if (keyword == "ncols") ins >> m_ncols;
        if (keyword == "nrows") ins >> m_nrows;
        if (keyword == "xllcorner") ins >> m_xllcorner;
        if (keyword == "yllcorner") ins >> m_yllcorner;
        if (keyword == "cellsize") ins >> m_cellsize;
        if (keyword == "nodata_value") ins >> m_nodata_value;
    }

    if (m_cellsize <= 0 || m_ncols <= 0 || m_nrows <= 0) return false;

    //allocate mem
    unsigned int size = m_nrows * m_ncols;
    m_values = new double[size];
    assert(m_values);

    //read values
    for (unsigned int i = 0; i < size; i++)
    {
        ins >> m_values[i];
    }

    //determine yulcorner
    //yulcorner=yllcorner-nrows*cellsize;

    return true;
}

//perform bilinear interpolation between (i,j), (i+1,j), (i,j+1), (i+1,j+1)
//x and y must be between 0 and 1
double bilinear_interpolation(AsciiGrid * grid, unsigned int i, unsigned int j, double x, double y)
{
    //get values
    double ul = grid->get_value(i, j);
    double ur = grid->get_value(i, j + 1);
    double lr = grid->get_value(i + 1, j + 1);
    double ll = grid->get_value(i + 1, j);

    //check the validness of the values...
    if (grid->is_valid(ul) == false && grid->is_valid(ur) == false && grid->is_valid(ll) == false &&
        grid->is_valid(lr) == false)
        return ul; //All values are invalid

    if (grid->is_valid(ul) && grid->is_valid(ur) == false)
    {
        ur = ul;
    }
    if (grid->is_valid(ul) == false && grid->is_valid(ur))
    {
        ul = ur;
    }
    if (grid->is_valid(lr) && grid->is_valid(ll) == false)
    {
        ll = lr;
    }
    if (grid->is_valid(lr) == false && grid->is_valid(ll))
    {
        lr = ll;
    }

    //interpolate in x dir
    double x1, x2;

    if (grid->is_valid(ul) == false && grid->is_valid(ur) == false)
        x1 = ul; //avoid interpolating invalid values...
    else
        x1 = (1 - x) * ul + x * ur;

    if (grid->is_valid(ll) == false && grid->is_valid(lr) == false)
        x2 = ll; //avoid interpolating invalid values...
    else
        x2 = (1 - x) * ll + x * lr;

    if (grid->is_valid(x1) && grid->is_valid(x2) == false)
    {
        x2 = x1;
    }
    if (grid->is_valid(x1) == false && grid->is_valid(x2))
    {
        x1 = x2;
    }

    //interpolate in y dir
    return (1 - y) * x1 + y * x2;
}

} //end of namespace cusg