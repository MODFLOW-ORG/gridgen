/*! \file GridVertices.h
\brief one vertex list for the cell corners of a grid
*/
#pragma once

#include <array>
#include <map>
#include <unordered_map>
#include <vector>

#include "Box.h"

namespace cusg
{

class Grid;
class ModflowGrid;

/*!\class GridVertices
 * \brief the cell corners of a grid, each stored once
 *
 * A corner is located from a base grid edge and a binary fraction of the
 * cell beyond it, so a corner shared by adjacent cells is one vertex with one
 * set of coordinates; rotation and offset are applied once per vertex.
 */
class GridVertices
{
public:
    // base grid edge, numerator, and power of two of the fraction beyond it
    typedef array<long long, 3> EdgeKey;

    GridVertices(Grid * grid);

    /// vertex coordinates
    const vector<Point2d> & vertices() const { return m_vertices; }

    /// vertex indices of the corners of a box, clockwise from the upper left
    const int * corners(const Box * box) const;

    /// corner coordinates of a box, clockwise from the upper left
    void get_corners(const Box * box, double * x, double * y) const;

    /// center of a box
    const Point2d & center(const Box * box) const;

private:
    struct Cell
    {
        int corner[4];
        Point2d center;
    };

    void add_box(const Box * box, int col, int srow, long long ix, long long iy, int depth);
    int vertex(const EdgeKey & kx, const EdgeKey & ky);
    double local_x(const EdgeKey & k) const;
    double local_y(const EdgeKey & k) const;
    Point2d to_world(double xl, double yl) const;

    vector<double> m_xe, m_ys;     // local edges from the west and south
    vector<double> m_delr, m_dely; // cell sizes from the west and south
    double m_xoffset, m_yoffset, m_cos, m_sin;

    vector<Point2d> m_vertices;
    map<pair<EdgeKey, EdgeKey>, int> m_index;
    unordered_map<const Box *, Cell> m_cells;
};

} // namespace cusg
