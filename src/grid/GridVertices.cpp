#include "GridVertices.h"

#include <cassert>
#include <cmath>

#include "ModflowGrid.h"
#include "QuadTree3D.h"

namespace cusg
{

namespace
{

// edge, numerator, and power of two in lowest terms, so a point has one key
GridVertices::EdgeKey reduce(long long edge, long long num, long long depth)
{
    while (depth > 0 && num % 2 == 0)
    {
        num /= 2;
        depth--;
    }
    if (num == 0) depth = 0;
    // the far edge of the cell is the near edge of the next one
    if (depth == 0 && num == 1)
    {
        edge++;
        num = 0;
    }
    return {edge, num, depth};
}

} // namespace

GridVertices::GridVertices(Grid * grid)
{
    QuadTree3D * qtree = dynamic_cast<QuadTree3D *>(grid);
    ModflowGrid * mfgrid = (qtree != NULL) ? qtree->getModflowGrid()
                                           : dynamic_cast<ModflowGrid *>(grid);
    assert(mfgrid);

    int nrow = mfgrid->nrow;
    int ncol = mfgrid->ncol;
    m_xe.assign(ncol + 1, 0.0);
    m_delr.assign(mfgrid->delr, mfgrid->delr + ncol);
    for (int c = 0; c < ncol; c++) m_xe[c + 1] = m_xe[c] + m_delr[c];
    m_ys.assign(nrow + 1, 0.0);
    m_dely.assign(nrow, 0.0);
    for (int s = 0; s < nrow; s++)
    {
        m_dely[s] = mfgrid->delc[nrow - 1 - s];
        m_ys[s + 1] = m_ys[s] + m_dely[s];
    }
    m_xoffset = mfgrid->xoffset;
    m_yoffset = mfgrid->yoffset;
    m_cos = cos(mfgrid->rotation);
    m_sin = sin(mfgrid->rotation);

    for (int k = 0; k < mfgrid->nlay; k++)
        for (int r = 0; r < nrow; r++)
            for (int c = 0; c < ncol; c++)
            {
                Box * root = grid->nodegroup[mfgrid->get_nodeid(k, r, c)];
                add_box(root, c, nrow - 1 - r, 0, 0, 0);
            }
}

void GridVertices::add_box(const Box * box, int col, int srow, long long ix,
                           long long iy, int depth)
{
    if (!box->isLeaf)
    {
        // children are NW, NE, SE, SW
        const int east[4] = {0, 1, 1, 0};
        const int north[4] = {1, 1, 0, 0};
        for (int i = 0; i < 4; i++)
            add_box(box->pChildren[i], col, srow, 2 * ix + east[i],
                    2 * iy + north[i], depth + 1);
    }

    EdgeKey west = reduce(col, ix, depth);
    EdgeKey east = reduce(col, ix + 1, depth);
    EdgeKey south = reduce(srow, iy, depth);
    EdgeKey north = reduce(srow, iy + 1, depth);

    Cell cell;
    cell.corner[0] = vertex(west, north);
    cell.corner[1] = vertex(east, north);
    cell.corner[2] = vertex(east, south);
    cell.corner[3] = vertex(west, south);
    cell.center = to_world(local_x(reduce(col, 2 * ix + 1, depth + 1)),
                           local_y(reduce(srow, 2 * iy + 1, depth + 1)));
    m_cells[box] = cell;
}

int GridVertices::vertex(const EdgeKey & kx, const EdgeKey & ky)
{
    pair<EdgeKey, EdgeKey> key(kx, ky);
    map<pair<EdgeKey, EdgeKey>, int>::iterator it = m_index.find(key);
    if (it != m_index.end()) return it->second;
    int id = m_vertices.size();
    m_vertices.push_back(to_world(local_x(kx), local_y(ky)));
    m_index[key] = id;
    return id;
}

double GridVertices::local_x(const EdgeKey & k) const
{
    if (k[1] == 0) return m_xe[k[0]];
    return m_xe[k[0]] + m_delr[k[0]] * ldexp((double)k[1], -(int)k[2]);
}

double GridVertices::local_y(const EdgeKey & k) const
{
    if (k[1] == 0) return m_ys[k[0]];
    return m_ys[k[0]] + m_dely[k[0]] * ldexp((double)k[1], -(int)k[2]);
}

Point2d GridVertices::to_world(double xl, double yl) const
{
    return Point2d(m_xoffset + (xl * m_cos - yl * m_sin),
                   m_yoffset + (xl * m_sin + yl * m_cos));
}

const int * GridVertices::corners(const Box * box) const
{
    return m_cells.at(box).corner;
}

void GridVertices::get_corners(const Box * box, double * x, double * y) const
{
    const int * ids = corners(box);
    for (int i = 0; i < 4; i++)
    {
        x[i] = m_vertices[ids[i]][0];
        y[i] = m_vertices[ids[i]][1];
    }
}

const Point2d & GridVertices::center(const Box * box) const
{
    return m_cells.at(box).center;
}

} // namespace cusg
