// Small grids built in memory for the unit tests
#pragma once

#include <vector>

#include "ModflowGrid.h"
#include "QuadTree3D.h"

namespace gridgen_test
{

using namespace cusg;

// flat layers: top at 0, each layer thk thick
inline ModflowGrid * make_mfgrid(int nlay, int nrow, int ncol, double delr,
                                 double delc, double thk = 10.0,
                                 double xoff = 0.0, double yoff = 0.0,
                                 double angle = 0.0)
{
    int ncpl = nrow * ncol;
    std::vector<double *> surfaces;
    for (int k = 0; k <= nlay; k++)
    {
        double * s = new double[ncpl];
        for (int i = 0; i < ncpl; i++) s[i] = -thk * k;
        surfaces.push_back(s);
    }
    ModflowGrid * mfgrid = new ModflowGrid(nlay, nrow, ncol, delr, delc,
                                           surfaces, xoff, yoff, angle);
    for (double * s : surfaces) delete[] s;
    return mfgrid;
}

struct RefinePoint
{
    double x;
    double y;
    int level;
    int layer;
};

// refine at points, optionally smooth, then number nodes zero based;
// the quadtree modifies mfgrid, so each quadtree needs its own
inline QuadTree3D * make_qtree(ModflowGrid * mfgrid,
                               const std::vector<RefinePoint> & points,
                               bool smooth = false)
{
    QuadTree3D * qtree = new QuadTree3D(mfgrid);
    for (const RefinePoint & rp : points)
    {
        Point2d pt(rp.x, rp.y);
        qtree->refine(pt, rp.level, rp.layer);
    }
    if (smooth) qtree->smooth_refinement(1, 1);
    qtree->set_one_based_node_numbering(false);
    qtree->number_nodes();
    qtree->refresh();
    // allocate the elevation arrays before computing them
    qtree->get_top();
    qtree->get_bot();
    qtree->recompute_z_of_cells();
    return qtree;
}

} // namespace gridgen_test
