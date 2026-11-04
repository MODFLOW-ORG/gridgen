// Small grids built in memory for the unit tests
#pragma once

#include <memory>
#include <vector>

#include "ModflowGrid.h"
#include "QuadTree3D.h"

namespace gridgen_test
{

using namespace cusg;

// flat layers: top at 0, each layer thk thick
inline std::unique_ptr<ModflowGrid> make_mfgrid(int nlay, int nrow, int ncol, double delr,
                                                double delc, double thk = 10.0, double xoff = 0.0,
                                                double yoff = 0.0, double angle = 0.0)
{
    int ncpl = nrow * ncol;
    std::vector<std::vector<double>> layers(nlay + 1, std::vector<double>(ncpl));
    std::vector<double *> surfaces;
    for (int k = 0; k <= nlay; k++)
    {
        for (int i = 0; i < ncpl; i++) layers[k][i] = -thk * k;
        surfaces.push_back(layers[k].data());
    }
    return std::make_unique<ModflowGrid>(nlay, nrow, ncol, delr, delc, surfaces, xoff, yoff, angle);
}

// a quadtree that owns the MODFLOW grid it refines
class OwningQuadTree : public QuadTree3D
{
public:
    explicit OwningQuadTree(std::unique_ptr<ModflowGrid> mfgrid)
        : QuadTree3D(mfgrid.get()), m_owned_mfgrid(std::move(mfgrid))
    {
    }

private:
    std::unique_ptr<ModflowGrid> m_owned_mfgrid;
};

struct RefinePoint
{
    double x;
    double y;
    int level;
    int layer;
};

// refine at points, optionally smooth, then number nodes zero based
inline void refine_qtree(QuadTree3D * qtree, const std::vector<RefinePoint> & points, bool smooth)
{
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
}

// a refined quadtree that owns mfgrid; the quadtree modifies mfgrid, so each
// quadtree needs its own
inline std::unique_ptr<QuadTree3D> make_qtree(std::unique_ptr<ModflowGrid> mfgrid,
                                              const std::vector<RefinePoint> & points,
                                              bool smooth = false)
{
    std::unique_ptr<QuadTree3D> qtree = std::make_unique<OwningQuadTree>(std::move(mfgrid));
    refine_qtree(qtree.get(), points, smooth);
    return qtree;
}

// a refined quadtree of mfgrid, which must outlive it
inline std::unique_ptr<QuadTree3D>
make_qtree(ModflowGrid * mfgrid, const std::vector<RefinePoint> & points, bool smooth = false)
{
    std::unique_ptr<QuadTree3D> qtree = std::make_unique<QuadTree3D>(mfgrid);
    refine_qtree(qtree.get(), points, smooth);
    return qtree;
}

} // namespace gridgen_test
