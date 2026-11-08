// Active domains read from shapefiles
#include <string>

#include <gtest/gtest.h>

#include "GridDomain.h"
#include "def_struct.h"
#include "grid_helpers.h"

using namespace cusg;
using gridgen_test::make_mfgrid;
using gridgen_test::make_qtree;

namespace
{

// active cells in layer, and all cells in it
void count_cells(QuadTree3D * qtree, int layer, int & active, int & total)
{
    active = total = 0;
    for (Box * b : qtree->node_group)
    {
        if (b == nullptr || !b->is_leaf || b->layer != layer) continue;
        total++;
        if (b->active) active++;
    }
}

} // namespace

TEST(GridDomain, PolygonEnclosingGridActivatesEveryCell)
{
    // test/shapefiles/enclosing_domain is a rectangle 5 km beyond grid01mfg;
    // applied to layer 2 of a grid with its dimensions, every cell in layer 2
    // is active and layer 1 is unchanged (#53)
    std::unique_ptr<QuadTree3D> qtree =
        make_qtree(make_mfgrid(2, 100, 58, 800.0, 800.0, 10.0, 543750.0, 2792250.0),
                   {{567000.0, 2830000.0, 2, 1}}, true);

    DefinitionBlock block;
    ActiveDomainBlock domain(block);
    domain.shapefile = std::string(GRIDGEN_TEST_DIR) + "/shapefiles/enclosing_domain";
    domain.feature_type = "polygon";

    qtree->node_group.disable_nodes_in_layer(1);
    int active = 0, total = 0;
    count_cells(qtree.get(), 1, active, total);
    ASSERT_EQ(active, 0);

    GridDomain().build(qtree.get(), 1, &domain);
    count_cells(qtree.get(), 1, active, total);
    EXPECT_GT(total, 5800);
    EXPECT_EQ(active, total);

    count_cells(qtree.get(), 0, active, total);
    EXPECT_EQ(active, total);
}
