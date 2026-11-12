// The MODFLOW 6 IHC, HWVA, and ANGLDEGX arrays and the MODFLOW-USG IVC array
#include <cmath>
#include <cstdlib>
#include <map>
#include <memory>
#include <set>
#include <sstream>
#include <utility>
#include <vector>

#include <gtest/gtest.h>

#include "GridVertices.h"
#include "def_parser.h"
#include "def_struct.h"
#include "grid_helpers.h"

using namespace cusg;
using gridgen_test::make_mfgrid;
using gridgen_test::make_qtree;

namespace
{

const double DEGREES_TO_RADIANS = 3.14159265358979323846 / 180.0;

// position of the m-n connection for each n-m connection, zero-based ia and ja
std::vector<int> reverse_positions(const CSRData & csr)
{
    std::map<std::pair<int, int>, int> position;
    for (int n = 0; n < csr.nodes; n++)
        for (int ipos = csr.ia[n]; ipos < csr.ia[n + 1]; ipos++)
            position[std::make_pair(n, csr.ja[ipos])] = ipos;
    std::vector<int> reverse(csr.nja);
    for (const auto & entry : position)
        reverse[entry.second] = position.at(std::make_pair(entry.first.second, entry.first.first));
    return reverse;
}

// IHC and IVC follow fldr, and the diagonal is 0 in each array
void expect_connection_types(const CSRData & csr)
{
    std::vector<int> ihc = csr.ihc();
    std::vector<int> ivc = csr.ivc();
    std::vector<double> angldegx = csr.angldegx(0.0);
    for (int ipos = 0; ipos < csr.nja; ipos++)
    {
        int direction = std::abs(csr.fldr[ipos]);
        EXPECT_EQ(ihc[ipos], (direction == 1 || direction == 2) ? 1 : 0) << "ipos " << ipos;
        EXPECT_EQ(ivc[ipos], direction == 3 ? 1 : 0) << "ipos " << ipos;
        if (direction == 0)
        {
            EXPECT_EQ(csr.hwva[ipos], 0.0);
            EXPECT_EQ(angldegx[ipos], 0.0);
        }
    }
}

// a quadtree of 4 x 4 cells of 100 in 2 layers of 10, with its lower left
// corner at (1000, 2000), refined twice in layer 1 and once in layer 2 at the
// point 150, 250 from that corner; points are refined before the rotation
std::unique_ptr<QuadTree3D> make_refined_qtree(double rotation_radians)
{
    return make_qtree(make_mfgrid(2, 4, 4, 100.0, 100.0, 10.0, 1000.0, 2000.0, rotation_radians),
                      {{1150.0, 2250.0, 2, 0}, {1150.0, 2250.0, 1, 1}}, true);
}

} // namespace

TEST(DisuArrays, ModflowGridWidthsAndAreas)
{
    // 2 layers of 2 rows x 3 columns, delr 10, delc 20
    std::unique_ptr<ModflowGrid> mfgrid = make_mfgrid(2, 2, 3, 10.0, 20.0, 7.5);
    CSRData csr = mfgrid->get_usg_csr_data();
    expect_connection_types(csr);
    for (int ipos = 0; ipos < csr.nja; ipos++)
    {
        int direction = std::abs(csr.fldr[ipos]);
        if (direction == 1) EXPECT_EQ(csr.hwva[ipos], 20.0);
        if (direction == 2) EXPECT_EQ(csr.hwva[ipos], 10.0);
        if (direction == 3) EXPECT_EQ(csr.hwva[ipos], 200.0);
    }
}

TEST(DisuArrays, QuadtreeWidthsAreSymmetric)
{
    std::unique_ptr<QuadTree3D> qtree = make_refined_qtree(0.0);
    CSRData csr = qtree->get_usg_csr_data();
    expect_connection_types(csr);
    std::vector<int> reverse = reverse_positions(csr);
    int different_sizes = 0;
    for (int n = 0; n < csr.nodes; n++)
    {
        Box * box = qtree->get_cell_box(n);
        for (int ipos = csr.ia[n] + 1; ipos < csr.ia[n + 1]; ipos++)
            if (std::abs(csr.fldr[ipos]) == 1 && qtree->get_cell_box(csr.ja[ipos])->dy != box->dy)
                different_sizes++;
        for (int ipos = csr.ia[n] + 1; ipos < csr.ia[n + 1]; ipos++)
        {
            Box * neighbor = qtree->get_cell_box(csr.ja[ipos]);
            int direction = std::abs(csr.fldr[ipos]);
            // MODFLOW 6 requires the two values of a connection to be equal
            EXPECT_EQ(csr.hwva[ipos], csr.hwva[reverse[ipos]]) << "n " << n;
            if (direction == 1)
                EXPECT_EQ(csr.hwva[ipos], std::min(box->dy, neighbor->dy)) << "n " << n;
            else if (direction == 2)
                EXPECT_EQ(csr.hwva[ipos], std::min(box->dx, neighbor->dx)) << "n " << n;
            else
                EXPECT_EQ(csr.hwva[ipos], csr.fahl[ipos]) << "n " << n;
            // the layers are 10 thick, so a horizontal face area is 10 times its width
            if (direction != 3) EXPECT_DOUBLE_EQ(csr.fahl[ipos], 10.0 * csr.hwva[ipos]);
        }
    }
    // the refinement gives x connections between cells of different sizes
    EXPECT_GT(different_sizes, 0);
}

TEST(DisuArrays, VerticalPassThroughAreas)
{
    // the cells of layer 2 under the refined cell of layer 1 are inactive, so
    // layer 1 connects to layer 3 there
    std::unique_ptr<ModflowGrid> mfgrid = make_mfgrid(3, 2, 2, 100.0, 100.0);
    std::unique_ptr<QuadTree3D> qtree = make_qtree(mfgrid.get(), {{50.0, 150.0, 1, 0}});
    qtree->node_group[mfgrid->get_box_index(1, 0, 0)]->active = false;
    qtree->number_nodes();
    qtree->refresh();
    qtree->vertical_pass_through = true;

    CSRData csr = qtree->get_usg_csr_data();
    expect_connection_types(csr);
    int pass_through = 0;
    for (int n = 0; n < csr.nodes; n++)
    {
        for (int ipos = csr.ia[n] + 1; ipos < csr.ia[n + 1]; ipos++)
        {
            if (std::abs(csr.fldr[ipos]) != 3) continue;
            EXPECT_EQ(csr.hwva[ipos], csr.fahl[ipos]) << "n " << n;
            int layers =
                std::abs(qtree->get_cell_box(n)->layer - qtree->get_cell_box(csr.ja[ipos])->layer);
            if (layers == 2) pass_through++;
        }
    }
    // 4 cells of layer 1 connected to the cell of layer 3 below them, both ways
    EXPECT_EQ(pass_through, 8);
}

TEST(DisuArrays, AnglesAreTheFaceNormals)
{
    // unrotated, rotated, rotated by a negative angle, and by more than 360
    // degrees; the angle of each face is within 90 degrees of the direction
    // from the cell center to the center of the connected cell
    for (double degrees : {0.0, 17.3, -30.0, 400.0})
    {
        std::unique_ptr<QuadTree3D> qtree = make_refined_qtree(degrees * DEGREES_TO_RADIANS);
        CSRData csr = qtree->get_usg_csr_data();
        // refined, so that cells of different sizes are connected
        EXPECT_GT(csr.nodes, 32) << degrees << " degrees";
        std::vector<double> angldegx = csr.angldegx(degrees);
        std::vector<int> reverse = reverse_positions(csr);
        GridVertices gv(qtree.get());
        std::set<long> angles;
        for (int n = 0; n < csr.nodes; n++)
        {
            const Point2d & center = gv.center(qtree->get_cell_box(n));
            for (int ipos = csr.ia[n] + 1; ipos < csr.ia[n + 1]; ipos++)
            {
                if (std::abs(csr.fldr[ipos]) == 3)
                {
                    EXPECT_EQ(angldegx[ipos], 0.0);
                    continue;
                }
                EXPECT_GE(angldegx[ipos], 0.0);
                EXPECT_LT(angldegx[ipos], 360.0);
                const Point2d & other = gv.center(qtree->get_cell_box(csr.ja[ipos]));
                double toward = std::atan2(other[1] - center[1], other[0] - center[0]);
                double normal = angldegx[ipos] * DEGREES_TO_RADIANS;
                EXPECT_GT(std::cos(normal - toward), 0.0)
                    << degrees << " degrees, n " << n << " m " << csr.ja[ipos];
                double opposite =
                    std::fmod(angldegx[ipos] - angldegx[reverse[ipos]] + 360.0, 360.0);
                EXPECT_NEAR(opposite, 180.0, 1e-9) << degrees << " degrees, n " << n;
                angles.insert(std::lround(angldegx[ipos] * 10.0));
            }
        }
        // the four directions, 90 degrees apart
        double first = std::fmod(std::fmod(degrees, 360.0) + 360.0, 360.0);
        std::set<long> expected;
        for (int quarter = 0; quarter < 4; quarter++)
            expected.insert(std::lround(std::fmod(first + 90.0 * quarter, 360.0) * 10.0));
        EXPECT_EQ(angles, expected) << degrees << " degrees";
    }
}

TEST(DisuArrays, RotationInDegreesIsKeptAsGiven)
{
    // 17.3 degrees does not survive a conversion to radians and back exactly
    std::istringstream dfn("BEGIN MODFLOW_GRID mfg\n"
                           "  ROTATION_ANGLE = 17.3\n"
                           "  X_OFFSET = 0.\n"
                           "  Y_OFFSET = 0.\n"
                           "  NLAY = 1\n"
                           "  NROW = 2\n"
                           "  NCOL = 2\n"
                           "  DELR = CONSTANT 10.\n"
                           "  DELC = CONSTANT 10.\n"
                           "  TOP = CONSTANT 10.\n"
                           "  BOTTOM LAYER 1 = CONSTANT 0.\n"
                           "END MODFLOW_GRID\n");
    DefinitionParser parser;
    ASSERT_TRUE(parser.parse(dfn));
    ModflowGridBlock * data = dynamic_cast<ModflowGridBlock *>(parser.get_block_data("mfg"));
    ASSERT_NE(data, nullptr);
    ASSERT_TRUE(data->initialize(parser.get_block_data()));
    ModflowGrid mfgrid(*data);
    EXPECT_EQ(mfgrid.rotation_degrees, 17.3);

    CSRData csr = mfgrid.get_usg_csr_data();
    std::vector<double> angldegx = csr.angldegx(mfgrid.rotation_degrees);
    for (int ipos = 0; ipos < csr.nja; ipos++)
    {
        if (csr.fldr[ipos] == 1) EXPECT_EQ(angldegx[ipos], 17.3);
        if (csr.fldr[ipos] == -1) EXPECT_DOUBLE_EQ(angldegx[ipos], 197.3);
        if (csr.fldr[ipos] == 2) EXPECT_DOUBLE_EQ(angldegx[ipos], 107.3);
        if (csr.fldr[ipos] == -2) EXPECT_DOUBLE_EQ(angldegx[ipos], 287.3);
    }
}
