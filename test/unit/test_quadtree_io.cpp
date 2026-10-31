// Quadtree numbering, reproducibility, rotation, and tree structure files
#include <cstdio>
#include <fstream>
#include <sstream>
#include <string>

#include <gtest/gtest.h>

#include "csr_checks.h"
#include "def_parser.h"
#include "def_struct.h"
#include "grid_helpers.h"

using namespace cusg;
using gridgen_test::expect_csr_equal;
using gridgen_test::expect_csr_valid;
using gridgen_test::make_mfgrid;
using gridgen_test::make_qtree;

namespace
{

// 2 layers of 3 x 3 cells of 100, refined twice at the center, smoothed
QuadTree3D * make_refined(double angle = 0.0)
{
    return make_qtree(make_mfgrid(2, 3, 3, 100.0, 100.0, 10.0, 0.0, 0.0, angle),
                      {{150.0, 150.0, 2, 0}}, true);
}

// write qtree to a tree structure file, read it back through a QUADTREE
// block, and return the rebuilt quadtree
QuadTree3D * round_trip(QuadTree3D * qtree, const std::string & name)
{
    std::string tsf = testing::TempDir() + name;
    {
        std::ofstream out(tsf.c_str());
        out << *qtree;
    }

    std::ostringstream dfn;
    dfn << "BEGIN MODFLOW_GRID mfg\n"
        << "  ROTATION_ANGLE = 0\n  X_OFFSET = 0\n  Y_OFFSET = 0\n"
        << "  LENGTH_UNIT = undefined\n  NLAY = 2\n  NROW = 3\n  NCOL = 3\n"
        << "  DELR = CONSTANT 100.\n  DELC = CONSTANT 100.\n"
        << "  TOP = CONSTANT 0.\n"
        << "  BOTTOM LAYER 1 = CONSTANT -10.\n"
        << "  BOTTOM LAYER 2 = CONSTANT -20.\n"
        << "END MODFLOW_GRID\n"
        << "BEGIN QUADTREE qtg\n"
        << "  MODFLOW_GRID = mfg\n"
        << "  STRUCTURE_FILE = OPEN/CLOSE " << tsf << "\n";
    for (int k = 0; k < 2; k++)
    {
        dfn << "  TOP LAYER " << k + 1 << " = CONSTANT " << -10.0 * k << "\n"
            << "  BOTTOM LAYER " << k + 1 << " = CONSTANT " << -10.0 * (k + 1) << "\n";
    }
    dfn << "END QUADTREE\n";

    std::istringstream in(dfn.str());
    DefinitionParser parser;
    QuadTree3D * rebuilt = nullptr;
    if (parser.parse(in))
    {
        quadtree_grid_raw_data * data =
            dynamic_cast<quadtree_grid_raw_data *>(parser.getBlockData("qtg"));
        if (data != nullptr && data->initialize(parser.getBlockData()))
        {
            rebuilt = dynamic_cast<QuadTree3D *>(data->build());
        }
    }
    std::remove(tsf.c_str());
    return rebuilt;
}

} // namespace

TEST(QuadTreeIO, BuildTwiceIsIdentical)
{
    CSRData a = make_refined()->get_usg_csr_data();
    CSRData b = make_refined()->get_usg_csr_data();
    expect_csr_equal(a, b);
    a.destory();
    b.destory();
}

TEST(QuadTreeIO, RotationDoesNotChangeConnectivity)
{
    CSRData a = make_refined()->get_usg_csr_data();
    CSRData b = make_refined(0.5236)->get_usg_csr_data();
    expect_csr_equal(a, b);
    a.destory();
    b.destory();
}

TEST(QuadTreeIO, NodesNumberedLayerByLayer)
{
    QuadTree3D * qtree = make_refined();
    int nodes = qtree->nodes();

    int * nodelay = qtree->get_nodelay();
    EXPECT_EQ(nodelay[0] + nodelay[1], nodes);

    for (int n = 0; n < nodes; n++)
    {
        int expected = (n < nodelay[0]) ? 0 : 1;
        EXPECT_EQ(qtree->get_nodeobj(n)->layer, expected) << "node " << n;
    }
}

TEST(QuadTreeIO, OneBasedNumberingOffsetsIaAndJa)
{
    QuadTree3D * qtree = make_refined();
    CSRData zero = qtree->get_usg_csr_data();
    expect_csr_valid(zero);

    qtree->set_one_based_node_numbering(true);
    CSRData one = qtree->get_usg_csr_data();

    ASSERT_EQ(zero.nja, one.nja);
    for (int n = 0; n <= zero.nodes; n++) EXPECT_EQ(one.ia[n], zero.ia[n] + 1);
    for (int i = 0; i < zero.nja; i++) EXPECT_EQ(one.ja[i], zero.ja[i] + 1);
    zero.destory();
    one.destory();
}

TEST(QuadTreeIO, TreeStructureFileRoundTrip)
{
    QuadTree3D * qtree = make_refined();
    qtree->set_one_based_node_numbering(true);
    QuadTree3D * rebuilt = round_trip(qtree, "gridgen_unit_roundtrip.tsf");
    ASSERT_NE(rebuilt, nullptr);
    EXPECT_TRUE(rebuilt->get_one_based_node_numbering());

    CSRData a = qtree->get_usg_csr_data();
    CSRData b = rebuilt->get_usg_csr_data();
    expect_csr_equal(a, b);
    a.destory();
    b.destory();
}

TEST(QuadTreeIO, TreeStructureFileRoundTripZeroBased)
{
    // the numbering is read back from the file (#50)
    QuadTree3D * qtree = make_refined();
    QuadTree3D * rebuilt = round_trip(qtree, "gridgen_unit_roundtrip0.tsf");
    ASSERT_NE(rebuilt, nullptr);
    EXPECT_FALSE(rebuilt->get_one_based_node_numbering());

    CSRData a = qtree->get_usg_csr_data();
    CSRData b = rebuilt->get_usg_csr_data();
    expect_csr_equal(a, b);
    a.destory();
    b.destory();
}

TEST(QuadTreeIO, TreeStructureFileIdentifiersAreOneBased)
{
    // the (layer, row, column) identifiers and quadrant numbers follow the
    // report whatever the node numbering; only the node numbers differ (#50)
    QuadTree3D * qtree = make_refined();
    std::ostringstream zero;
    zero << *qtree;
    qtree->set_one_based_node_numbering(true);
    std::ostringstream one;
    one << *qtree;

    std::istringstream zin(zero.str());
    std::istringstream oin(one.str());
    std::string zline, oline;
    std::getline(zin, zline);
    std::getline(oin, oline);
    ASSERT_EQ(zline, oline);

    int first_active = -1;
    while (std::getline(zin, zline) && std::getline(oin, oline))
    {
        SCOPED_TRACE(zline);
        std::size_t zcomma = zline.find(',');
        std::size_t ocomma = oline.find(',');
        ASSERT_NE(zcomma, std::string::npos);
        EXPECT_EQ(zline.substr(zcomma), oline.substr(ocomma));

        int number = std::stoi(zline.substr(0, zcomma));
        if (number >= 0 && first_active < 0) first_active = number;

        int lay = 0, row = 0, col = 0;
        char code[32] = "";
        int nread =
            std::sscanf(zline.c_str() + zcomma, ", (%d,%d,%d) %31s", &lay, &row, &col, code);
        ASSERT_GE(nread, 3);
        EXPECT_GE(lay, 1);
        EXPECT_GE(row, 1);
        EXPECT_GE(col, 1);
        for (const char * c = code; *c; c++)
        {
            EXPECT_TRUE(*c >= '1' && *c <= '4') << "quadrant " << *c;
        }
    }
    EXPECT_EQ(first_active, 0);
}
