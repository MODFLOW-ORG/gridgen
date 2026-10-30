// Elevations interpolated from ESRI ASCII grids
#include <fstream>
#include <string>

#include <gtest/gtest.h>

#include "grid_helpers.h"

using namespace cusg;
using gridgen_test::make_mfgrid;
using gridgen_test::make_qtree;

namespace
{

// a 6 x 6 ASCII grid of 10 m cells from (0, 0), every value 10 except the
// cells from row r0 to r1 and column c0 to c1 (from the top left), which are
// nodata
std::string write_asc(const std::string & name, int r0, int r1, int c0, int c1)
{
    std::string path = testing::TempDir() + name;
    std::ofstream out(path.c_str());
    out << "ncols 6\nnrows 6\nxllcorner 0.0\nyllcorner 0.0\ncellsize 10.0\n"
        << "NODATA_value -9999\n";
    for (int r = 0; r < 6; r++)
    {
        for (int c = 0; c < 6; c++)
        {
            bool nodata = r >= r0 && r <= r1 && c >= c0 && c <= c1;
            out << (nodata ? "-9999" : "10") << (c < 5 ? " " : "\n");
        }
    }
    return path;
}

// a 4 x 4 grid of 15 m cells over the ASCII grid, its bottom read from path
QuadTree3D * qtree_over(const std::string & path)
{
    QuadTree3D * qtree = make_qtree(make_mfgrid(1, 4, 4, 15.0, 15.0), {});
    qtree->setBotOperators(0, QuadTree3D::Z_VALUE_LINEAR_INTERPOLATION_OPEARTOR);
    qtree->setBotSource(0, path);
    return qtree;
}

} // namespace

TEST(AsciiGrid, NodataIsLeftOutOfInterpolation)
{
    // each 15 m cell overlaps some valid 10 m cells, so every bottom is the
    // valid value, both area weighted and bilinear (#54)
    std::string path = write_asc("gridgen_unit_nodata.asc", 2, 2, 2, 2);
    for (bool area_weighted : {true, false})
    {
        SCOPED_TRACE(area_weighted ? "area weighted" : "bilinear");
        QuadTree3D * qtree = qtree_over(path);
        qtree->setBotAreaWeighted(0, area_weighted);
        double * bot = qtree->get_bot();
        for (int n = 0; n < qtree->nodes(); n++) EXPECT_DOUBLE_EQ(bot[n], 10.0) << n;
    }
}

TEST(AsciiGrid, NodataAreaIsLeftOutOfAverage)
{
    // the top left 15 m cell overlaps ASCII cells (row 0, column 0) by 100 m2
    // and (0, 1), (1, 0), and (1, 1) by 50, 50, and 25 m2; with values equal
    // to the column number plus 1 and (1, 1) nodata, the average over the
    // valid area is (100 * 1 + 50 * 2 + 50 * 1) / 200 = 1.25 (#54)
    std::string path = testing::TempDir() + "gridgen_unit_nodata_weights.asc";
    {
        std::ofstream out(path.c_str());
        out << "ncols 6\nnrows 6\nxllcorner 0.0\nyllcorner 0.0\ncellsize 10.0\n"
            << "NODATA_value -9999\n";
        for (int r = 0; r < 6; r++)
        {
            for (int c = 0; c < 6; c++)
            {
                if (r == 1 && c == 1) out << "-9999";
                else out << c + 1;
                out << (c < 5 ? " " : "\n");
            }
        }
    }
    QuadTree3D * qtree = qtree_over(path);
    double * bot = qtree->get_bot();
    Box * b = qtree->get_nodeobj(0);
    ASSERT_DOUBLE_EQ(b->x, 7.5);
    ASSERT_DOUBLE_EQ(b->y, 52.5);
    EXPECT_DOUBLE_EQ(bot[0], 1.25);
}

TEST(AsciiGrid, CellOverOnlyNodataIsAnError)
{
    // the cell from x 15 to 30 and y 30 to 45 lies inside the nodata block
    // from x 10 to 40 and y 20 to 50 (#54)
    std::string path = write_asc("gridgen_unit_nodata_block.asc", 1, 3, 1, 3);
    for (bool area_weighted : {true, false})
    {
        SCOPED_TRACE(area_weighted ? "area weighted" : "bilinear");
        QuadTree3D * qtree = qtree_over(path);
        qtree->setBotAreaWeighted(0, area_weighted);
        EXPECT_EXIT(qtree->get_bot(), testing::ExitedWithCode(1), "");
    }
}

TEST(AsciiGrid, CellBeyondGridTakesEdgeValue)
{
    // an ASCII grid covering x 0 to 30 of the 60 m grid; cells beyond it take
    // the value at its edge rather than the nodata value or an error (#54)
    std::string path = testing::TempDir() + "gridgen_unit_half.asc";
    {
        std::ofstream out(path.c_str());
        out << "ncols 3\nnrows 6\nxllcorner 0.0\nyllcorner 0.0\ncellsize 10.0\n"
            << "NODATA_value -9999\n";
        for (int r = 0; r < 6; r++) out << "10 10 10\n";
    }
    for (bool area_weighted : {true, false})
    {
        SCOPED_TRACE(area_weighted ? "area weighted" : "bilinear");
        QuadTree3D * qtree = qtree_over(path);
        qtree->setBotAreaWeighted(0, area_weighted);
        double * bot = qtree->get_bot();
        for (int n = 0; n < qtree->nodes(); n++) EXPECT_DOUBLE_EQ(bot[n], 10.0) << n;
    }
}

TEST(AsciiGrid, AverageCoversCellThatIsNotSquare)
{
    // ASCII values are the column number plus 1, plus 10 times the row number
    // from the top, so a cell covering whole 10 m ASCII cells averages to
    // x / 10 + 55.5 - y at its center; checked for cells taller than wide and
    // wider than tall (#73)
    std::string path = testing::TempDir() + "gridgen_unit_xy.asc";
    {
        std::ofstream out(path.c_str());
        out << "ncols 6\nnrows 6\nxllcorner 0.0\nyllcorner 0.0\ncellsize 10.0\n"
            << "NODATA_value -9999\n";
        for (int r = 0; r < 6; r++)
        {
            for (int c = 0; c < 6; c++) out << c + 1 + 10 * r << (c < 5 ? " " : "\n");
        }
    }
    struct Shape
    {
        int nrow, ncol;
        double delr, delc;
    };
    for (const Shape & s : {Shape{3, 6, 10.0, 20.0}, Shape{6, 3, 20.0, 10.0}})
    {
        SCOPED_TRACE("delr " + std::to_string(s.delr) + ", delc " + std::to_string(s.delc));
        QuadTree3D * qtree = make_qtree(make_mfgrid(1, s.nrow, s.ncol, s.delr, s.delc), {});
        qtree->setBotOperators(0, QuadTree3D::Z_VALUE_LINEAR_INTERPOLATION_OPEARTOR);
        qtree->setBotSource(0, path);
        double * bot = qtree->get_bot();
        for (int n = 0; n < qtree->nodes(); n++)
        {
            Box * b = qtree->get_nodeobj(n);
            EXPECT_NEAR(bot[n], b->x / 10.0 + 55.5 - b->y, 1e-9) << "cell " << n;
        }
    }
}
