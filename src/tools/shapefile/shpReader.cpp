
#include "shpReader.h"
#include "GridgenError.h"
#include "Warnings.h"
#include "intersection.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <cassert>
#include <numeric>
using namespace std;

namespace cusg
{

using uint = unsigned int;

void ShapeReader::destroy()
{
    for (auto & polygon : m_polygon_list) polygon.destroy();
    for (auto & arc : m_arc_list) arc.destroy();
    m_polygon_list.clear();
    m_arc_list.clear();
    m_pt_list.clear();

    m_intFlds.clear();
    m_dbFlds.clear();
    m_strFlds.clear();
    m_logFlds.clear();
}

//
//read field attributes from DBF file according to id
bool ShapeReader::readAttributes(DBFHandle & dbfHandle, int iShape)
{
    for (auto & field : m_intFlds)
    {
        field.values.insert(
            make_pair(iShape, DBFReadIntegerAttribute(dbfHandle, iShape, field.id)));
    }
    for (auto & field : m_dbFlds)
    {
        field.values.insert(make_pair(iShape, DBFReadDoubleAttribute(dbfHandle, iShape, field.id)));
    }
    for (auto & field : m_strFlds)
    {
        field.values.insert(
            make_pair(iShape, string(DBFReadStringAttribute(dbfHandle, iShape, field.id))));
    }
    for (auto & field : m_logFlds)
    {
        field.values.insert(
            make_pair(iShape, string(DBFReadLogicalAttribute(dbfHandle, iShape, field.id))));
    }
    return true;
}

//
// read from files (both shp and db files)
//
bool ShapeReader::read(const string & filename)
{

    if (already_read) return true;

    //start to read
    SHPHandle hSHP;
    hSHP = SHPOpen(filename.c_str(), "rb");

    //the caller reports a shapefile that cannot be opened
    if (hSHP == nullptr) return false;
    m_filename = filename;

    /************************************************************/
    DBFHandle hDBF = DBFOpen(filename.c_str(), "rb");

    //check
    if (hDBF == nullptr)
    {
        SHPClose(hSHP);
        return false;
    }

    ////check more

    //used for reading from db
    int nWidth, nDecimals;

    //position of these fields in database...
    for (int i = 0; i < DBFGetFieldCount(hDBF); i++)
    {
        char szTitle[12];
        DBFFieldType eType = DBFGetFieldInfo(hDBF, i, szTitle, &nWidth, &nDecimals);
        if (eType != FTInteger && eType != FTDouble && eType != FTString && eType != FTLogical)
            continue; //not good...

        //collect the field name and type info
        if (eType == FTInteger)
        {
            ShapeIntField field;
            field.id = i;
            field.name = szTitle;
            m_intFlds.push_back(field);
        }
        else if (eType == FTDouble)
        {
            ShapeDoubleField field;
            field.id = i;
            field.name = szTitle;
            m_dbFlds.push_back(field);
        }
        else if (eType == FTString)
        {
            ShapeStringField field;
            field.id = i;
            field.name = szTitle;
            m_strFlds.push_back(field);
        }
        else if (eType == FTLogical)
        {
            ShapeLogField field;
            field.id = i;
            field.name = szTitle;
            m_logFlds.push_back(field);
        }
    }

    ////check if we can get the field ids for the information we wanted

    /************************************************************/

    //Get info
    int nShapeType, nEntities;
    double adfMinBound[4], adfMaxBound[4];
    SHPGetInfo(hSHP, &nEntities, &nShapeType, adfMinBound, adfMaxBound);

    //reserve space
    m_polygon_list.reserve(nEntities);

    //	Skim over the list of shapes, printing all the vertices
    int nunsupported = 0;
    for (uint i = 0; i < (uint)nEntities; i++)
    {

        //
        SHPObject * psShape = SHPReadObject(hSHP, i);

        if (psShape == nullptr)
        {
            SHPClose(hSHP);
            DBFClose(hDBF);
            throw GridgenError("shape " + to_string(i + 1) + " in shapefile " + filename +
                               " cannot be read");
        }

        //the format requires the first part to start at the first vertex
        if (psShape->nParts > 0 && psShape->panPartStart[0] != 0)
        {
            SHPDestroyObject(psShape);
            SHPClose(hSHP);
            DBFClose(hDBF);
            throw GridgenError("shape " + to_string(i + 1) + " in shapefile " + filename +
                               " cannot be read; its first part does not start at its first "
                               "vertex");
        }

        switch (psShape->nSHPType)
        {
            case SHPT_POLYGON:
            case SHPT_POLYGONZ:
            case SHPT_POLYGONM:
            {
                read_polygon(psShape, i);
                break;
            }
            case SHPT_ARC:
            case SHPT_ARCZ:
            case SHPT_ARCM:
            {
                readArc(psShape, i);
                break;
            }
            case SHPT_POINT:
            case SHPT_POINTM:
            case SHPT_POINTZ:
            case SHPT_MULTIPOINT:
            case SHPT_MULTIPOINTM:
            case SHPT_MULTIPOINTZ:
            {
                readPts(psShape, i);
                break;
            }
            case SHPT_MULTIPATCH:
            {
                readMultiPatch(psShape, i);
                break;
            }
            default:
            {
                nunsupported++;
                break;
            }
        }

        //read the field values
        readAttributes(hDBF, i);

        SHPDestroyObject(psShape);
    }

    SHPClose(hSHP);
    DBFClose(hDBF);

    if (nunsupported > 0)
        warn(to_string(nunsupported) + " of the " + to_string(nEntities) + " shapes in shapefile " +
             filename + " have no geometry or a type that is not supported, and are ignored.");

    //done reading the shp file

    //
    //reading from database now.
    //

    already_read = true;
    return true;
}

//whether p is on segment ab, within a small fraction of its length
static bool onSegment(const Point2d & p, const Point2d & a, const Point2d & b)
{
    double dx = b[0] - a[0], dy = b[1] - a[1];
    double len2 = dx * dx + dy * dy;
    if (len2 == 0) return p[0] == a[0] && p[1] == a[1]; //the closing vertex repeats the first
    double cross = dx * (p[1] - a[1]) - dy * (p[0] - a[0]);
    double dot = dx * (p[0] - a[0]) + dy * (p[1] - a[1]);
    double tol = 1e-9 * len2;
    return fabs(cross) <= tol && dot >= -tol && dot <= len2 + tol;
}

//whether ring a is inside ring b, which it does not cross; the first vertex
//or edge midpoint of a that is not on b decides, since a may touch b
static bool ringInRing(const vector<Point2d> & a, const vector<Point2d> & b)
{
    size_t size_a = a.size(), size_b = b.size();
    for (size_t j = 0; j < 2 * size_a; j++)
    {
        const Point2d & v = a[j / 2];
        Point2d p = (j % 2 == 0) ? v
                                 : Point2d((v[0] + a[(j / 2 + 1) % size_a][0]) / 2,
                                           (v[1] + a[(j / 2 + 1) % size_a][1]) / 2);
        bool on_b = false, inside = false;
        for (size_t k = 0; k < size_b && !on_b; k++)
        {
            const Point2d & b0 = b[k];
            const Point2d & b1 = b[(k + 1) % size_b];
            on_b = onSegment(p, b0, b1);
            if (CrossesRay(b0.get(), b1.get(), p.get())) inside = !inside;
        }
        if (!on_b) return inside;
    }
    return false; //a lies on b
}

using Box = array<double, 4>; //xmin, xmax, ymin, ymax

//whether box outer holds box inner
static bool boxHolds(const Box & outer, const Box & inner)
{
    return outer[0] <= inner[0] && inner[1] <= outer[1] && outer[2] <= inner[2] &&
           inner[3] <= outer[3];
}

//a tree of boxes that finds the boxes holding a given box without testing
//every box; each node bounds a run of the boxes, split in two at the median
//center along its longer side
class BoxTree
{
public:
    explicit BoxTree(const vector<Box> & boxes) : boxes(boxes), order(boxes.size())
    {
        iota(order.begin(), order.end(), 0);
        if (!boxes.empty()) build(0, order.size());
    }

    //the indices of the boxes that hold box query, in ascending order
    vector<size_t> holding(const Box & query) const
    {
        vector<size_t> found, stack;
        if (!nodes.empty()) stack.push_back(0);
        while (!stack.empty())
        {
            const Node & node = nodes[stack.back()];
            stack.pop_back();
            //no box in a node can hold a box its bounds do not
            if (!boxHolds(node.bounds, query)) continue;
            if (node.left < 0)
            {
                for (size_t j = node.begin; j < node.end; j++)
                    if (boxHolds(boxes[order[j]], query)) found.push_back(order[j]);
            }
            else
            {
                stack.push_back(node.left);
                stack.push_back(node.right);
            }
        }
        sort(found.begin(), found.end());
        return found;
    }

private:
    struct Node
    {
        Box bounds;
        size_t begin, end; //the run of order the node bounds
        int left, right;   //the child nodes, or -1 for a leaf
    };

    static const size_t leaf_size = 8;
    const vector<Box> & boxes;
    vector<size_t> order;
    vector<Node> nodes;

    //add the node bounding the boxes order[begin, end) and the nodes below it
    int build(size_t begin, size_t end)
    {
        //fmin and fmax skip a coordinate that is not a number, so such a box
        //does not hide the boxes beside it
        Box bounds = {HUGE_VAL, -HUGE_VAL, HUGE_VAL, -HUGE_VAL};
        for (size_t j = begin; j < end; j++)
        {
            const Box & box = boxes[order[j]];
            bounds = {fmin(bounds[0], box[0]), fmax(bounds[1], box[1]), fmin(bounds[2], box[2]),
                      fmax(bounds[3], box[3])};
        }
        int index = static_cast<int>(nodes.size());
        nodes.push_back({bounds, begin, end, -1, -1});
        if (end - begin <= leaf_size) return index;

        int axis = bounds[1] - bounds[0] >= bounds[3] - bounds[2] ? 0 : 2;
        size_t middle = begin + (end - begin) / 2;
        //a coordinate that is not a number sorts first, so the order is strict
        auto center = [&](size_t i)
        {
            double sum = boxes[i][axis] + boxes[i][axis + 1];
            return isnan(sum) ? -HUGE_VAL : sum;
        };
        nth_element(order.begin() + begin, order.begin() + middle, order.begin() + end,
                    [&](size_t a, size_t b) { return center(a) < center(b); });
        int left = build(begin, middle);
        int right = build(middle, end);
        nodes[index].left = left;
        nodes[index].right = right;
        return index;
    }
};

void ShapeReader::read_polygon(SHPObject * psShape, int iShape)
{
    //the vertices of each part
    vector<vector<Point2d>> rings;
    int iPart = 0;
    for (int j = 0; j < psShape->nVertices; j++)
    {
        //a part with no vertices starts at the same vertex as the next part
        bool start = rings.empty();
        while (iPart < psShape->nParts && psShape->panPartStart[iPart] <= j)
        {
            start = start || psShape->panPartStart[iPart] == j;
            iPart++;
        }
        if (start) rings.emplace_back();

        //Points may duplicate....
        double & x = psShape->padfX[j];
        double & y = psShape->padfY[j];
        rings.back().emplace_back(x, y);

        if (x < bbox[0]) bbox[0] = x;
        if (x > bbox[1]) bbox[1] = x;
        if (y < bbox[2]) bbox[2] = y;
        if (y > bbox[3]) bbox[3] = y;
    }

    //a ring inside an odd number of the others is a hole of the smallest ring
    //containing it, and any other ring is an outer boundary, whatever the
    //order and orientation of the rings
    size_t nrings = rings.size();
    vector<int> depth(nrings, 0), parent(nrings, -1);
    vector<Box> box(nrings);
    for (size_t i = 0; i < nrings; i++)
    {
        box[i] = {rings[i][0][0], rings[i][0][0], rings[i][0][1], rings[i][0][1]};
        for (const Point2d & p : rings[i])
        {
            box[i][0] = min(box[i][0], p[0]);
            box[i][1] = max(box[i][1], p[0]);
            box[i][2] = min(box[i][2], p[1]);
            box[i][3] = max(box[i][3], p[1]);
        }
    }
    auto extent = [&](size_t i) { return (box[i][1] - box[i][0]) * (box[i][3] - box[i][2]); };
    BoxTree tree(box);
    for (size_t i = 0; i < nrings; i++)
    {
        //ring k can contain ring i only if its box holds the box of ring i
        for (size_t k : tree.holding(box[i]))
        {
            if (k == i || !ringInRing(rings[i], rings[k])) continue;
            depth[i]++;
            if (parent[i] < 0 || extent(k) < extent(parent[i])) parent[i] = static_cast<int>(k);
        }
    }

    //outer rings first, so each hole finds the polygon of its ring; a hole
    //whose ring is not an outer ring, which happens only when rings cross, is
    //an outer ring
    vector<int> polygon_of(nrings, -1);
    for (int hole = 0; hole <= 1; hole++)
    {
        for (size_t i = 0; i < nrings; i++)
        {
            if (depth[i] % 2 != hole) continue;
            int owner = hole && parent[i] >= 0 ? polygon_of[parent[i]] : -1;
            Ring ring(owner >= 0 ? Ring::PIN : Ring::POUT);
            ring.beginPoly();
            for (const Point2d & p : rings[i]) ring.addVertex(p[0], p[1]);
            ring.endPoly();
            if (owner >= 0)
                m_polygon_list[owner].push_back(ring);
            else
            {
                polygon_of[i] = static_cast<int>(m_polygon_list.size());
                add_polygon(ring, iShape);
            }
        }
    }
}
void ShapeReader::readMultiPatch(SHPObject * psShape, int iShape)
{
    //outer and first rings start a polygon, and inner rings are its holes; a
    //plain ring is a hole after another ring
    bool started = false;
    for (int i = 0; i < psShape->nParts; i++)
    {
        int stId = psShape->panPartStart[i];
        int endId = psShape->nVertices - 1;
        if (i < psShape->nParts - 1)
        {
            endId = psShape->panPartStart[i + 1] - 1;
        }
        int type = psShape->panPartType[i];
        if (type != SHPP_OUTERRING && type != SHPP_INNERRING && type != SHPP_FIRSTRING &&
            type != SHPP_RING)
        {
            warn("Part " + to_string(i + 1) + " of multipatch shape " + to_string(iShape + 1) +
                 " in shapefile " + m_filename + " is a " + SHPPartTypeName(type) +
                 ", which is not supported, and is skipped.");
            continue;
        }

        if (type == SHPP_INNERRING && !started)
        {
            warn("Part " + to_string(i + 1) + " of multipatch shape " + to_string(iShape + 1) +
                 " in shapefile " + m_filename +
                 " is an inner ring that follows no outer ring, and is skipped.");
            continue;
        }

        bool outer = !started || type == SHPP_OUTERRING || type == SHPP_FIRSTRING;
        Ring ring(outer ? Ring::POUT : Ring::PIN);
        ring.beginPoly();
        for (int j = stId; j <= endId; j++)
        {
            //Points may duplicate....
            double & x = psShape->padfX[j];
            double & y = psShape->padfY[j];
            ring.addVertex(x, y);

            if (x < bbox[0]) bbox[0] = x;
            if (x > bbox[1]) bbox[1] = x;
            if (y < bbox[2]) bbox[2] = y;
            if (y > bbox[3]) bbox[3] = y;
        }
        ring.endPoly();
        add_polygon(ring, iShape);
        started = true;
    }
}

void ShapeReader::readArc(SHPObject * psShape, int iShape)
{
    //read each part
    int iPart = 1;

    //
    GIS_polyline arc;
    arc.beginPoly();

    for (int j = 0; j < psShape->nVertices; j++)
    {
        string pszPartType;

        if (j == 0 && psShape->nParts > 0) pszPartType = SHPPartTypeName(psShape->panPartType[0]);

        //this defines a new loop...
        if (iPart < psShape->nParts && psShape->panPartStart[iPart] == j)
        {
            pszPartType = SHPPartTypeName(psShape->panPartType[iPart]);
            iPart++;
            arc.endPoly();
            arc.m_id = iShape;
            arc.computeGeodesicDist();
            m_arc_list.push_back(arc);

            //start a new polyline
            arc = GIS_polyline();
            arc.beginPoly();
        }

        //Points may duplicate....
        double & x = psShape->padfX[j];
        double & y = psShape->padfY[j];
        arc.addVertex(x, y);

        if (x < bbox[0]) bbox[0] = x;
        if (x > bbox[1]) bbox[1] = x;
        if (y < bbox[2]) bbox[2] = y;
        if (y > bbox[3]) bbox[3] = y;
    }

    arc.endPoly();
    arc.m_id = iShape;
    arc.computeGeodesicDist();
    m_arc_list.push_back(arc);
}

void ShapeReader::readPts(SHPObject * psShape, int iShape)
{
    m_pt_list.reserve(psShape->nVertices);
    for (int j = 0; j < psShape->nVertices; j++)
    {
        //Points may duplicate....
        double & x = psShape->padfX[j];
        double & y = psShape->padfY[j];
        if (x < bbox[0]) bbox[0] = x;
        if (x > bbox[1]) bbox[1] = x;
        if (y < bbox[2]) bbox[2] = y;
        if (y > bbox[3]) bbox[3] = y;

        GIS_Point2d pt(x, y);
        pt.m_id = iShape;
        m_pt_list.push_back(pt);
    }
}

void ShapeReader::add_polygon(Ring & ring, int iShape)
{
    if (ring.getType() == Ring::POUT)
    {
        //create a new polygon
        GIS_polygon polygon;
        polygon.m_id = iShape;
        polygon.push_back(ring);
        m_polygon_list.push_back(polygon);
    }
    else
    { //Ring::PIN
        //add to the last polygon
        m_polygon_list.back().push_back(ring);
    }
}

} //end namespace cusg
