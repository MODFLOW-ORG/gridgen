
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

void ShapeReader::destroy()
{
    for (auto & polygon : m_polygon_list) polygon.destroy();
    for (auto & arc : m_arc_list) arc.destroy();
    m_polygon_list.clear();
    m_arc_list.clear();
    m_pt_list.clear();

    m_int_fields.clear();
    m_double_fields.clear();
    m_str_fields.clear();
    m_logical_fields.clear();
}

//
//read field attributes from DBF file according to id
bool ShapeReader::read_attributes(DBFHandle & dbf_handle, int shape_index)
{
    for (auto & field : m_int_fields)
    {
        field.values.insert(
            make_pair(shape_index, DBFReadIntegerAttribute(dbf_handle, shape_index, field.id)));
    }
    for (auto & field : m_double_fields)
    {
        field.values.insert(
            make_pair(shape_index, DBFReadDoubleAttribute(dbf_handle, shape_index, field.id)));
    }
    for (auto & field : m_str_fields)
    {
        field.values.insert(make_pair(
            shape_index, string(DBFReadStringAttribute(dbf_handle, shape_index, field.id))));
    }
    for (auto & field : m_logical_fields)
    {
        field.values.insert(make_pair(
            shape_index, string(DBFReadLogicalAttribute(dbf_handle, shape_index, field.id))));
    }
    return true;
}

//
// read from files (both shp and db files)
//
bool ShapeReader::read(const string & filename)
{

    if (m_already_read) return true;

    //start to read
    SHPHandle shp_handle;
    shp_handle = SHPOpen(filename.c_str(), "rb");

    //the caller reports a shapefile that cannot be opened
    if (shp_handle == nullptr) return false;
    m_filename = filename;

    /************************************************************/
    DBFHandle dbf_handle = DBFOpen(filename.c_str(), "rb");

    //check
    if (dbf_handle == nullptr)
    {
        SHPClose(shp_handle);
        return false;
    }

    ////check more

    //used for reading from db
    int width, decimals;

    //position of these fields in database...
    for (int i = 0; i < DBFGetFieldCount(dbf_handle); i++)
    {
        char title[12];
        DBFFieldType field_type = DBFGetFieldInfo(dbf_handle, i, title, &width, &decimals);
        if (field_type != FTInteger && field_type != FTDouble && field_type != FTString &&
            field_type != FTLogical)
            continue; //not good...

        //collect the field name and type info
        if (field_type == FTInteger)
        {
            ShapeIntField field;
            field.id = i;
            field.name = title;
            m_int_fields.push_back(field);
        }
        else if (field_type == FTDouble)
        {
            ShapeDoubleField field;
            field.id = i;
            field.name = title;
            m_double_fields.push_back(field);
        }
        else if (field_type == FTString)
        {
            ShapeStringField field;
            field.id = i;
            field.name = title;
            m_str_fields.push_back(field);
        }
        else if (field_type == FTLogical)
        {
            ShapeLogField field;
            field.id = i;
            field.name = title;
            m_logical_fields.push_back(field);
        }
    }

    ////check if we can get the field ids for the information we wanted

    /************************************************************/

    //Get info
    int shape_type, entity_count;
    double min_bound[4], max_bound[4];
    SHPGetInfo(shp_handle, &entity_count, &shape_type, min_bound, max_bound);

    //reserve space
    m_polygon_list.reserve(entity_count);

    //	Skim over the list of shapes, printing all the vertices
    int nunsupported = 0;
    for (unsigned int i = 0; i < (unsigned int)entity_count; i++)
    {

        //
        SHPObject * shape = SHPReadObject(shp_handle, i);

        if (shape == nullptr)
        {
            SHPClose(shp_handle);
            DBFClose(dbf_handle);
            throw GridgenError("shape " + to_string(i + 1) + " in shapefile " + filename +
                               " cannot be read");
        }

        //the format requires the first part to start at the first vertex
        if (shape->nParts > 0 && shape->panPartStart[0] != 0)
        {
            SHPDestroyObject(shape);
            SHPClose(shp_handle);
            DBFClose(dbf_handle);
            throw GridgenError("shape " + to_string(i + 1) + " in shapefile " + filename +
                               " cannot be read; its first part does not start at its first "
                               "vertex");
        }

        switch (shape->nSHPType)
        {
            case SHPT_POLYGON:
            case SHPT_POLYGONZ:
            case SHPT_POLYGONM:
            {
                read_polygon(shape, i);
                break;
            }
            case SHPT_ARC:
            case SHPT_ARCZ:
            case SHPT_ARCM:
            {
                read_arc(shape, i);
                break;
            }
            case SHPT_POINT:
            case SHPT_POINTM:
            case SHPT_POINTZ:
            case SHPT_MULTIPOINT:
            case SHPT_MULTIPOINTM:
            case SHPT_MULTIPOINTZ:
            {
                read_pts(shape, i);
                break;
            }
            case SHPT_MULTIPATCH:
            {
                read_multi_patch(shape, i);
                break;
            }
            default:
            {
                nunsupported++;
                break;
            }
        }

        //read the field values
        read_attributes(dbf_handle, i);

        SHPDestroyObject(shape);
    }

    SHPClose(shp_handle);
    DBFClose(dbf_handle);

    if (nunsupported > 0)
        warn(to_string(nunsupported) + " of the " + to_string(entity_count) +
             " shapes in shapefile " + filename +
             " have no geometry or a type that is not supported, and are ignored.");

    //done reading the shp file

    //
    //reading from database now.
    //

    m_already_read = true;
    return true;
}

//whether p is on segment ab, within a small fraction of its length
static bool on_segment(const Point2d & p, const Point2d & a, const Point2d & b)
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
static bool ring_in_ring(const vector<Point2d> & a, const vector<Point2d> & b)
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
            on_b = on_segment(p, b0, b1);
            if (crosses_ray(b0.get(), b1.get(), p.get())) inside = !inside;
        }
        if (!on_b) return inside;
    }
    return false; //a lies on b
}

using Box = array<double, 4>; //xmin, xmax, ymin, ymax

//whether box outer holds box inner
static bool box_holds(const Box & outer, const Box & inner)
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
    explicit BoxTree(const vector<Box> & boxes) : m_boxes(boxes), m_order(boxes.size())
    {
        iota(m_order.begin(), m_order.end(), 0);
        if (!boxes.empty()) build(0, m_order.size());
    }

    //the indices of the boxes that hold box query, in ascending order
    vector<size_t> holding(const Box & query) const
    {
        vector<size_t> found, stack;
        if (!m_nodes.empty()) stack.push_back(0);
        while (!stack.empty())
        {
            const Node & node = m_nodes[stack.back()];
            stack.pop_back();
            //no box in a node can hold a box its bounds do not
            if (!box_holds(node.bounds, query)) continue;
            if (node.left < 0)
            {
                for (size_t j = node.begin; j < node.end; j++)
                    if (box_holds(m_boxes[m_order[j]], query)) found.push_back(m_order[j]);
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

    static const size_t LEAF_SIZE = 8;
    const vector<Box> & m_boxes;
    vector<size_t> m_order;
    vector<Node> m_nodes;

    //add the node bounding the boxes order[begin, end) and the nodes below it
    int build(size_t begin, size_t end)
    {
        //fmin and fmax skip a coordinate that is not a number, so such a box
        //does not hide the boxes beside it
        Box bounds = {HUGE_VAL, -HUGE_VAL, HUGE_VAL, -HUGE_VAL};
        for (size_t j = begin; j < end; j++)
        {
            const Box & box = m_boxes[m_order[j]];
            bounds = {fmin(bounds[0], box[0]), fmax(bounds[1], box[1]), fmin(bounds[2], box[2]),
                      fmax(bounds[3], box[3])};
        }
        int index = static_cast<int>(m_nodes.size());
        m_nodes.push_back({bounds, begin, end, -1, -1});
        if (end - begin <= LEAF_SIZE) return index;

        int axis = bounds[1] - bounds[0] >= bounds[3] - bounds[2] ? 0 : 2;
        size_t middle = begin + (end - begin) / 2;
        //a coordinate that is not a number sorts first, so the order is strict
        auto center = [&](size_t i)
        {
            double sum = m_boxes[i][axis] + m_boxes[i][axis + 1];
            return isnan(sum) ? -HUGE_VAL : sum;
        };
        nth_element(m_order.begin() + begin, m_order.begin() + middle, m_order.begin() + end,
                    [&](size_t a, size_t b) { return center(a) < center(b); });
        int left = build(begin, middle);
        int right = build(middle, end);
        m_nodes[index].left = left;
        m_nodes[index].right = right;
        return index;
    }
};

void ShapeReader::read_polygon(SHPObject * shape, int shape_index)
{
    //the vertices of each part
    vector<vector<Point2d>> rings;
    int part_index = 0;
    for (int j = 0; j < shape->nVertices; j++)
    {
        //a part with no vertices starts at the same vertex as the next part
        bool start = rings.empty();
        while (part_index < shape->nParts && shape->panPartStart[part_index] <= j)
        {
            start = start || shape->panPartStart[part_index] == j;
            part_index++;
        }
        if (start) rings.emplace_back();

        //Points may duplicate....
        double & x = shape->padfX[j];
        double & y = shape->padfY[j];
        rings.back().emplace_back(x, y);

        if (x < m_bbox[0]) m_bbox[0] = x;
        if (x > m_bbox[1]) m_bbox[1] = x;
        if (y < m_bbox[2]) m_bbox[2] = y;
        if (y > m_bbox[3]) m_bbox[3] = y;
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
            if (k == i || !ring_in_ring(rings[i], rings[k])) continue;
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
            ring.begin_poly();
            for (const Point2d & p : rings[i]) ring.add_vertex(p[0], p[1]);
            ring.end_poly();
            if (owner >= 0)
                m_polygon_list[owner].push_back(ring);
            else
            {
                polygon_of[i] = static_cast<int>(m_polygon_list.size());
                add_polygon(ring, shape_index);
            }
        }
    }
}
void ShapeReader::read_multi_patch(SHPObject * shape, int shape_index)
{
    //outer and first rings start a polygon, and inner rings are its holes; a
    //plain ring is a hole after another ring
    bool started = false;
    for (int i = 0; i < shape->nParts; i++)
    {
        int st_id = shape->panPartStart[i];
        int end_id = shape->nVertices - 1;
        if (i < shape->nParts - 1)
        {
            end_id = shape->panPartStart[i + 1] - 1;
        }
        int type = shape->panPartType[i];
        if (type != SHPP_OUTERRING && type != SHPP_INNERRING && type != SHPP_FIRSTRING &&
            type != SHPP_RING)
        {
            warn("Part " + to_string(i + 1) + " of multipatch shape " + to_string(shape_index + 1) +
                 " in shapefile " + m_filename + " is a " + SHPPartTypeName(type) +
                 ", which is not supported, and is skipped.");
            continue;
        }

        if (type == SHPP_INNERRING && !started)
        {
            warn("Part " + to_string(i + 1) + " of multipatch shape " + to_string(shape_index + 1) +
                 " in shapefile " + m_filename +
                 " is an inner ring that follows no outer ring, and is skipped.");
            continue;
        }

        bool outer = !started || type == SHPP_OUTERRING || type == SHPP_FIRSTRING;
        Ring ring(outer ? Ring::POUT : Ring::PIN);
        ring.begin_poly();
        for (int j = st_id; j <= end_id; j++)
        {
            //Points may duplicate....
            double & x = shape->padfX[j];
            double & y = shape->padfY[j];
            ring.add_vertex(x, y);

            if (x < m_bbox[0]) m_bbox[0] = x;
            if (x > m_bbox[1]) m_bbox[1] = x;
            if (y < m_bbox[2]) m_bbox[2] = y;
            if (y > m_bbox[3]) m_bbox[3] = y;
        }
        ring.end_poly();
        add_polygon(ring, shape_index);
        started = true;
    }
}

void ShapeReader::read_arc(SHPObject * shape, int shape_index)
{
    //read each part
    int part_index = 1;

    //
    GisPolyline arc;
    arc.begin_poly();

    for (int j = 0; j < shape->nVertices; j++)
    {
        string part_type;

        if (j == 0 && shape->nParts > 0) part_type = SHPPartTypeName(shape->panPartType[0]);

        //this defines a new loop...
        if (part_index < shape->nParts && shape->panPartStart[part_index] == j)
        {
            part_type = SHPPartTypeName(shape->panPartType[part_index]);
            part_index++;
            arc.end_poly();
            arc.m_id = shape_index;
            arc.compute_geodesic_dist();
            m_arc_list.push_back(arc);

            //start a new polyline
            arc = GisPolyline();
            arc.begin_poly();
        }

        //Points may duplicate....
        double & x = shape->padfX[j];
        double & y = shape->padfY[j];
        arc.add_vertex(x, y);

        if (x < m_bbox[0]) m_bbox[0] = x;
        if (x > m_bbox[1]) m_bbox[1] = x;
        if (y < m_bbox[2]) m_bbox[2] = y;
        if (y > m_bbox[3]) m_bbox[3] = y;
    }

    arc.end_poly();
    arc.m_id = shape_index;
    arc.compute_geodesic_dist();
    m_arc_list.push_back(arc);
}

void ShapeReader::read_pts(SHPObject * shape, int shape_index)
{
    m_pt_list.reserve(shape->nVertices);
    for (int j = 0; j < shape->nVertices; j++)
    {
        //Points may duplicate....
        double & x = shape->padfX[j];
        double & y = shape->padfY[j];
        if (x < m_bbox[0]) m_bbox[0] = x;
        if (x > m_bbox[1]) m_bbox[1] = x;
        if (y < m_bbox[2]) m_bbox[2] = y;
        if (y > m_bbox[3]) m_bbox[3] = y;

        GisPoint2d pt(x, y);
        pt.m_id = shape_index;
        m_pt_list.push_back(pt);
    }
}

void ShapeReader::add_polygon(Ring & ring, int shape_index)
{
    if (ring.get_type() == Ring::POUT)
    {
        //create a new polygon
        GisPolygon polygon;
        polygon.m_id = shape_index;
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
