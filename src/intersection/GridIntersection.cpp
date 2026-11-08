/*! \file GridIntersection.cpp

\brief a class for computing the intersection between the grid and a geometric object

\author <a href="http://masc.cs.gmu.edu">MASC group</a>, George Mason University
\author <a href="http://profile.usgs.gov/langevin/">Christian Langevin</a>, USGS
\bug    No known bugs.
*/

#include "GridIntersection.h"
#include "shapefile/shpReader.h"
#include <algorithm>
#include <functional>
#include <vector>
#include <string>
#include "intersecReadWrite.h"
#include "RingClipper.h"
#include "PolygonChecker.h"
#include "GridgenError.h"
#include "Warnings.h"
using namespace std;

namespace cusg
{

//GridIntersection

void GridIntersection::build(GridIntersectionBlock & data)
{
    //check the refinement rules
    map<string, ShapeReader> shp_readers;

    //cast
    QuadTree3D * qtgrid = dynamic_cast<QuadTree3D *>(data.grid);
    ModflowGrid * mfgrid = dynamic_cast<ModflowGrid *>(data.grid);
    ModflowGrid2D * mf2dqtgrid = dynamic_cast<ModflowGrid2D *>(data.grid);

    //as least one of the grids needs to be valid
    assert(qtgrid || qtgrid || mf2dqtgrid);

    //rotate
    ModflowGrid * rmfgrid = mfgrid;
    if (rmfgrid == nullptr)
    {
        if (qtgrid != nullptr) rmfgrid = qtgrid->get_modflow_grid();
    }
    if (rmfgrid == nullptr)
        throw GridgenError("the grid of " + data.label() +
                           " is neither a MODFLOW grid nor a quadtree grid");
    double cx, cy, angle;
    rmfgrid->get_rotation(cx, cy, angle);

    //read from shapefile
    //ShapeReader shp_reader;
    string block = "GRID_INTERSECTION " + data.name;
    m_shp_reader.destroy();
    if (!m_shp_reader.read(data.shapefile))
        throw GridgenError("cannot open shapefile " + data.shapefile + " (" + block +
                           "); its .shp or .dbf file is missing or cannot be read");
    int nfeatures = 0, noutside = 0;

    //only support
    if (data.feature_type == "line")
    {
        vector<GisPolyline> & arcs = m_shp_reader.get_arc_list();
        //subdivide the qtree using the features
        nfeatures = arcs.size();
        for (int j = 0; j < nfeatures; j++)
        {
            int previous_seg_inter_num = m_seg_interaction_info.size();
            //intersect returns true even for a segment outside the grid, so count the cells found
            size_t ncells = m_intersection.size();

            int seg_size = arcs[j].get_size() - 1; //get_size returns the number of vertices
            for (int is = 0; is < seg_size; is++)
            {
                GisLineSeg2d line = arcs[j].get_line_segment(is);
                //rotate
                if (angle != 0) line.rotate(cx, cy, -angle);

                if (qtgrid != nullptr)
                    intersect_segment(qtgrid, line);
                else if (mfgrid != nullptr)
                    intersect_segment(mfgrid, line);
                else if (mf2dqtgrid != nullptr)
                    intersect_segment(mf2dqtgrid, line);
            } //end for i
            if (m_intersection.size() == ncells) noutside++;

            //merge the intersection segments that are inside the same box
            int after_seg_inter_num = m_seg_interaction_info.size();
            int k = previous_seg_inter_num;
            while (k < after_seg_inter_num)
            {
                PolylineCellIntersectionInfo polyline_info;
                polyline_info.box = m_seg_interaction_info[k].box;
                polyline_info.m_curve_id = m_seg_interaction_info[k].m_curve_id;
                polyline_info.m_start_geo_dist = m_seg_interaction_info[k].m_start_geo_dist;
                polyline_info.m_end_geo_dist = m_seg_interaction_info[k].m_end_geo_dist;
                polyline_info.m_len = m_seg_interaction_info[k].m_len;
                polyline_info.m_polyline = new Polyline();
                polyline_info.m_polyline->begin_poly();

                Point2d & start_point = m_seg_interaction_info[k].m_start_point;
                Point2d & end_point = m_seg_interaction_info[k].m_end_point;
                if (angle != 0)
                {
                    start_point.rotate(cx, cy, angle);
                    end_point.rotate(cx, cy, angle); //rotate back
                }
                polyline_info.m_polyline->add_vertex(start_point[0], start_point[1]);
                polyline_info.m_polyline->add_vertex(end_point[0], end_point[1]);

                while (k < (after_seg_inter_num - 1) &&
                       m_seg_interaction_info[k].m_curve_id ==
                           m_seg_interaction_info[k + 1].m_curve_id &&
                       m_seg_interaction_info[k].box == m_seg_interaction_info[k + 1].box &&
                       m_seg_interaction_info[k].m_end_point.almost_equ(
                           m_seg_interaction_info[k + 1].m_start_point))
                {
                    end_point = m_seg_interaction_info[k + 1].m_end_point;
                    if (angle != 0)
                    {
                        end_point.rotate(cx, cy, angle); //rotate back
                    }
                    polyline_info.m_polyline->add_vertex(end_point[0], end_point[1]);
                    polyline_info.m_end_geo_dist = m_seg_interaction_info[k + 1].m_end_geo_dist;
                    polyline_info.m_len += m_seg_interaction_info[k + 1].m_len;
                    k++;
                }

                polyline_info.m_polyline->end_poly();
                m_polyline_intersection_info.push_back(polyline_info);

                k++;
            }

        } //end for j
    } //end if
    else if (data.feature_type == "polygon")
    {
        int previous_size = m_polygon_intersection_info.size();
        vector<GisPolygon> & polygons = m_shp_reader.get_poly_list();
        //subdivide the qtree using the features
        nfeatures = polygons.size();
        for (int j = 0; j < nfeatures; j++)
        {

            if (j % 500 == 0) cout << "- Processing polygon: " << j << "    ......" << endl;

            //rotate
            if (angle != 0) polygons[j].rotate(cx, cy, -angle);

            vector<GisPolygon *> parts;
            check_hole(polygons[j], parts);

            //loop over all the polygons
            bool intersects = false;
            for (auto & part : parts)
            {
                GisPolygon polygon_part = *part;

                RingClipper * ring_clip = convert_polygon_to_clip(*part);

                if (qtgrid != nullptr)
                    intersects |= intersect(qtgrid, polygon_part, ring_clip); //polygons[j]);
                else if (mfgrid != nullptr)
                    intersects |= intersect(mfgrid, polygon_part, ring_clip); //polygons[j]);
                else if (mf2dqtgrid != nullptr)
                    intersects |= intersect(mf2dqtgrid, polygon_part, ring_clip); //polygons[j]);

                ring_clip->destroy();
                delete ring_clip;
            }
            //the rings are those of polygons[j], which the reader frees
            for (auto polygon_part : parts) delete polygon_part;
            if (!intersects) noutside++;
        }
        int after_size = m_polygon_intersection_info.size();
        if (angle != 0)
        {
            for (int k = previous_size; k < after_size; k++)
            {
                m_polygon_intersection_info[k].m_inter_polygon->rotate(cx, cy, angle);
            }
        }
    }
    else if (data.feature_type == "point")
    {
        int previous_size = m_pt_interaction_info.size();
        vector<GisPoint2d> & pts = m_shp_reader.get_point_list();
        //subdivide the qtree using the features
        nfeatures = pts.size();
        for (int j = 0; j < nfeatures; j++)
        {

            //rotate
            if (angle != 0) pts[j].rotate(cx, cy, -angle);

            bool intersects = false;
            if (qtgrid != nullptr)
                intersects = intersect(qtgrid, pts[j]);
            else if (mfgrid != nullptr)
                intersects = intersect(mfgrid, pts[j]);
            else if (mf2dqtgrid != nullptr)
                intersects = intersect(mf2dqtgrid, pts[j]);
            if (!intersects) noutside++;
        }
        //rotate back
        int after_size = this->m_pt_interaction_info.size();
        if (angle != 0)
        {
            for (int k = previous_size; k < after_size; k++)
            {
                m_pt_interaction_info[k].point.rotate(cx, cy, angle);
            }
        }

    } //end if

    if (nfeatures == 0)
        warn("Shapefile " + data.shapefile + " (" + block + ") has no " + data.feature_type +
             " features; it has no intersections with the grid.");
    warn_features_outside(noutside, nfeatures, data.shapefile, block, m_layer,
                          "no intersections are written for them");
}

//shape file
bool GridIntersection::save_shp(string & shp_name, vector<string> & attributes)
{
    if (shp_name.size() != 0)
    {

        if (m_pt_interaction_info.size() != 0)
        {
            write_pnt_intersect(this, shp_name, attributes, WRITE_SHP);
        }
        else if (m_seg_interaction_info.size() != 0)
        {
            write_line_intersect(this, shp_name, attributes, WRITE_SHP);
        }
        else if (m_polygon_intersection_info.size() != 0)
        {
            write_poly_intersect(this, shp_name, attributes, WRITE_SHP);
        }
        else
        {
            warn("No feature intersects the grid, so " + shp_name + " is not written.");
            return true;
        }
        cout << "- Save shape data to file: " << shp_name << endl;
    }

    return true;
}

//save txt file
bool GridIntersection::save(string & filename, vector<string> & attributes)
{
    if (filename.size() != 0)
    {

        if (m_pt_interaction_info.size() != 0)
        {
            write_pnt_intersect(this, filename, attributes, WRITE_TXT);
        }
        else if (m_seg_interaction_info.size() != 0)
        {
            write_line_intersect(this, filename, attributes, WRITE_TXT);
        }
        else if (m_polygon_intersection_info.size() != 0)
        {
            write_poly_intersect(this, filename, attributes, WRITE_TXT);
        }
        else
        {
            warn("No feature intersects the grid, so " + filename + " is not written.");
            return true;
        }
        cout << "- Save txt data to file: " << filename << endl;
    }
    else
    {
        throw GridgenError("no OUTPUT_FILE is given for the grid intersection");
    }

    return true;
}

//
//
// ModflowGrid2D
//
//

bool GridIntersection::intersect(ModflowGrid2D * grid, const Point2d & pt)
{
    Index2d id2d = grid->get_index(pt);

    if (id2d[0] == INT_MAX || id2d[1] == INT_MAX || id2d[0] == INT_MIN || id2d[1] == INT_MIN)
        return false;

    //get the box and put it in m_intersection
    int id = grid->get_box_index(id2d[0], id2d[1]);
    Box * box = grid->node_group[id];
    m_intersection.push_back(box);

    ////shpPntId.insert(make_pair(box, m_cur_id));

    return true;
}

bool GridIntersection::clip_line(ModflowGrid2D * grid, const LineSeg2d & line,
                                 LineSeg2d & clipped_line)
{
    double cx = grid->xoffset + grid->width / 2;
    double cy = grid->yoffset + grid->height / 2;

    Box box(cx, cy, grid->width, grid->height);

    bool found = intersect(&box, line);
    if (found == false) return false; //there is no intersection...

    found = get_intersection(&box, line, clipped_line);
    if (found == false) return false; //there is no intersection...

    if (box.in(clipped_line.pt[0][0], clipped_line.pt[0][1]) == false)
    {
        //move clipped_line.pt[0] toward clipped_line.pt[1] until the point is inside the box
        Point2d mid((clipped_line.pt[0][0] + clipped_line.pt[1][0]) / 2,
                    (clipped_line.pt[0][1] + clipped_line.pt[1][1]) / 2);
        bracket_to_box_boundary(box, clipped_line.pt[0], mid);
        clipped_line.pt[0] = mid;
    }

    if (box.in(clipped_line.pt[1][0], clipped_line.pt[1][1]) == false)
    {
        //move clipped_line.pt[1] toward clipped_line.pt[0]
        Point2d mid((clipped_line.pt[0][0] + clipped_line.pt[1][0]) / 2,
                    (clipped_line.pt[0][1] + clipped_line.pt[1][1]) / 2);
        bracket_to_box_boundary(box, clipped_line.pt[1], mid);
        clipped_line.pt[1] = mid;
    }

    return true;
}

bool GridIntersection::clip_line(QuadTree3D * grid, const LineSeg2d & line,
                                 LineSeg2d & clipped_line)
{
    return clip_line((ModflowGrid2D *)grid->get_modflow_grid(), line, clipped_line);
}

Point2d GridIntersection::grid_center(ModflowGrid2D * grid)
{
    return Point2d(grid->xoffset + grid->width / 2, grid->yoffset + grid->height / 2);
}

Point2d GridIntersection::grid_center(QuadTree3D * grid)
{
    return grid_center((ModflowGrid2D *)grid->get_modflow_grid());
}

void GridIntersection::bracket_to_box_boundary(Box & box, Point2d & out, Point2d & in) const
{
    //make sure that in point is in and out point is out of box
    assert(box.in(in[0], in[1]));
    assert(box.in(out[0], out[1]) == false);
    //
    while (true)
    {
        Point2d mid((in[0] + out[0]) / 2, (in[1] + out[1]) / 2);
        if (equal(mid.get(), in.get())) break;
        if (equal(mid.get(), out.get())) break;

        if (box.in(mid[0], mid[1]))
        {
            in = mid;
        }
        else
            out = mid;
    } //end while
}

//
//
// ModflowGrid
//
//
bool GridIntersection::intersect(ModflowGrid * grid, const Point2d & pt) //, bool noBound)
{
    Index2d id2d = grid->get_index(pt);
    if (id2d[0] == INT_MAX || id2d[1] == INT_MAX || id2d[0] == INT_MIN || id2d[1] == INT_MIN)
        return false;

    //get the box and put it in m_intersection
    int id = grid->get_box_index(this->m_layer, id2d[0], id2d[1]);

    Box * box = grid->node_group[id];
    m_intersection.push_back(box);

    return true;
}

//
//
//QuadTree3d
//
//

bool GridIntersection::intersect(QuadTree3D * grid, const Point2d & pt) //, bool noBoud)
{
    ModflowGrid * mf = grid->get_modflow_grid();

    Index2d id2d = mf->get_index(pt);

    if (id2d[0] == INT_MAX || id2d[1] == INT_MAX || id2d[0] == INT_MIN || id2d[1] == INT_MIN)
        return false;

    //get the box and put it in m_intersection
    int id = mf->get_box_index(m_layer, id2d[0], id2d[1]);
    Box * box = grid->node_group[id];

    box = box->find(pt[0], pt[1]); //find a leave that contains pt

    //make sure box is valid
    assert(box);

    //remember
    m_intersection.push_back(box);

    ////shpPntId.insert(make_pair(box, m_cur_id));//add by guilin

    //remember the boxes that are exactly on the boundary of the polygon
    if (box->on_bound(pt[0], pt[1])) this->m_on_bound.push_back(box);

    return true;
}

//check intersection between a box and a line
//TODO: this can be improved
bool GridIntersection::intersect(Box * b, const LineSeg2d & line)
{
    //check if the end points are in
    if (b->in(line.pt[0][0], line.pt[0][1])) return true;
    if (b->in(line.pt[1][0], line.pt[1][1])) return true;

    //create each edge of b
    const Point2d ul(b->x - b->dx / 2, b->y + b->dy / 2); //upper left corner
    const Point2d ur(b->x + b->dx / 2, b->y + b->dy / 2); //upper right corner
    const Point2d lr(b->x + b->dx / 2, b->y - b->dy / 2); //lower right corner
    const Point2d ll(b->x - b->dx / 2, b->y - b->dy / 2); //lower left corner

    const double * start = line.pt[0].get();
    const double * end = line.pt[1].get();

    //now, check intersections between the line segment and the edges
    if (seg_seg_int(ul.get(), ur.get(), start, end)) return true;
    if (seg_seg_int(ur.get(), lr.get(), start, end)) return true;
    if (seg_seg_int(lr.get(), ll.get(), start, end)) return true;
    if (seg_seg_int(ll.get(), ul.get(), start, end)) return true;

    //JML: I am not sure if the code below is really necessary.
    if (area_sign(ul.get(), ur.get(), start) == 0)
        if (between(ul.get(), ur.get(), start)) return true;
    if (area_sign(ur.get(), lr.get(), start) == 0)
        if (between(ur.get(), lr.get(), start)) return true;
    if (area_sign(lr.get(), ll.get(), start) == 0)
        if (between(lr.get(), ll.get(), start)) return true;
    if (area_sign(ll.get(), ul.get(), start) == 0)
        if (between(ll.get(), ul.get(), start)) return true;

    if (area_sign(ul.get(), ur.get(), end) == 0)
        if (between(ul.get(), ur.get(), end)) return true;
    if (area_sign(ur.get(), lr.get(), end) == 0)
        if (between(ur.get(), lr.get(), end)) return true;
    if (area_sign(lr.get(), ll.get(), end) == 0)
        if (between(lr.get(), ll.get(), end)) return true;
    if (area_sign(ll.get(), ul.get(), end) == 0)
        if (between(ll.get(), ul.get(), end)) return true;

    //otherwise
    return false;
}

void GridIntersection::build_intersection_info(Box * box, const GisPolygon & polygon,
                                               vector<PolygonCellIntersectionInfo> & info_result,
                                               RingClipper * ring_clip)
{
    vector<GisPolygon *> result_polygons;
    intersect_polygon(box, polygon, result_polygons, ring_clip);
    for (auto polygon_part : result_polygons)
    {
        PolygonCellIntersectionInfo info;
        info.m_inter_polygon = polygon_part;
        info.box = box;
        info.m_poly_id = polygon.m_id;
        info.m_area = info.m_inter_polygon->get_area();

        info_result.push_back(info);
    }
}

SegmentCellIntersectionInfo GridIntersection::build_intersection_info(Box * box,
                                                                      const GisLineSeg2d & line)
{
    SegmentCellIntersectionInfo info;
    info.box = box;
    info.m_seg_id = line.m_seg_id;
    info.m_curve_id = line.m_curve_id;

    LineSeg2d clipped_segment;
    bool found = get_intersection(box, (LineSeg2d)line, clipped_segment);

    if (found == false)
    {
        warn("The part of a line segment inside a cell it crosses cannot be found, so it is "
             "written with zero length.");
        return info;
    }

    info.m_len = (clipped_segment.pt[0] - clipped_segment.pt[1]).norm();
    info.m_start_geo_dist = (clipped_segment.pt[0] - line.pt[0]).norm() + line.m_start_geo_dist;
    info.m_end_geo_dist = (clipped_segment.pt[1] - line.pt[0]).norm() + line.m_start_geo_dist;

    if (info.m_start_geo_dist > info.m_end_geo_dist)
    {
        swap(info.m_end_geo_dist, info.m_start_geo_dist);
        info.m_start_point = clipped_segment.pt[1];
        info.m_end_point = clipped_segment.pt[0];
    }
    else
    {
        info.m_start_point = clipped_segment.pt[0];
        info.m_end_point = clipped_segment.pt[1];
    }

    return info;
}

PointCellIntersectionInfo GridIntersection::build_intersection_info(Box * box,
                                                                    const GisPoint2d & pt)
{
    PointCellIntersectionInfo info;
    info.box = box;
    info.m_point_id = pt.m_id;
    info.point = pt;
    return info;
}

bool GridIntersection::get_intersection(Box * b, const LineSeg2d & line,
                                        LineSeg2d & clipped_segment)
{

    //check if the end points are in
    bool end_in_box[2] = {false};
    if (b->in(line.pt[0][0], line.pt[0][1])) end_in_box[0] = true;
    if (b->in(line.pt[1][0], line.pt[1][1])) end_in_box[1] = true;

    if (end_in_box[0] && end_in_box[1])
    {
        clipped_segment = line;
        return true;
    }

    //create each edge of b
    Point2d ul(b->x - b->dx / 2, b->y + b->dy / 2); //upper left corner
    Point2d ur(b->x + b->dx / 2, b->y + b->dy / 2); //upper right corner
    Point2d lr(b->x + b->dx / 2, b->y - b->dy / 2); //lower right corner
    Point2d ll(b->x - b->dx / 2, b->y - b->dy / 2); //lower left corner

    //now, check intersections between the line segment and the edges
    char edge_code[4];
    Point2d p[4];

    edge_code[0] = seg_seg_int(ul.get(), ur.get(), line.pt[0].get(), line.pt[1].get(), p[0].get());
    edge_code[1] = seg_seg_int(ur.get(), lr.get(), line.pt[0].get(), line.pt[1].get(), p[1].get());
    edge_code[2] = seg_seg_int(lr.get(), ll.get(), line.pt[0].get(), line.pt[1].get(), p[2].get());
    edge_code[3] = seg_seg_int(ll.get(), ul.get(), line.pt[0].get(), line.pt[1].get(), p[3].get());

    //no intersection....
    if ((edge_code[0] || edge_code[1] || edge_code[2] || edge_code[3]) == false) return false;

    if (end_in_box[0] || end_in_box[1])
    {
        float max_len = 0;
        for (int i = 0; i < 4; i++)
        {
            if (edge_code[i] == '0') continue;
            if (edge_code[i] == 'e')
            {
                //compute the union between the edge of box and the line
                Point2d seg[2];
                const double * u[2] = {nullptr, nullptr};

                switch (i)
                {
                    case 0:
                        seg[0] = ul;
                        seg[1] = ur;
                        break;
                    case 1:
                        seg[0] = ur;
                        seg[1] = lr;
                        break;
                    case 2:
                        seg[0] = lr;
                        seg[1] = ll;
                        break;
                    case 3:
                        seg[0] = ll;
                        seg[1] = ul;
                        break;
                }

                bool overlaps = segment_union(seg[0].get(), seg[1].get(), line.pt[0].get(),
                                              line.pt[1].get(), u);
                if (overlaps)
                {
                    clipped_segment.pt[0].set(u[0][0], u[0][1]);
                    clipped_segment.pt[1].set(u[1][0], u[1][1]);
                }

                return overlaps;
            }
            //found the intersection...
            const Point2d & pt1 = p[i];
            const Point2d & pt2 = (end_in_box[0]) ? line.pt[0] : line.pt[1];

            double len = (pt1 - pt2).normsqr();
            if (len > max_len)
            {
                clipped_segment.pt[0] = pt1;
                clipped_segment.pt[1] = pt2;
                max_len = len;
            }
        }
        //done

        return true;
    }

    //count
    float max_len = 0;
    for (int i = 0; i < 4; i++)
    {
        if (edge_code[i] == '0') continue;

        if (edge_code[i] == 'e')
            continue; //this can't happen, it should be handled in the previous loop

        for (int j = i; j < 4; j++)
        {
            if (edge_code[j] == '0') continue;
            if (edge_code[j] == 'e') continue;

            //found the intersection...
            double len = (p[i] - p[j]).normsqr();
            if (len > max_len)
            {
                clipped_segment.pt[0] = p[i];
                clipped_segment.pt[1] = p[j];
                max_len = len;
            }
        }
    }

    return true;
}

double GridIntersection::intersection_length(Box * b, const LineSeg2d & line)
{
    //check if the end points are in
    LineSeg2d clipped_segment;
    bool found = get_intersection(b, line, clipped_segment);
    if (found) return (clipped_segment.pt[0] - clipped_segment.pt[1]).norm();
    return 0; //no intersection...
}

//
// check if a point is inside a ring
// ray shooting is slow. Can be improved.
//
bool GridIntersection::intersect(const Point2d & pt, const Ring & ring)
{

    int max_trial = 200;

    Vector2d offset(ring.get_head()->get_pos().get());
    float ray_length = ring.get_radius() * 4;
    Point2d near = pt - offset;
    Point2d far = near;

    while (max_trial > 0)
    {
        //get a random point
        Vector2d vec(drand48(), drand48());
        if (drand48() < 0.5) vec[0] = -vec[0];
        if (drand48() < 0.5) vec[1] = -vec[1];
        far = pt - offset + vec.normalize() * ray_length;

        //
        LinkedVertex * ptr = ring.get_head();
        bool degenerate = false;
        int int_count = 0;

        do
        {
            LinkedVertex * next = ptr->get_next();
            double x[2];
            Point2d edge_start = ptr->get_pos() - offset;
            Point2d edge_end = next->get_pos() - offset;

            char code = seg_seg_int(near.get(), far.get(), edge_start.get(), edge_end.get(), x);

            if (code == '1')
                int_count++;
            else if (code != '0')
            {
                degenerate = true;
                break;
            }

            ptr = next;
        } while (ptr != ring.get_head());

        if (degenerate)
        {
            max_trial--;
            continue;
        }
        return (int_count % 2) == 1;
    }

    return false;
}

//
// check if a point is inside a polygon
// ray shooting is slow. Can be improved.
//
bool GridIntersection::intersect(const Point2d & pt, const Polygon & pgon)
{

    int max_trial = 200;

    Vector2d offset(pgon.front().get_head()->get_pos().get());
    const double * pgon_box = pgon.get_bounding_box();
    float ray_length = max(pgon_box[1] - pgon_box[0], pgon_box[3] - pgon_box[2]) * 4;

    Point2d near = pt - offset;
    Point2d far = near;

    while (max_trial > 0)
    {
        //get a random point
        Vector2d vec(drand48(), drand48());
        if (drand48() < 0.5) vec[0] = -vec[0];
        if (drand48() < 0.5) vec[1] = -vec[1];
        far = pt - offset + vec.normalize() * ray_length;

        //count the number of intersection
        bool degenerate = false; //flag raised if degenerate intersection is found.
        int int_count = 0;       //# of intersections

        //loop through each ring
        for (const auto & ring : pgon)
        {
            LinkedVertex * ptr = ring.get_head();

            do
            {
                LinkedVertex * next = ptr->get_next();
                double x[2];
                Point2d edge_start = ptr->get_pos() - offset;
                Point2d edge_end = next->get_pos() - offset;

                char code = seg_seg_int(near.get(), far.get(), edge_start.get(), edge_end.get(), x);

                if (code == '1')
                    int_count++;
                else if (code != '0')
                {
                    degenerate = true;
                    break;
                }

                ptr = next;
            } while (ptr != ring.get_head());

            if (degenerate) break; //something is wrong
        }

        if (degenerate)
        {
            max_trial--;
            continue;
        }

        return (int_count % 2) == 1;
    }

    return false;
}

//void GridIntersection::intersect_polygon(Box* b, const GisPolygon& polygon, vector<GisPolygon*>& result_polygons)
void GridIntersection::intersect_polygon(Box * b, const GisPolygon & polygon,
                                         vector<GisPolygon *> & result_polygons,
                                         RingClipper * initial_clip)
{

    double xmin = min(b->x - b->dx / 2, b->x + b->dx / 2);
    double xmax = max(b->x - b->dx / 2, b->x + b->dx / 2);
    double ymin = min(b->y - b->dy / 2, b->y + b->dy / 2);
    double ymax = max(b->y - b->dy / 2, b->y + b->dy / 2);

    /*****************************************************************/
    RingClipper * ring_clip = initial_clip;
    if (initial_clip == nullptr)
    {
        ring_clip = convert_polygon_to_clip(polygon);
    }

    list<RingClipper *> clips;
    ring_clip->clip(xmin, xmax, ymin, ymax, clips);

    for (auto & clip_part : clips)
    {
        GisPolygon * npolygons = new GisPolygon();
        convert_clip_to_polygon(clip_part, npolygons);
        result_polygons.push_back(npolygons);

        clip_part->destroy();
        delete clip_part;
    }

    ring_clip->reset_status();

    if (initial_clip == nullptr)
    {
        ring_clip->destroy();
        delete ring_clip;
    }
    else
    {
        ring_clip->reset_status();
    }
}

//the side of box that other shares an edge with, or -1 if they are not neighbors
int GridIntersection::shared_edge_side(Box * box, Box * other)
{
    for (int i = 0; i < 4; i++)
    {
        list<Box *> neighbors;
        box->get_neighbors(neighbors, i);
        if (std::find(neighbors.begin(), neighbors.end(), other) != neighbors.end()) return i;
    }

    return -1;
}

//check box edge collinearly intersect with line
bool GridIntersection::line_runs_along_edge(Box * box, int dir, const LineSeg2d & line)
{
    double edge_start[2], edge_end[2];
    double line_start[2] = {line.pt[0][0], line.pt[0][1]};
    double line_end[2] = {line.pt[1][0], line.pt[1][1]};
    double p[2];
    if (dir == 0)
    {
        edge_start[0] = box->x - box->dx / 2;
        edge_end[0] = box->x + box->dx / 2;
        edge_start[1] = box->y + box->dy / 2;
        edge_end[1] = box->y + box->dy / 2;
    }
    else if (dir == 1)
    {
        edge_start[0] = box->x + box->dx / 2;
        edge_end[0] = box->x + box->dx / 2;
        edge_start[1] = box->y + box->dy / 2;
        edge_end[1] = box->y - box->dy / 2;
    }
    else if (dir == 2)
    {
        edge_start[0] = box->x + box->dx / 2;
        edge_end[0] = box->x - box->dx / 2;
        edge_start[1] = box->y - box->dy / 2;
        edge_end[1] = box->y - box->dy / 2;
    }
    else
    {
        assert(dir == 3);
        edge_start[0] = box->x - box->dx / 2;
        edge_end[0] = box->x - box->dx / 2;
        edge_start[1] = box->y + box->dy / 2;
        edge_end[1] = box->y - box->dy / 2;
    }
    //check collinear
    if (seg_seg_int(edge_start, edge_end, line_start, line_end, p) == 'e') return true;

    return false;
}

///collect all cells in the grid from b
void GridIntersection::flood_grid(Box * start_box, set<Box *> & visited, list<Box *> & flooded)
{
    list<Box *> open;
    open.push_back(start_box);

    while (open.empty() == false) //not empty
    {
        Box * b = open.front();
        open.pop_front();
        flooded.push_back(b);

        //collect neighboring cells
        list<Box *> neighbors[4];

        for (int i = 0; i < 4; i++)
        {
            b->get_neighbors(neighbors[i], i);
        }

        //get all unvisited cells
        for (auto & side : neighbors)
        {
            for (list<Box *>::iterator n = side.begin(); n != side.end(); n++)
            {
                Box * neighbor = *n;
                if (visited.find(neighbor) != visited.end()) continue; //visited
                visited.insert(neighbor);
                open.push_back(neighbor);
            } //end for n
        } //end for i
    } //end while

} //end GridIntersection::flood_grid
}
