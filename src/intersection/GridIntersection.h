/*! \file GridIntersection.h

\brief a class for computing the intersection between the grid and a geometric object

\author <a href="http://masc.cs.gmu.edu">MASC group</a>, George Mason University
\author <a href="http://profile.usgs.gov/langevin/">Christian Langevin</a>, USGS
\bug    No known bugs.
*/

#pragma once

#include "ModflowGrid.h"
#include "QuadTree3D.h"
#include "intersection.h"
#include "lineseg.h"
#include "Point.h"
#include "polygon.h"
#include "polyline.h"
#include "def_struct.h"
#include "shpReader.h"
#include "RingClipper.h"
#include "PolygonChecker.h"
#include "Warnings.h"
#include <map>
using namespace std;

namespace cusg
{

inline double get_sq_length(const LineSeg2d & line)
{
    double dx2 = (line.pt[0][0] - line.pt[1][0]) * (line.pt[0][0] - line.pt[1][0]);
    double dy2 = (line.pt[0][1] - line.pt[1][1]) * (line.pt[0][1] - line.pt[1][1]);
    return dx2 + dy2;
}

struct PointCellIntersectionInfo
{
    int m_point_id;
    Box * box;
    GisPoint2d point;
};

struct PolylineCellIntersectionInfo
{
    int m_curve_id;
    double m_len;
    double m_start_geo_dist;
    double m_end_geo_dist;
    Box * box;
    Polyline * m_polyline; //the cut polyline
};

struct SegmentCellIntersectionInfo
{
    SegmentCellIntersectionInfo()
    {
        m_curve_id = m_seg_id = -1;
        m_len = m_start_geo_dist = m_end_geo_dist = 0;
        box = nullptr;
    }

    int m_curve_id; //id of the curve that contains this segment
    int m_seg_id;   //id of the segment
    double m_len;   //length of the segment inside the box
    double m_start_geo_dist;
    double m_end_geo_dist;
    Box * box;
    Point2d m_start_point; // the starting point
    Point2d m_end_point;   // the ending point
};

struct PolygonCellIntersectionInfo
{
    PolygonCellIntersectionInfo()
    {
        m_poly_id = -1;
        m_area = 0;
        box = nullptr;
        m_inter_polygon = nullptr;
    }

    int m_poly_id;
    double m_area; //intersection area
    Box * box;
    GisPolygon * m_inter_polygon;
};

class GridIntersection
{
public:
    GridIntersection(int max_level, int at_layer)
    {
        m_max_level = max_level;
        m_layer = at_layer;
        m_b_include_poly_boundary = true;
        m_b_save_intersection_info = false;
    }
    ~GridIntersection()
    {
        for (auto & pit : m_polygon_intersection_info)
        {
            pit.m_inter_polygon->destroy();

            delete pit.m_inter_polygon;
            pit.m_inter_polygon = nullptr;
        }
        for (auto & lit : m_polyline_intersection_info)
        {
            lit.m_polyline->destroy();

            delete lit.m_polyline;
            lit.m_polyline = nullptr;
        }
    }

    ///turn on/off m_b_include_poly_boundary to enable/disable storing
    ///the intersected boxes between the grid and the given polygon
    void set_include_poly_boundary_flag(bool flag) { m_b_include_poly_boundary = flag; }

    ///turn on/off m_b_save_intersection_info to enable/disable storing
    ///the detailed intersection information
    void set_save_intersection_info(bool flag) { m_b_save_intersection_info = flag; }

    //build intersection from grid_intersection raw data
    void build(GridIntersectionBlock & data);

    //save intersection information to file
    bool save(string & filename, vector<string> & attributes);
    void save(ostream & out, vector<string> & attributes);

    ///get intersected boxes
    vector<Box *> & get_intersections() { return m_intersection; }

    //save the shape file
    bool save_shp(string & shp_name, vector<string> & attributes);

    //
    bool intersect(QuadTree3D * grid, double x, double y) { return intersect(grid, Point2d(x, y)); }

    bool intersect_segment(QuadTree3D * grid, double x1, double y1, double x2, double y2)
    {
        Point2d p1(x1, y1);
        Point2d p2(x2, y2);
        LineSeg2d segment(p1, p2);

        return intersect_segment(grid, segment);
    }

    bool intersect(QuadTree3D * grid, vector<double> & coordinates)
    {
        //create ring
        Ring ring(Ring::POUT);
        int size = coordinates.size();
        ring.begin_poly();
        for (int i = 0; i < size; i++)
        {
            double x = coordinates[i];
            double y = coordinates[++i];
            ring.add_vertex(x, y);
        }
        ring.end_poly();

        return intersect(grid, ring);
    }

    Box * get_intersections(int i) { return m_intersection[i]; }

    // in all intersection functions, true is returned if
    // there are intersections; otherwise false is returned

    //
    // ModflowGrid2D
    //
    bool intersect(ModflowGrid2D * grid, const Point2d & pt);

    //
    // ModflowGrid, on_bound=true means ignore the boxes only intersecting at box boundary
    //
    bool intersect(ModflowGrid * grid, const Point2d & pt); //, bool noBound=false);

    //
    //QuadTree3d, on_bound=true means ignore the boxes only intersecting at box boundary
    //
    bool intersect(QuadTree3D * grid, const Point2d & pt); //, bool noBound=false);

    //add the cells that a line segment crosses to the intersection; false if an end
    //of the segment inside the grid is in no cell, and true otherwise, including
    //when the segment misses the grid
    template <typename T> bool intersect_segment(T * grid, const LineSeg2d & line);

    //intersect a generic Grid with a polygon..
    template <typename T> bool intersect(T * grid, const Ring & ring); //, bool noBound=false);
    template <typename T>
    bool intersect(T * grid, const Polygon & polygon); //, bool noBound=false);

    //intersect a generic Grid with a ring..
    template <typename T>
    bool intersect(T * grid, const Ring & ring, list<Box *> & bd_intersection,
                   list<Box *> & enclosed_intersection);

    //intersect a generic Grid with a polygon..
    template <typename T>
    bool intersect(T * grid, const Polygon & polygon, list<Box *> & bd_intersection,
                   list<Box *> & enclosed_intersection);

    bool intersect(Box * b, const LineSeg2d & line);
    double intersection_length(Box * b, const LineSeg2d & line);
    void intersect_polygon(Box * b, const GisPolygon & polygon,
                           vector<GisPolygon *> & result_polygons,
                           RingClipper * ring_clip = nullptr);

    //retrieve detailed information
    bool get_intersection(Box * b, const LineSeg2d & line, LineSeg2d & clipped_segment);

    //get the number of intersections...
    int get_intersection_size() { return m_intersection.size(); }

    //intersection for GIS based primitives
    template <typename T> bool intersect(T * grid, /*const*/ GisPoint2d & pt);
    template <typename T> bool intersect_segment(T * grid, /*const*/ GisLineSeg2d & line);
    template <typename T>
    bool intersect(T * grid, /*const*/ GisPolygon & polygon, RingClipper * ring_clip = nullptr);

    /// clip off a line segment if its end points are out side of the grid
    /// return false when the entire segment is outside of the grid
    bool clip_line(ModflowGrid2D * grid, const LineSeg2d & line, LineSeg2d & clipped_line);
    bool clip_line(QuadTree3D * grid, const LineSeg2d & line, LineSeg2d & clipped_line);

    //center of the grid, in the coordinates of clip_line
    Point2d grid_center(ModflowGrid2D * grid);
    Point2d grid_center(QuadTree3D * grid);

    //cells of the whole layer if the polygon encloses the grid, when no
    //polygon edge crosses it
    template <typename T, typename P>
    void find_cells_inside(T * grid, const P & polygon, set<Box *> & visited,
                           list<Box *> & enclosed_intersection);

    //utility functions
    int shared_edge_side(Box * box, Box * other);
    //check box edge collinearly intersect with line
    bool line_runs_along_edge(Box * box, int dir, const LineSeg2d & line);

protected:
    //check if the given pt is inside the ring
    bool intersect(const Point2d & pt, const Ring & ring);

    //check if the given pt is inside the polygon
    bool intersect(const Point2d & pt, const Polygon & polygon);

    //create information for intersection
    //the given box should be intersecting with the given geometry already
    void build_intersection_info(Box * box, const GisPolygon & polygon,
                                 vector<PolygonCellIntersectionInfo> & info_result,
                                 RingClipper * ring_clip = nullptr);

    SegmentCellIntersectionInfo build_intersection_info(Box * box, const GisLineSeg2d & line);
    PointCellIntersectionInfo build_intersection_info(Box * box, const GisPoint2d & pt);

    void flood_grid(Box * start_box, set<Box *> & visited, list<Box *> & flooded);

    //given a box and two points, one inside and one outside.
    //close the gap between in and out so that their distance is smaller than
    //a given threshold
    void bracket_to_box_boundary(Box & box, Point2d & out, Point2d & in) const;

    int m_max_level;
    int m_layer;
    vector<Box *> m_intersection; /// a list of intersected boxes

    bool m_b_include_poly_boundary;

    bool m_b_save_intersection_info;

    vector<PolygonCellIntersectionInfo> m_polygon_intersection_info;
    vector<SegmentCellIntersectionInfo> m_seg_interaction_info;
    vector<PointCellIntersectionInfo> m_pt_interaction_info;
    vector<PolylineCellIntersectionInfo> m_polyline_intersection_info;

    map<Box *, list<LinkedVertex *>>
        m_box_vertex_map; //mapping a box to a list of polygon vertices/edges inside the box

    //
    vector<Box *> m_on_bound;

public:
    ShapeReader m_shp_reader;

    vector<Box *> test_boxes;

    vector<PolygonCellIntersectionInfo> & get_polygon_intersection_info()
    {
        return m_polygon_intersection_info;
    }

    vector<SegmentCellIntersectionInfo> & get_segment_intersection_info()
    {
        return m_seg_interaction_info;
    }

    vector<PointCellIntersectionInfo> & get_point_intersection_info()
    {
        return m_pt_interaction_info;
    }

    vector<PolylineCellIntersectionInfo> & get_polyline_intersection_info()
    {
        return m_polyline_intersection_info;
    }
};

//
//
// implementation of templated functions
//
//

//intersection a generic Grid with a polygon..
template <typename T>
bool GridIntersection::intersect(T * grid, const Ring & ring) //, bool noBound)
{
    list<Box *> bd_intersection;
    list<Box *> enclosed_intersection;
    bool found = intersect(grid, ring, bd_intersection, enclosed_intersection);

    if (found)
    {
        //add the newly discovered intersections
        if (m_b_include_poly_boundary)
            m_intersection.insert(m_intersection.end(), bd_intersection.begin(),
                                  bd_intersection.end());

        m_intersection.insert(m_intersection.end(), enclosed_intersection.begin(),
                              enclosed_intersection.end());
    }

    return found;
}

//intersection a generic Grid with a polygon..
template <typename T>
bool GridIntersection::intersect(T * grid, const Polygon & polygon) //, bool noBound)
{
    list<Box *> bd_intersection;
    list<Box *> enclosed_intersection;
    bool found = intersect(grid, polygon, bd_intersection, enclosed_intersection);
    if (found)
    {
        //add the newly discovered intersections
        if (m_b_include_poly_boundary)
            m_intersection.insert(m_intersection.end(), bd_intersection.begin(),
                                  bd_intersection.end());
        m_intersection.insert(m_intersection.end(), enclosed_intersection.begin(),
                              enclosed_intersection.end());
    }
    return found;
}

//intersection a Grid with a ring..
template <typename T>
bool GridIntersection::intersect(T * grid, const Ring & ring, list<Box *> & bd_intersection,
                                 list<Box *> & enclosed_intersection)
{
    vector<Box *> backup;

    backup.swap(m_intersection);
    set<Box *> visited;

    //get all boundary cells
    LinkedVertex * ptr = ring.get_head();
    do
    {
        LinkedVertex * next = ptr->get_next();

        LineSeg2d line(ptr->get_pos(), next->get_pos());
        bool found = intersect_segment(grid, line);
        assert(found);

        //record the relationship between the boxes and this segment
        if (m_b_save_intersection_info)
        {
            for (auto & i : m_intersection) m_box_vertex_map[i].push_back(ptr);
        }

        //boundary cells in the order found, the set only skips duplicates
        for (auto & box : m_intersection)
        {
            if (visited.insert(box).second) bd_intersection.push_back(box);
        }
        m_intersection.clear();
        ptr = next;
    } while (ptr != ring.get_head());

    //no edge crosses the grid, so the polygon encloses all of the layer or none of it
    if (visited.empty())
    {
        find_cells_inside(grid, ring, visited, enclosed_intersection);
        backup.swap(m_intersection);
        return !enclosed_intersection.empty();
    }

    //get all enclosed cells
    list<Box *> open = bd_intersection;
    while (open.empty() == false) //not empty
    {
        Box * b = open.front();
        open.pop_front();

        //collect neighboring cells
        list<Box *> neighbors[4];

        for (int i = 0; i < 4; i++)
        {
            b->get_neighbors(neighbors[i], i);
        }

        //check if cells are inside the polygon
        for (auto & side : neighbors)
        {
            for (list<Box *>::iterator n = side.begin(); n != side.end(); n++)
            {
                Box * neighbor = *n;

                if (visited.find(neighbor) != visited.end()) continue; //visited
                Point2d box_center(neighbor->x, neighbor->y);          //center of the box
                if (intersect(box_center, ring)) //if the center of the box is inside...
                {
                    open.push_back(neighbor);
                    visited.insert(neighbor);
                    enclosed_intersection.push_back(neighbor);
                }
            } //end for n
        } //end for i

    } //end while

    //restore
    backup.swap(m_intersection);

    return true;
}

template <typename T, typename P>
void GridIntersection::find_cells_inside(T * grid, const P & polygon, set<Box *> & visited,
                                         list<Box *> & enclosed_intersection)
{
    Point2d center = grid_center(grid);
    if (!intersect(center, polygon)) return;

    //flood the layer from the cell at the center
    vector<Box *> center_cells;
    center_cells.swap(m_intersection);
    bool found = intersect(grid, center);
    center_cells.swap(m_intersection);
    if (!found || center_cells.empty()) return;
    visited.insert(center_cells.back());
    flood_grid(center_cells.back(), visited, enclosed_intersection);
}

//intersection a Grid with a polygon..
template <typename T>
bool GridIntersection::intersect(T * grid, const Polygon & polygon, list<Box *> & bd_intersection,
                                 list<Box *> & enclosed_intersection)
{
    vector<Box *> backup;

    backup.swap(m_intersection);
    set<Box *> visited;

    //get all boundary cells
    for (const auto & ring : polygon)
    {
        LinkedVertex * ptr = ring.get_head();

        do
        {
            LinkedVertex * next =
                ptr->get_next(); //ptr is the current vertex and next is the next vertex in this ring

            LineSeg2d line(ptr->get_pos(), next->get_pos());
            bool found = intersect_segment(grid, line);

            assert(found);

            //record the relationship between the boxes and this segment
            if (m_b_save_intersection_info)
            {
                for (auto & box : m_intersection) m_box_vertex_map[box].push_back(ptr);
            }

            //boundary cells in the order found, the set only skips duplicates
            for (auto & box : m_intersection)
            {
                if (visited.insert(box).second) bd_intersection.push_back(box);
            }
            m_intersection.clear();
            ptr = next;
        } while (ptr != ring.get_head());
    }

    //no edge crosses the grid, so the polygon encloses all of the layer or none of it
    if (visited.empty())
    {
        find_cells_inside(grid, polygon, visited, enclosed_intersection);
        backup.swap(m_intersection);
        return !enclosed_intersection.empty();
    }

    //get all enclosed cells
    list<Box *> open = bd_intersection;
    while (open.empty() == false) //not empty
    {
        Box * b = open.front();
        open.pop_front();

        //collect neighboring cells
        list<Box *> neighbors[4];

        for (int i = 0; i < 4; i++)
        {
            b->get_neighbors(neighbors[i], i);
        }

        //check if cells are inside the polygon
        for (auto & side : neighbors)
        {
            for (list<Box *>::iterator n = side.begin(); n != side.end(); n++)
            {
                Box * neighbor = *n;

                if (visited.find(neighbor) != visited.end()) continue; //visited
                Point2d box_center(neighbor->x, neighbor->y);          //center of the box
                if (intersect(box_center, polygon)) //if the center of the box is inside...
                {
                    open.push_back(neighbor);
                    visited.insert(neighbor);
                    enclosed_intersection.push_back(neighbor);
                }
            } //end for n
        } //end for i

    } //end while

    //restore
    backup.swap(m_intersection);

    return true;
}

//modify here to solve the problem of non-symmetrical refinement

template <typename T> bool GridIntersection::intersect_segment(T * grid, const LineSeg2d & line)
{
    LineSeg2d clipped_line;
    bool clipped = clip_line(grid, line, clipped_line);
    if (clipped == false)
        return true; //line segment is entirely outside, so there is no intersection
    if (get_sq_length(clipped_line) < 1e-10) return true;

    //get the first box that encloses one of the end pt
    bool found = intersect(grid, clipped_line.pt[0]);
    if (found == false)
    {
        warn("A line segment has an end inside the grid that is in no cell, so the segment is "
             "not intersected with the grid.");
        return false;
    }

    //first box
    Box * start_box = m_intersection.back();

    found = intersect(grid, clipped_line.pt[1]);
    if (found == false)
    {
        warn("A line segment has an end inside the grid that is in no cell, so the segment is "
             "not intersected with the grid.");
        return false;
    }

    //last box
    Box * end_box = m_intersection.back();

    //the line segment is short and contained in the same box
    if (start_box == end_box)
    {
        m_intersection.pop_back();
        return true; //done
    }

    //check if start_box and end_box share an edge that the line runs along
    int neighbor_id = shared_edge_side(start_box, end_box);
    if (neighbor_id != -1)
    {
        //check if the box edge collinearly intersect with the line
        if (line_runs_along_edge(start_box, neighbor_id, line))
        {
            m_intersection.pop_back();
            return true; //done
        }
    }

    //visit the neighbors the line crosses, starting from both end boxes
    list<Box *> open;
    open.push_back(start_box);
    open.push_back(end_box);

    set<Box *> visited;
    visited.insert(start_box);
    visited.insert(end_box);

    while (open.empty() == false) //not empty
    {
        Box * b = open.front();
        open.pop_front();

        //get the neighbors of the box b
        list<Box *> neighbors[4];

        for (int i = 0; i < 4; i++)
        {
            b->get_neighbors(neighbors[i], i);
        }

        //process neighbors
        for (int i = 0; i < 4; i++)
        {
            for (list<Box *>::iterator n = neighbors[i].begin(); n != neighbors[i].end(); n++)
            {
                Box * neighbor = *n;

                if (visited.find(neighbor) != visited.end()) continue; //visited

                if (intersect(neighbor, line))
                {
                    //////////////////////////////////////////////////////////////////////////
                    //add by Guilin
                    //if collinear intersection
                    //a line along the shared edge is assigned to the cell b
                    if (line_runs_along_edge(b, i, line)) continue;

                    double epsilon = min(neighbor->dx, neighbor->dy) / 1000;
                    if (intersection_length(neighbor, line) > epsilon) //long enough...
                    {

                        m_intersection.push_back(neighbor);
                    }

                    open.push_back(neighbor);
                    visited.insert(neighbor);
                }
            }
        }

    } //end while

    return true;
}

/**************************************************************************************************/
//intersection for GIS based primitives
template <typename T> bool GridIntersection::intersect(T * grid, /*const*/ GisPoint2d & pt)
{
    ////rotate
    ////rotate

    vector<Box *> previous;
    m_intersection.swap(previous);
    bool found = intersect(grid, (const Point2d &)pt);

    //process information
    if (this->m_b_save_intersection_info && found)
    {
        //store intersection info here...
        for (auto & i : m_intersection)
        {
            m_pt_interaction_info.push_back(build_intersection_info(i, pt));
        }
    }

    previous.insert(previous.end(), m_intersection.begin(), m_intersection.end());
    m_intersection.swap(previous);
    return found;
}

template <typename T>
bool GridIntersection::intersect_segment(T * grid, /*const*/ GisLineSeg2d & line)
{
    ////rotate
    ////rotate

    vector<Box *> previous;
    m_intersection.swap(previous);
    bool found = intersect_segment(grid, (const LineSeg2d &)line);

    //process information
    if (this->m_b_save_intersection_info && found)
    {
        //store intersection info here...
        for (auto & i : m_intersection)
        {
            m_seg_interaction_info.push_back(build_intersection_info(i, line));
        } //end for i

    } //end if

    previous.insert(previous.end(), m_intersection.begin(), m_intersection.end());
    m_intersection.swap(previous);
    return found;
}

template <typename T>
bool GridIntersection::intersect(T * grid, /*const*/ GisPolygon & polygon, RingClipper * ring_clip)
{

    list<Box *> bd_intersection;
    list<Box *> enclosed_intersection;
    bool found =
        intersect(grid, (const Ring &)polygon.front(), bd_intersection, enclosed_intersection);

    if (found)
    {
        //handle the enclosed boxes
        m_intersection.insert(m_intersection.end(), enclosed_intersection.begin(),
                              enclosed_intersection.end());

        //process intersection information
        if (this->m_b_save_intersection_info)
        {
            //store intersection info here...
            for (auto box : enclosed_intersection)
            {
                //store intersection info here...
                vector<PolygonCellIntersectionInfo> infos;
                build_intersection_info(box, polygon, infos, ring_clip);
                m_polygon_intersection_info.insert(m_polygon_intersection_info.end(), infos.begin(),
                                                   infos.end());
            }
        }

        //handle the boundary boxes if needed
        if (m_b_include_poly_boundary)
        {
            m_intersection.insert(m_intersection.end(), bd_intersection.begin(),
                                  bd_intersection.end());
            //process information
            if (this->m_b_save_intersection_info)
            {
                //store intersection info here...
                for (auto & i : bd_intersection)
                {

                    //store intersection info here...
                    vector<PolygonCellIntersectionInfo> infos;
                    build_intersection_info(i, polygon, infos, ring_clip);
                    m_polygon_intersection_info.insert(m_polygon_intersection_info.end(),
                                                       infos.begin(), infos.end());
                }
            }
        } //end if(m_b_include_poly_boundary)
    }

    return found;
}

} //end of namespace
