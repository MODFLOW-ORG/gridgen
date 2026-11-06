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
#include "BoxPlyClip.h"
#include "PolygonChecker.h"
#include "GridgenError.h"
#include "Warnings.h"
using namespace std;

namespace cusg
{

//GridIntersection

void GridIntersection::build(grid_intersection_raw_data & data)
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
        if (qtgrid != nullptr) rmfgrid = qtgrid->getModflowGrid();
    }
    if (rmfgrid == nullptr)
    {
        cerr << "error: abnormal grid !" << endl;
    }
    double cx, cy, angle;
    rmfgrid->getRotatePara(cx, cy, angle);

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
        vector<GIS_plyline> & arcs = m_shp_reader.getArcList();
        //subdivide the qtree using the features
        nfeatures = arcs.size();
        for (int j = 0; j < nfeatures; j++)
        {
            int previousSegInterNum = m_seg_interaction_info.size();
            //intersect returns true even for a segment outside the grid, so count the cells found
            size_t ncells = m_intersection.size();

            int seg_size = arcs[j].getSize() - 1; //getSize returns the number of vertices
            for (int is = 0; is < seg_size; is++)
            {
                GIS_LineSeg2d line = arcs[j].getLineSegment(is);
                //rotate
                if (angle != 0) line.rotate(cx, cy, -angle);

                if (qtgrid != nullptr)
                    intersect(qtgrid, line);
                else if (mfgrid != nullptr)
                    intersect(mfgrid, line);
                else if (mf2dqtgrid != nullptr)
                    intersect(mf2dqtgrid, line);
            } //end for i
            if (m_intersection.size() == ncells) noutside++;

            //merge the intersection segments that are inside the same box
            int afterSegInterNum = m_seg_interaction_info.size();
            int k = previousSegInterNum;
            while (k < afterSegInterNum)
            {
                PolylineCellIntersectionInfo plyInfo;
                plyInfo.box = m_seg_interaction_info[k].box;
                plyInfo.m_curveID = m_seg_interaction_info[k].m_curveID;
                plyInfo.m_startGeoDist = m_seg_interaction_info[k].m_startGeoDist;
                plyInfo.m_endGeoDist = m_seg_interaction_info[k].m_endGeoDist;
                plyInfo.m_len = m_seg_interaction_info[k].m_len;
                plyInfo.m_plyline = new c_plyline();
                plyInfo.m_plyline->beginPoly();

                Point2d & start_point = m_seg_interaction_info[k].m_startPoint;
                Point2d & end_point = m_seg_interaction_info[k].m_endPoint;
                if (angle != 0)
                {
                    start_point.rotate(cx, cy, angle);
                    end_point.rotate(cx, cy, angle); //rotate back
                }
                plyInfo.m_plyline->addVertex(start_point[0], start_point[1]);
                plyInfo.m_plyline->addVertex(end_point[0], end_point[1]);

                while (k < (afterSegInterNum - 1) &&
                       m_seg_interaction_info[k].m_curveID ==
                           m_seg_interaction_info[k + 1].m_curveID &&
                       m_seg_interaction_info[k].box == m_seg_interaction_info[k + 1].box &&
                       m_seg_interaction_info[k].m_endPoint.almost_equ(
                           m_seg_interaction_info[k + 1].m_startPoint))
                {
                    end_point = m_seg_interaction_info[k + 1].m_endPoint;
                    if (angle != 0)
                    {
                        end_point.rotate(cx, cy, angle); //rotate back
                    }
                    plyInfo.m_plyline->addVertex(end_point[0], end_point[1]);
                    plyInfo.m_endGeoDist = m_seg_interaction_info[k + 1].m_endGeoDist;
                    plyInfo.m_len += m_seg_interaction_info[k + 1].m_len;
                    k++;
                }

                plyInfo.m_plyline->endPoly();
                m_plyline_intersection_info.push_back(plyInfo);

                k++;
            }

        } //end for j
    } //end if
    else if (data.feature_type == "polygon")
    {
        int previous_size = m_ply_interaction_info.size();
        vector<GIS_polygon> & polygons = m_shp_reader.getPolyList();
        //subdivide the qtree using the features
        nfeatures = polygons.size();
        for (int j = 0; j < nfeatures; j++)
        {

            if (j % 500 == 0) cout << "- Processing polygon: " << j << "    ......" << endl;

            //rotate
            if (angle != 0) polygons[j].rotate(cx, cy, -angle);

            vector<GIS_polygon *> parts;
            checkHole(polygons[j], parts);

            //loop over all the polygons
            bool intersects = false;
            for (auto & part : parts)
            {
                GIS_polygon polygon_part = *part;

                BoxPlyClip * ply_clip = convertCPolygonToClip(*part);

                if (qtgrid != nullptr)
                    intersects |= intersect(qtgrid, polygon_part, ply_clip); //polygons[j]);
                else if (mfgrid != nullptr)
                    intersects |= intersect(mfgrid, polygon_part, ply_clip); //polygons[j]);
                else if (mf2dqtgrid != nullptr)
                    intersects |= intersect(mf2dqtgrid, polygon_part, ply_clip); //polygons[j]);

                ply_clip->destroy();
                delete ply_clip;
            }
            //the rings are those of polygons[j], which the reader frees
            for (auto polygon_part : parts) delete polygon_part;
            if (!intersects) noutside++;
        }
        int afterSize = m_ply_interaction_info.size();
        if (angle != 0)
        {
            for (int k = previous_size; k < afterSize; k++)
            {
                m_ply_interaction_info[k].m_interPolygon->rotate(cx, cy, angle);
            }
        }
    }
    else if (data.feature_type == "point")
    {
        int previous_size = m_pt_interaction_info.size();
        vector<GIS_Point2d> & pts = m_shp_reader.getPointList();
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
        int afterSize = this->m_pt_interaction_info.size();
        if (angle != 0)
        {
            for (int k = previous_size; k < afterSize; k++)
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
bool GridIntersection::saveShp(string & shpName, vector<string> & attributes)
{
    if (shpName.size() != 0)
    {

        if (m_pt_interaction_info.size() != 0)
        {
            writePntIntersect(this, shpName, attributes, WRITE_SHP);
        }
        else if (m_seg_interaction_info.size() != 0)
        {
            writeLineIntersect(this, shpName, attributes, WRITE_SHP);
        }
        else if (m_ply_interaction_info.size() != 0)
        {
            writePolyIntersect(this, shpName, attributes, WRITE_SHP);
        }
        else
        {
            warn("No feature intersects the grid, so " + shpName + " is not written.");
            return true;
        }
        cout << "- Save shape data to file: " << shpName << endl;
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
            writePntIntersect(this, filename, attributes, WRITE_TXT);
        }
        else if (m_seg_interaction_info.size() != 0)
        {
            writeLineIntersect(this, filename, attributes, WRITE_TXT);
        }
        else if (m_ply_interaction_info.size() != 0)
        {
            writePolyIntersect(this, filename, attributes, WRITE_TXT);
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
        cerr << "output_file is not defined!" << endl;
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
    Index2d id2d = grid->getIndex(pt);

    if (id2d[0] == INT_MAX || id2d[1] == INT_MAX || id2d[0] == INT_MIN || id2d[1] == INT_MIN)
        return false;

    //get the box and put it in m_intersection
    int id = grid->get_nodeid(id2d[0], id2d[1]);
    Box * box = grid->nodegroup[id];
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

    found = getIntersection(&box, line, clipped_line);
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
    return clip_line((ModflowGrid2D *)grid->getModflowGrid(), line, clipped_line);
}

Point2d GridIntersection::grid_center(ModflowGrid2D * grid)
{
    return Point2d(grid->xoffset + grid->width / 2, grid->yoffset + grid->height / 2);
}

Point2d GridIntersection::grid_center(QuadTree3D * grid)
{
    return grid_center((ModflowGrid2D *)grid->getModflowGrid());
}

void GridIntersection::bracket_to_box_boundary(Box & box, Point2d & out, Point2d & in) const
{
    //make sure that in point is in and out point is out of box
    if (box.in(in[0], in[1]) == false)
    {
        cout << "holy" << endl;
    }

    assert(box.in(in[0], in[1]));
    assert(box.in(out[0], out[1]) == false);
    //
    while (true)
    {
        Point2d mid((in[0] + out[0]) / 2, (in[1] + out[1]) / 2);
        if (Equal(mid.get(), in.get())) break;
        if (Equal(mid.get(), out.get())) break;

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
    Index2d id2d = grid->getIndex(pt);
    if (id2d[0] == INT_MAX || id2d[1] == INT_MAX || id2d[0] == INT_MIN || id2d[1] == INT_MIN)
        return false;

    //get the box and put it in m_intersection
    int id = grid->get_nodeid(this->m_layer, id2d[0], id2d[1]);

    Box * box = grid->nodegroup[id];
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
    ModflowGrid * mf = grid->getModflowGrid();

    Index2d id2d = mf->getIndex(pt);

    if (id2d[0] == INT_MAX || id2d[1] == INT_MAX || id2d[0] == INT_MIN || id2d[1] == INT_MIN)
        return false;

    //get the box and put it in m_intersection
    int id = mf->get_nodeid(m_layer, id2d[0], id2d[1]);
    Box * box = grid->nodegroup[id];

    box = box->find(pt[0], pt[1]); //find a leave that contains pt

    //make sure box is valid
    assert(box);

    //remember
    m_intersection.push_back(box);

    ////shpPntId.insert(make_pair(box, m_cur_id));//add by guilin

    //remember the boxes that are exactly on the boundary of the polygon
    if (box->onBound(pt[0], pt[1])) this->m_onBound.push_back(box);

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
    if (SegSegInt(ul.get(), ur.get(), start, end)) return true;
    if (SegSegInt(ur.get(), lr.get(), start, end)) return true;
    if (SegSegInt(lr.get(), ll.get(), start, end)) return true;
    if (SegSegInt(ll.get(), ul.get(), start, end)) return true;

    //JML: I am not sure if the code below is really necessary.
    if (AreaSign(ul.get(), ur.get(), start) == 0)
        if (Between(ul.get(), ur.get(), start)) return true;
    if (AreaSign(ur.get(), lr.get(), start) == 0)
        if (Between(ur.get(), lr.get(), start)) return true;
    if (AreaSign(lr.get(), ll.get(), start) == 0)
        if (Between(lr.get(), ll.get(), start)) return true;
    if (AreaSign(ll.get(), ul.get(), start) == 0)
        if (Between(ll.get(), ul.get(), start)) return true;

    if (AreaSign(ul.get(), ur.get(), end) == 0)
        if (Between(ul.get(), ur.get(), end)) return true;
    if (AreaSign(ur.get(), lr.get(), end) == 0)
        if (Between(ur.get(), lr.get(), end)) return true;
    if (AreaSign(lr.get(), ll.get(), end) == 0)
        if (Between(lr.get(), ll.get(), end)) return true;
    if (AreaSign(ll.get(), ul.get(), end) == 0)
        if (Between(ll.get(), ul.get(), end)) return true;

    //otherwise
    return false;
}

void GridIntersection::buildIntersectionInfo(Box * box, const GIS_polygon & ply,
                                             vector<PolygonCellIntersectionInfo> & infoResult,
                                             BoxPlyClip * ply_clip)
{
    vector<GIS_polygon *> result_polygons;
    intersect_polygon(box, ply, result_polygons, ply_clip);
    for (auto polygon_part : result_polygons)
    {
        PolygonCellIntersectionInfo info;
        info.m_interPolygon = polygon_part;
        info.box = box;
        info.m_polyID = ply.m_id;
        info.m_area = info.m_interPolygon->getArea();

        infoResult.push_back(info);
    }
}

SegmentCellIntersectionInfo GridIntersection::buildIntersectionInfo(Box * box,
                                                                    const GIS_LineSeg2d & line)
{
    SegmentCellIntersectionInfo info;
    info.box = box;
    info.m_segID = line.m_seg_id;
    info.m_curveID = line.m_curve_id;

    LineSeg2d clipped_segment;
    bool found = getIntersection(box, (LineSeg2d)line, clipped_segment);

    if (found == false)
    {
        cerr
            << "! Warning: GridIntersection::buildIntersectionInfo(Box,const GIS_LineSeg2d): Cannot find intersection information"
            << endl;
        return info;
    }

    info.m_len = (clipped_segment.pt[0] - clipped_segment.pt[1]).norm();
    info.m_startGeoDist = (clipped_segment.pt[0] - line.pt[0]).norm() + line.m_start_geo_dist;
    info.m_endGeoDist = (clipped_segment.pt[1] - line.pt[0]).norm() + line.m_start_geo_dist;

    if (info.m_startGeoDist > info.m_endGeoDist)
    {
        swap(info.m_endGeoDist, info.m_startGeoDist);
        info.m_startPoint = clipped_segment.pt[1];
        info.m_endPoint = clipped_segment.pt[0];
    }
    else
    {
        info.m_startPoint = clipped_segment.pt[0];
        info.m_endPoint = clipped_segment.pt[1];
    }

    return info;
}

PointCellIntersectionInfo GridIntersection::buildIntersectionInfo(Box * box, const GIS_Point2d & pt)
{
    PointCellIntersectionInfo info;
    info.box = box;
    info.m_pointID = pt.m_id;
    info.point = pt;
    return info;
}

bool GridIntersection::getIntersection(Box * b, const LineSeg2d & line, LineSeg2d & clipped_segment)
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

    edge_code[0] = SegSegInt(ul.get(), ur.get(), line.pt[0].get(), line.pt[1].get(), p[0].get());
    edge_code[1] = SegSegInt(ur.get(), lr.get(), line.pt[0].get(), line.pt[1].get(), p[1].get());
    edge_code[2] = SegSegInt(lr.get(), ll.get(), line.pt[0].get(), line.pt[1].get(), p[2].get());
    edge_code[3] = SegSegInt(ll.get(), ul.get(), line.pt[0].get(), line.pt[1].get(), p[3].get());

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

                bool overlaps =
                    Union(seg[0].get(), seg[1].get(), line.pt[0].get(), line.pt[1].get(), u);
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

double GridIntersection::intersect_len(Box * b, const LineSeg2d & line)
{
    //check if the end points are in
    LineSeg2d clipped_segment;
    bool found = getIntersection(b, line, clipped_segment);
    if (found) return (clipped_segment.pt[0] - clipped_segment.pt[1]).norm();
    return 0; //no intersection...
}

//
// check if a point is inside a ply
// ray shooting is slow. Can be improved.
//
bool GridIntersection::intersect(const Point2d & pt, const c_ply & ply)
{

    int max_trial = 200;

    Vector2d offset(ply.getHead()->getPos().get());
    float ray_length = ply.getRadius() * 4;
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
        ply_vertex * ptr = ply.getHead();
        bool degenerate = false;
        int int_count = 0;

        do
        {
            ply_vertex * next = ptr->getNext();
            double x[2];
            Point2d edge_start = ptr->getPos() - offset;
            Point2d edge_end = next->getPos() - offset;

            char code = SegSegInt(near.get(), far.get(), edge_start.get(), edge_end.get(), x);

            if (code == '1')
                int_count++;
            else if (code != '0')
            {
                degenerate = true;
                break;
            }

            ptr = next;
        } while (ptr != ply.getHead());

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
bool GridIntersection::intersect(const Point2d & pt, const c_polygon & pgon)
{

    int max_trial = 200;

    Vector2d offset(pgon.front().getHead()->getPos().get());
    const double * pgon_box = pgon.getBBox();
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

        //loop through each ply
        for (const auto & ply : pgon)
        {
            ply_vertex * ptr = ply.getHead();

            do
            {
                ply_vertex * next = ptr->getNext();
                double x[2];
                Point2d edge_start = ptr->getPos() - offset;
                Point2d edge_end = next->getPos() - offset;

                char code = SegSegInt(near.get(), far.get(), edge_start.get(), edge_end.get(), x);

                if (code == '1')
                    int_count++;
                else if (code != '0')
                {
                    degenerate = true;
                    break;
                }

                ptr = next;
            } while (ptr != ply.getHead());

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

//void GridIntersection::intersect_polygon(Box* b, const GIS_polygon& polygon, vector<GIS_polygon*>& result_polygons)
void GridIntersection::intersect_polygon(Box * b, const GIS_polygon & polygon,
                                         vector<GIS_polygon *> & result_polygons,
                                         BoxPlyClip * initial_clip)
{

    double xmin = min(b->x - b->dx / 2, b->x + b->dx / 2);
    double xmax = max(b->x - b->dx / 2, b->x + b->dx / 2);
    double ymin = min(b->y - b->dy / 2, b->y + b->dy / 2);
    double ymax = max(b->y - b->dy / 2, b->y + b->dy / 2);

    /*****************************************************************/
    BoxPlyClip * ply_clip = initial_clip;
    if (initial_clip == nullptr)
    {
        ply_clip = convertCPolygonToClip(polygon);
    }

    list<BoxPlyClip *> clips;
    ply_clip->clip(xmin, xmax, ymin, ymax, clips);

    for (auto & clip_part : clips)
    {
        GIS_polygon * npolygons = new GIS_polygon();
        convertClipToPolygon(clip_part, npolygons);
        result_polygons.push_back(npolygons);

        clip_part->destroy();
        delete clip_part;
    }

    ply_clip->resetStatus();

    if (initial_clip == nullptr)
    {
        ply_clip->destroy();
        delete ply_clip;
    }
    else
    {
        ply_clip->resetStatus();
    }
}

//the side of box that other shares an edge with, or -1 if they are not neighbors
int GridIntersection::checkNeighbor(Box * box, Box * other)
{
    for (int i = 0; i < 4; i++)
    {
        list<Box *> neighbors;
        box->getNeighbors(neighbors, i);
        if (std::find(neighbors.begin(), neighbors.end(), other) != neighbors.end()) return i;
    }

    return -1;
}

//check box edge collinearly intersect with line
bool GridIntersection::checkEdgeCollinear(Box * box, int dir, const LineSeg2d & line)
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
    if (SegSegInt(edge_start, edge_end, line_start, line_end, p) == 'e') return true;

    return false;
}

///collect all cells in the grid from b
void GridIntersection::floodGrid(Box * start_box, set<Box *> & visited, list<Box *> & flooded)
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
            b->getNeighbors(neighbors[i], i);
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

} //end GridIntersection::floodGrid
}
