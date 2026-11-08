/*! \file GridDomain.cpp

\brief a class for computing the active/inactive domain

\author <a href="http://masc.cs.gmu.edu">MASC group</a>, George Mason University
\author <a href="http://profile.usgs.gov/langevin/">Christian Langevin</a>, USGS
\bug    No known bugs.
*/

#include "GridDomain.h"
#include "shapefile/shpReader.h"
#include "GridIntersection.h"
#include "GridgenError.h"
#include "Warnings.h"

namespace cusg
{

void GridDomain::build(Grid * grid, int layer, ActiveDomainBlock * data)
{
    //check the refinement rules
    map<string, ShapeReader> shp_readers;
    vector<Box *> domain_boxes;

    //cast
    QuadTree3D * qtgrid = dynamic_cast<QuadTree3D *>(grid);
    ModflowGrid * mfgrid = dynamic_cast<ModflowGrid *>(grid);
    ModflowGrid2D * mf2dqtgrid = dynamic_cast<ModflowGrid2D *>(grid);

    //rotate
    ModflowGrid * rmfgrid = mfgrid;
    if (rmfgrid == nullptr)
    {
        if (qtgrid != nullptr) rmfgrid = qtgrid->getModflowGrid();
    }
    if (rmfgrid == nullptr)
        throw GridgenError("the grid of " + data->label() +
                           " is neither a MODFLOW grid nor a quadtree grid");
    double cx, cy, angle;
    rmfgrid->get_rotation(cx, cy, angle);

    //read from shapefile
    string block = "ACTIVE_DOMAIN " + data->name;
    ShapeReader shp_reader;
    if (!shp_reader.read(data->shapefile))
        throw GridgenError("cannot open shapefile " + data->shapefile + " (" + block +
                           "); its .shp or .dbf file is missing or cannot be read");
    int nfeatures = 0, noutside = 0;

    //start to find the domain
    GridIntersection gi(0, layer);
    if (data->feature_type == "line")
    {
        vector<GIS_polyline> & arcs = shp_reader.getArcList();
        //subdivide the qtree using the features
        nfeatures = arcs.size();
        for (int j = 0; j < nfeatures; j++)
        {
            if (angle != 0) arcs[j].rotate(cx, cy, -angle);

            //intersect returns true even for a segment outside the grid, so count the cells found
            int ncells = gi.getIntersectionSize();
            LinkedVertex * ptr = arcs[j].getHead();
            do
            {
                LinkedVertex * next = ptr->getNext();
                if (next != nullptr)
                {
                    LineSeg2d line(ptr->getPos(), next->getPos());

                    if (qtgrid != nullptr)
                        gi.intersect_segment(qtgrid, line);
                    else if (mfgrid != nullptr)
                        gi.intersect_segment(mfgrid, line);
                    else if (mf2dqtgrid != nullptr)
                        gi.intersect_segment(mf2dqtgrid, line);
                }
                ptr = next;
            } while (ptr != nullptr);
            if (gi.getIntersectionSize() == ncells) noutside++;
            if (angle != 0) arcs[j].rotate(cx, cy, angle);
        }
    } //end if
    else if (data->feature_type == "polygon")
    {
        gi.set_include_poly_boundary_flag(data->include_boundary);
        vector<GIS_polygon> & polygons = shp_reader.getPolyList();
        nfeatures = polygons.size();

        for (int j = 0; j < nfeatures; j++)
        {
            //rotate
            if (angle != 0)
            {
                polygons[j].rotate(cx, cy, -angle);
            }

            //make sure that bounding box is computed
            polygons[j].buildBoxAndCenter();

            //determine intersection
            bool intersects = false;
            if (qtgrid != nullptr)
                intersects = gi.intersect(qtgrid, (c_polygon &)polygons[j]);
            else if (mfgrid != nullptr)
                intersects = gi.intersect(mfgrid, polygons[j]);
            else if (mf2dqtgrid != nullptr)
                intersects = gi.intersect(mf2dqtgrid, polygons[j]);
            if (!intersects) noutside++;

            //rotate back
            if (angle != 0)
            {
                polygons[j].rotate(cx, cy, angle);
            }
        }
    }
    else if (data->feature_type == "point")
    {
        vector<GIS_Point2d> & pts = shp_reader.getPointList();

        //subdivide the qtree using the features
        nfeatures = pts.size();
        for (int j = 0; j < nfeatures; j++)
        {
            //rotate
            if (angle != 0) pts[j].rotate(cx, cy, -angle);

            bool intersects = false;
            if (qtgrid != nullptr)
                intersects = gi.intersect(qtgrid, pts[j]);
            else if (mfgrid != nullptr)
                intersects = gi.intersect(mfgrid, pts[j]);
            else if (mf2dqtgrid != nullptr)
                intersects = gi.intersect(mf2dqtgrid, pts[j]);
            if (!intersects) noutside++;

            //rotate back
            if (angle != 0) pts[j].rotate(cx, cy, angle);
        }
    } //end if

    if (nfeatures == 0)
        warn("Shapefile " + data->shapefile + " (" + block + ") has no " + data->feature_type +
             " features; no cells in layer " + to_string(layer + 1) + " are made active by it.");
    warn_features_outside(noutside, nfeatures, data->shapefile, block, layer,
                          "no cells are made active by them");

    vector<Box *> & boxes = gi.getIntersections();
    domain_boxes.insert(domain_boxes.end(), boxes.begin(), boxes.end());

    //now go through the intersections and disable the boxes
    for (auto & box : domain_boxes)
    {
        box->active = true;
    }
}

} //end namespace cusg
