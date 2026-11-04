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

void GridDomain::build(Grid * grid, int layer, active_domain_raw_data * data)
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
    {
        cerr << "error: abnormal grid !" << endl;
    }
    double cx, cy, angle;
    rmfgrid->getRotatePara(cx, cy, angle);

    //read from shapefile
    string block = "ACTIVE_DOMAIN " + data->name;
    ShapeReader shp_reader;
    if (!shp_reader.read(data->shapefile))
        throw GridgenError("cannot open shapefile " + data->shapefile + " (" + block + ")");
    int nfeatures = 0, noutside = 0;

    //start to find the domain
    GridIntersection gi(0, layer);
    if (data->feature_type == "line")
    {
        vector<GIS_plyline> & arcs = shp_reader.getArcList();
        //subdivide the qtree using the features
        int arcssize = arcs.size();
        nfeatures = arcssize;
        for (int j = 0; j < arcssize; j++)
        {
            if (angle != 0) arcs[j].rotate(cx, cy, -angle);

            bool touched = false;
            ply_vertex * ptr = arcs[j].getHead();
            do
            {
                ply_vertex * next = ptr->getNext();
                if (next != nullptr)
                {
                    LineSeg2d line(ptr->getPos(), next->getPos());

                    if (qtgrid != nullptr)
                        touched |= gi.intersect(qtgrid, line);
                    else if (mfgrid != nullptr)
                        touched |= gi.intersect(mfgrid, line);
                    else if (mf2dqtgrid != nullptr)
                        touched |= gi.intersect(mf2dqtgrid, line);
                }
                ptr = next;
            } while (ptr != nullptr);
            if (!touched) noutside++;
            if (angle != 0) arcs[j].rotate(cx, cy, angle);
        }
    } //end if
    else if (data->feature_type == "polygon")
    {
        gi.set_include_poly_boundary_flag(data->include_boundary);
        vector<GIS_polygon> & polygons = shp_reader.getPolyList();
        int pgonsize = polygons.size();
        nfeatures = pgonsize;

        for (int j = 0; j < pgonsize; j++)
        {
            //rotate
            if (angle != 0)
            {
                polygons[j].rotate(cx, cy, -angle);
            }

            //make sure that bounding box is computed
            polygons[j].buildBoxAndCenter();

            //determine intersection
            bool touched = false;
            if (qtgrid != nullptr)
                touched = gi.intersect(qtgrid, (c_polygon &)polygons[j]);
            else if (mfgrid != nullptr)
                touched = gi.intersect(mfgrid, polygons[j]);
            else if (mf2dqtgrid != nullptr)
                touched = gi.intersect(mf2dqtgrid, polygons[j]);
            if (!touched) noutside++;

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
        int psize = pts.size();
        nfeatures = psize;
        for (int j = 0; j < psize; j++)
        {
            //rotate
            if (angle != 0) pts[j].rotate(cx, cy, -angle);

            bool touched = false;
            if (qtgrid != nullptr)
                touched = gi.intersect(qtgrid, pts[j]);
            else if (mfgrid != nullptr)
                touched = gi.intersect(mfgrid, pts[j]);
            else if (mf2dqtgrid != nullptr)
                touched = gi.intersect(mf2dqtgrid, pts[j]);
            if (!touched) noutside++;

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
