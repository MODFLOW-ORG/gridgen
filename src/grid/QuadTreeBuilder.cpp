/*! \file QuadTreeBuilder.cpp

\brief QuadTree Builder

\author <a href="masc.cs.gmu.edu/">MASC group</a>, George Mason University,
        <a href="profile.usgs.gov/langevin/">Christian Langevin</a>, USGS
\bug    No known bugs.

 */

//   $Id: $
#include <locale>
#include <string>
#include "QuadTreeBuilder.h"
#include "GridgenError.h"
#include "Warnings.h"
#include <sstream>
#include "shapefile/shpReader.h"
#include "GridDomain.h"

namespace cusg
{

inline string tolower(const string & s)
{
    string lower;
    for (char i : s)
    {
        char c = ::tolower(i);
        lower.push_back(c);
    }
    return lower;
}

inline int getField(string field_name, ShapeReader & shpReader)
{
    int id = 0;
    vector<ShapeIntField> & int_fields = shpReader.getIntFields();
    for (vector<ShapeIntField>::iterator sit = int_fields.begin(); sit != int_fields.end();
         ++sit, ++id)
    {
        if (field_name == sit->name) return id;
    }

    return -1;
}

void Quadtree_Builder::build(QuadTree3D * qtree, quadtree_builder_raw_data & rawdata)
{
    //rotate qtree
    ModflowGrid * rmfgrid = qtree->getModflowGrid();
    double cx, cy, angle;
    rmfgrid->getRotatePara(cx, cy, angle);

    //set z-value interpolation methods
    for (int i = 0; i < rawdata.modflow_data->nlay; i++)
    {
        string operation_top = tolower(rawdata.top[i].first);

        if (operation_top == "replicate")
            qtree->setTopOperators(i, QuadTree3D::Z_VALUE_REPLICATE_OPERATOR);
        else if (operation_top == "interpolate")
            qtree->setTopOperators(i, QuadTree3D::Z_VALUE_LINEAR_INTERPOLATION_OPERATOR);
        else if (operation_top == "asciigrid") //ASCII file
        {
            qtree->setTopOperators(i, QuadTree3D::Z_VALUE_LINEAR_INTERPOLATION_OPERATOR);
            qtree->setTopSource(i, rawdata.top[i].second);

            if (rawdata.top_area_weighted.find(i) != rawdata.top_area_weighted.end())
            {
                qtree->setTopAreaWeighted(i, rawdata.top_area_weighted[i]);
            }
        }
        else
        {
            throw GridgenError("Unknown top " + to_string(i) +
                               " elevation operator: " + rawdata.top[i].first);
        }

        string operation_bot = tolower(rawdata.bot[i].first);
        if (operation_bot == "replicate")
            qtree->setBotOperators(i, QuadTree3D::Z_VALUE_REPLICATE_OPERATOR);
        else if (operation_bot == "interpolate")
            qtree->setBotOperators(i, QuadTree3D::Z_VALUE_LINEAR_INTERPOLATION_OPERATOR);
        else if (operation_bot == "asciigrid") //ASCII file
        {
            qtree->setBotOperators(i, QuadTree3D::Z_VALUE_LINEAR_INTERPOLATION_OPERATOR);
            qtree->setBotSource(i, rawdata.bot[i].second);

            if (rawdata.bot_area_weighted.find(i) != rawdata.bot_area_weighted.end())
            {
                qtree->setBotAreaWeighted(i, rawdata.bot_area_weighted[i]);
            }
        }
        else
        {
            throw GridgenError("Unknown bottom " + to_string(i) + " elevation operator: (" +
                               rawdata.bot[i].first + ")");
        }
    }

    //check the refinement rules
    map<string, ShapeReader> shp_readers;
    for (int i = 0; i < rawdata.modflow_data->nlay; i++)
    {
        cout << "  * Refine grid at layer " << i + 1 << endl;

        list<refinement_raw_data *> & refinements = rawdata.refine_by_features[i];
        for (auto refinement : refinements)
        {
            string block = "REFINEMENT_FEATURES " + refinement->name;
            ShapeReader & shp_reader = shp_readers[refinement->shapefile];
            if (!shp_reader.read(refinement->shapefile))
                throw GridgenError("cannot open shapefile " + refinement->shapefile + " (" + block +
                                   "); its .shp or .dbf file is missing or cannot be read");

            //prepare for attribute based refinement and store the specified using field_index
            int field_index = -1; //index of refinement_level
            bool bByAttribute = (refinement->refinement_level_by_attribute.empty() == false);
            vector<ShapeIntField> & int_fields = shp_reader.getIntFields();

            if (bByAttribute)
            {
                field_index = getField(refinement->refinement_level_by_attribute, shp_reader);
                if (field_index == -1)
                    throw GridgenError("shapefile " + refinement->shapefile + " (" + block +
                                       ") has no integer field named " +
                                       refinement->refinement_level_by_attribute +
                                       " for REFINEMENT_LEVEL_BY_ATTRIBUTE");
            }

            // refine using the given feature (refinement)
            string feature_type = tolower(refinement->featuretype);
            int nfeatures = 0, noutside = 0;
            if (feature_type == "line")
            {
                vector<GIS_plyline> & arcs = shp_reader.getArcList();
                //subdivide the qtree using the features
                nfeatures = arcs.size();
                for (int j = 0; j < nfeatures; j++)
                {
                    //rotate
                    if (angle != 0) arcs[j].rotate(cx, cy, -angle);

                    int level = (bByAttribute) ? int_fields[field_index].values[arcs[j].m_id]
                                               : refinement->refinement_level;
                    if (!qtree->refine((c_plyline &)arcs[j], level, i)) noutside++;

                    //rotate back
                    if (angle != 0) arcs[j].rotate(cx, cy, angle);
                }
            } //end if
            else if (feature_type == "polygon")
            {
                vector<GIS_polygon> & polygons = shp_reader.getPolyList();
                //subdivide the qtree using the features
                nfeatures = polygons.size();

                for (int j = 0; j < nfeatures; j++)
                {
                    //make sure that bounding box is computed
                    polygons[j].buildBoxAndCenter();

                    //rotate
                    if (angle != 0) polygons[j].rotate(cx, cy, -angle);

                    int level = (bByAttribute) ? int_fields[field_index].values[polygons[j].m_id]
                                               : refinement->refinement_level;
                    if (!qtree->refine(polygons[j], level, i)) noutside++; //i is the layer

                    //rotate back
                    if (angle != 0) polygons[j].rotate(cx, cy, angle);
                }
            }
            else if (feature_type == "point")
            {
                vector<GIS_Point2d> & pts = shp_reader.getPointList();
                //subdivide the qtree using the features
                nfeatures = pts.size();
                for (int j = 0; j < nfeatures; j++)
                {
                    //rotate
                    if (angle != 0) pts[j].rotate(cx, cy, -angle);

                    int level = (bByAttribute) ? int_fields[field_index].values[pts[j].m_id]
                                               : refinement->refinement_level;
                    if (!qtree->refine(pts[j], level, i)) noutside++;

                    //rotate back
                    if (angle != 0) pts[j].rotate(cx, cy, angle);
                }
            }
            else
            {
                throw GridgenError("FEATURE_TYPE " + refinement->featuretype + " (" + block +
                                   ") is not line, polygon, or point");
            }

            if (nfeatures == 0)
                warn("Shapefile " + refinement->shapefile + " (" + block + ") has no " +
                     feature_type + " features; the grid is not refined for it.");
            warn_features_outside(noutside, nfeatures, refinement->shapefile, block, i,
                                  "the grid is not refined for them");
        } //end for f

    } //end for i

    if (rawdata.smoothing)
    {
        cout << "  * Smooth grid: horizontal=" << rawdata.smoothing_level_horizontal
             << " vertical=" << rawdata.smoothing_level_vertical << endl;
        qtree->smooth_refinement(rawdata.smoothing_level_horizontal,
                                 rawdata.smoothing_level_vertical);
    }

    //if active domain is specified
    GridDomain gd;
    for (int i = 0; i < rawdata.modflow_data->nlay; i++)
    {
        list<active_domain_raw_data *> & refinements = rawdata.active_domain_data[i];
        if (refinements.empty()) continue;

        cout << "  * Determine the active domain for layer " << i + 1 << endl;

        qtree->nodegroup.disable_nodes_in_layer(i);
        for (auto & f : refinements)
        {
            gd.build(qtree, i, f);
        }
    }

    //
    // ask the grid to renumber nodes
    // node number are used to retrieve the leaf nodoes
    //
    qtree->set_one_based_node_numbering(rawdata.one_based_node_numbering);
    qtree->number_nodes();
    qtree->refresh();

    //
    // build surface elevation
    //
    cout << "  * Build surface elevation" << endl;
    double * top_elev = qtree->get_top();
    double * bot_elev = qtree->get_bot();
    qtree->recompute_z_of_cells();

    if (rawdata.auto_alignment)
    {
        //perform auto alignment here
        qtree->align_elev();
    }

    // count number of active cells with negative thickness
    analyze_elevation(qtree, top_elev, bot_elev);

} // void Quadtree_Builder::build

void Quadtree_Builder::analyze_elevation(QuadTree3D * qtree, double * top_elev, double * bot_elev)
{
    int nodesize = qtree->nodes();
    int negative_count = 0;
    int max_negative_node = -1;
    double max_negative_thickness = -FLT_MAX;

    //for each leaf node
    for (int nodenumber = 0; nodenumber < nodesize; nodenumber++)
    {
        Box * box = qtree->get_nodeobj(nodenumber);
        assert(box->isLeaf);                //make sure this node is a leaf node...
        if (box->active == false) continue; //disabled...

        //check to top and bottom
        if (top_elev[nodenumber] < bot_elev[nodenumber]) //negative thickness..
        {
            double negative_thickness = bot_elev[nodenumber] - top_elev[nodenumber];
            negative_count++;
            if (negative_thickness > max_negative_thickness)
                max_negative_thickness = negative_thickness;
            max_negative_node = nodenumber;
        }
    } //end for nodenumber

    if (negative_count > 0)
    {
        stringstream msg;
        msg << "There are " << negative_count << " cells with negative thickness. Node #"
            << max_negative_node
            << " has the largest negative thickness=" << -max_negative_thickness;
        warn(msg.str());
    }
    //
}
//end Quadtree_Builder::analyze_elevation()

} //namespace cusg
