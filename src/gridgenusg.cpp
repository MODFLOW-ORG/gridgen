
#include "def_struct.h"
#include "def_parser.h"
#include "def_exporter.h"
#include "ModflowGrid.h"
#include "QuadTree3D.h"
#include "shapefile/shpReader.h"

#include "GridIntersection.h"
#include "GridToShp.h"
#include "GridToVTK.h"

#include "vert_pass.h"
#include "GridToVTK.h"
#include "gridgen_version.h"

using namespace cusg;

#define GRIDGEN_COPYRIGHT                                                                          \
    "This program is public domain and is released on the condition that neither the U.S. Geological Survey nor the United States Government may be held liable for any damages resulting from their authorized or unauthorized use."

int main(int argc, char ** argv)
{
    cout << "GRIDGEN Version " << GRIDGEN_VERSION;
    if (strlen(GRIDGEN_DATE) > 0) cout << " " << GRIDGEN_DATE;
    cout << "\n";
    cout << "A program for generating unstructured grids.\n";
    cout << "\n" << GRIDGEN_COPYRIGHT << "\n" << endl;

    if (argc < 3)
    {
        cerr << "usage: " << argv[0] << " action *.dfn" << endl;
        return 1;
    }

    string action;
    string dfn_file;

    for (int i = 1; i < argc; i++)
    {
        // -g opened the OpenGL viewer, which is no longer built
        if (strcmp(argv[i], "-g") == 0)
            continue;
        else if (action.empty())
            action = argv[i];
        else
            dfn_file = argv[i];
    }

    DefinitionParser def_parser;
    def_parser.parse(dfn_file);

    block_raw_data * data = def_parser.getBlockData(action);

    if (data == NULL)
    {
        cerr << "! Error: There is no block named: " << action << endl;
        return 1;
    }

    bool rinit = data->initialize(def_parser.getBlockData());
    if (rinit == false)
    {
        cerr << "! Error: Failed to initialize block: " << action << endl;
        return 1;
    }

    Grid * grid = NULL;

    //try to case the data to various blocks
    grid_raw_data * grid_data = dynamic_cast<grid_raw_data *>(data);
    quadtree_builder_raw_data * qtree_builder_data =
        dynamic_cast<quadtree_builder_raw_data *>(data);
    gridtoshape_raw_data * grid2shp = dynamic_cast<gridtoshape_raw_data *>(data);
    grid_to_usgdata_raw_data * grid2usg = dynamic_cast<grid_to_usgdata_raw_data *>(data);
    grid_intersection_raw_data * gridintersection =
        dynamic_cast<grid_intersection_raw_data *>(data);
    gridtovtk_raw_data * grid2vtk = dynamic_cast<gridtovtk_raw_data *>(data);

    //check if this is a type of grid
    if (grid_data != NULL)
    {
        grid = grid_data->build();
        if (grid == NULL)
        {
            cerr << "! Error: Failed to build block: " << action << endl;
            return 1;
        }
    }
    else if (grid2shp != NULL)
    {
        grid = grid2shp->grid;
    }
    else if (grid2usg != NULL)
    {
        grid = grid2usg->grid;
    }
    else if (gridintersection != NULL)
    {
        grid = gridintersection->grid;
    }
    else if (grid2vtk != NULL)
    {
        grid = grid2vtk->grid;
    }

    QuadTree3D * qtree_grid = dynamic_cast<QuadTree3D *>(grid);
    ModflowGrid * modflow_grid = dynamic_cast<ModflowGrid *>(grid);

    //check if save files
    if (qtree_builder_data != NULL)
    {
        DefinitionExporter def_exporter;
        const string & outputfile = qtree_builder_data->grid_definition_file;

        def_exporter.save(outputfile, *qtree_builder_data);
    }
    else if (grid2shp != NULL)
    {
        cout << "- Saving grid: " << grid2shp->grid_name << " to shapefile: " << grid2shp->shapefile
             << endl;
        //save to shapefile...
        if (qtree_grid != NULL)
        {
            writeQuadtree2Shp(qtree_grid, grid2shp->shapefile, grid2shp->feature_type,
                              grid2shp->concide);
        }
        else if (modflow_grid != NULL)
        {
            writeModflowGrid2Shp(modflow_grid, grid2shp->shapefile, grid2shp->feature_type,
                                 grid2shp->without_inactive, grid2shp->one_based_node_numbering,
                                 grid2shp->concide);
        }
        else
        {
            cerr << "! Error: Unknown Grid type" << endl;
        }
    }
    else if (grid2vtk != NULL)
    {
        cout << "- Saving grid: " << grid2vtk->grid_name << " to vtkfile: " << grid2vtk->vtkfile
             << endl;
        //save to shapefile...
        if (qtree_grid != NULL)
        {
            if (grid2vtk->b_share_vertex)
            {
                writeQuadtreeGrid2VTK_shareV(qtree_grid, grid2vtk->vtkfile);
            }
            else
            {
                writeQuadtreeGrid2VTK(qtree_grid, grid2vtk->vtkfile);
            }
        }
        else if (modflow_grid != NULL)
        {
            if (grid2vtk->b_share_vertex)
            {
                writeModflowGrid2VTK_shareV(modflow_grid, grid2vtk->vtkfile,
                                            grid2vtk->one_based_node_numbering);
            }
            else
            {
                writeModflowGrid2VTK(modflow_grid, grid2vtk->vtkfile,
                                     grid2vtk->one_based_node_numbering);
            }
        }
        else
        {
            cerr << "! Error: Unknown Grid type" << endl;
        }
    }
    else if (grid2usg != NULL)
    {
        DefinitionExporter def_exporter;
        def_exporter.save(*grid2usg);
        def_exporter.save_node_coordinate(grid2usg->usg_data_prefix, grid2usg->grid);
    }
    else if (gridintersection != NULL)
    {
        GridIntersection gi(0, gridintersection->layer);
        gi.set_save_intersection_info(true);
        gi.build(*gridintersection);
        //save results
        gi.save(gridintersection->output_file, gridintersection->attributes);
        gi.saveShp(gridintersection->output_shapefile, gridintersection->attributes);
    }

    return 0;
}
