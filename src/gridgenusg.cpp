
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
#include "GridgenError.h"

using namespace cusg;

#define GRIDGEN_COPYRIGHT                                                                          \
    "This program is public domain and is released on the condition that neither the U.S. Geological Survey nor the United States Government may be held liable for any damages resulting from their authorized or unauthorized use."

//run the action from the definition file and return the exit code
static int run(const string & action, const string & dfn_file)
{
    DefinitionParser def_parser;
    def_parser.parse(dfn_file);

    DefinitionBlock * data = def_parser.get_block_data(action);

    if (data == nullptr)
        throw GridgenError("definition file " + dfn_file + " has no block named " + action);

    bool initialized = data->initialize(def_parser.get_block_data());
    if (initialized == false) throw GridgenError(data->label() + " cannot be read");
    data->warn_if_not_action();

    Grid * grid = nullptr;

    //try to case the data to various blocks
    GridBlock * grid_data = dynamic_cast<GridBlock *>(data);
    QuadtreeBuilderBlock * qtree_builder_data = dynamic_cast<QuadtreeBuilderBlock *>(data);
    GridToShapefileBlock * grid2shp = dynamic_cast<GridToShapefileBlock *>(data);
    GridToUsgdataBlock * grid2usg = dynamic_cast<GridToUsgdataBlock *>(data);
    GridIntersectionBlock * grid_intersection = dynamic_cast<GridIntersectionBlock *>(data);
    GridToVtkfileBlock * grid2vtk = dynamic_cast<GridToVtkfileBlock *>(data);

    //check if this is a type of grid
    if (grid_data != nullptr)
    {
        grid = grid_data->build();
        if (grid == nullptr) throw GridgenError(data->label() + " cannot be built");
    }
    else if (grid2shp != nullptr)
    {
        grid = grid2shp->grid;
    }
    else if (grid2usg != nullptr)
    {
        grid = grid2usg->grid;
    }
    else if (grid_intersection != nullptr)
    {
        grid = grid_intersection->grid;
    }
    else if (grid2vtk != nullptr)
    {
        grid = grid2vtk->grid;
    }

    QuadTree3D * qtree_grid = dynamic_cast<QuadTree3D *>(grid);
    ModflowGrid * modflow_grid = dynamic_cast<ModflowGrid *>(grid);

    //check if save files
    if (qtree_builder_data != nullptr)
    {
        DefinitionExporter def_exporter;
        const string & output_file = qtree_builder_data->grid_definition_file;

        def_exporter.save(output_file, *qtree_builder_data);
    }
    else if (grid2shp != nullptr)
    {
        cout << "- Saving grid: " << grid2shp->grid_name << " to shapefile: " << grid2shp->shapefile
             << endl;
        if (qtree_grid != nullptr)
        {
            write_shapefile(qtree_grid, grid2shp->shapefile, grid2shp->feature_type,
                            grid2shp->close_polygons);
        }
        else if (modflow_grid != nullptr)
        {
            write_shapefile(modflow_grid, grid2shp->shapefile, grid2shp->feature_type,
                            grid2shp->without_inactive, grid2shp->one_based_node_numbering,
                            grid2shp->close_polygons);
        }
        else
        {
            throw GridgenError("the grid to export is neither a MODFLOW grid nor a quadtree grid");
        }
    }
    else if (grid2vtk != nullptr)
    {
        cout << "- Saving grid: " << grid2vtk->grid_name << " to vtkfile: " << grid2vtk->vtkfile
             << endl;
        if (qtree_grid != nullptr)
        {
            write_vtk(qtree_grid, grid2vtk->vtkfile, grid2vtk->b_share_vertex);
        }
        else if (modflow_grid != nullptr)
        {
            write_vtk(modflow_grid, grid2vtk->vtkfile, grid2vtk->b_share_vertex,
                      grid2vtk->without_inactive, grid2vtk->one_based_node_numbering);
        }
        else
        {
            throw GridgenError("the grid to export is neither a MODFLOW grid nor a quadtree grid");
        }
    }
    else if (grid2usg != nullptr)
    {
        DefinitionExporter def_exporter;
        def_exporter.save(*grid2usg);
        def_exporter.save_node_coordinate(grid2usg->usg_data_prefix, grid2usg->grid);
    }
    else if (grid_intersection != nullptr)
    {
        GridIntersection gi(0, grid_intersection->layer);
        gi.set_save_intersection_info(true);
        gi.build(*grid_intersection);
        //save results
        gi.save(grid_intersection->output_file, grid_intersection->attributes);
        gi.save_shp(grid_intersection->output_shapefile, grid_intersection->attributes);
    }

    return 0;
}

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

    try
    {
        return run(action, dfn_file);
    }
    catch (const GridgenError & e)
    {
        cerr << "! Error: " << e.what() << endl;
        return 1;
    }
    catch (const exception & e)
    {
        cerr << "! Error: " << e.what() << endl;
        return 1;
    }
}
