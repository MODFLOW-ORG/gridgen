grid_sources = files(
    'Box.cpp',
    'ModflowGrid.cpp',
    'Octree.cpp',
    'Cube.cpp',
    'ModflowGrid2D.cpp',
    'QuadTree3D.cpp',
    'GridDomain.cpp',
    'NodeGroup.cpp',
    'QuadTreeBuilder.cpp',
)

grid_lib = static_library(
                          'grid_lib', 
                          grid_sources,
                          include_directories : [
                                                 inc_basic, 
                                                 inc_shapefile,
                                                 inc_shapelib,
                                                 inc_vpass,
                                                 inc_tools, 
                                                 inc_grid,
                                                 inc_intersection,
                                                 ],
)