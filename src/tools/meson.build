tools_sources = files(
    'ascii_grid.cpp',
    'def_struct.cpp',
    'GridToVTK.cpp',
    'def_exporter.cpp',
    'GhostNode.cpp',
    'reader_treestruct.cpp',
    'def_parser.cpp',
    'GridToShp.cpp',
    'shapefile' / 'intersecReadWrite.cpp',
    'shapefile' / 'shpGridToLine.cpp',
    'shapefile' / 'shpReader.cpp',
    'shapefile' / 'shapelib' / 'dbfopen.cpp',
    'shapefile' / 'shapelib' / 'safileio.cpp',
    'shapefile' / 'shapelib' / 'shpopen.cpp',
    'VPass' / 'vert_pass.cpp',
)

tools_lib = static_library(
                           'tools_lib', 
                            tools_sources,
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
