basic_sources = files(
    'polygon.cpp',
    'Util.cpp',
    'Basic.cpp',
    'polyline.cpp',
    'Point.h',
)

basic_lib = static_library('basic_lib', 
                           basic_sources,
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
