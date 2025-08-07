intersection_sources = files(
    'BoxPlyClip.cpp',
    'GridIntersection.cpp',
    'PolygonChecker.cpp',
    'intersection.h',
)

intersection_lib = static_library(
                                  'intersection_lib', 
                                  intersection_sources,
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
