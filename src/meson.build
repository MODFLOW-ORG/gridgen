inc_basic = include_directories('basic')
inc_tools = include_directories('tools')
inc_shapelib = include_directories('tools/shapefile/shapelib')
inc_shapefile = include_directories('tools/shapefile')
inc_vpass = include_directories('tools/VPass')
inc_grid = include_directories('grid')
inc_intersection = include_directories('intersection')

subdir('basic')
subdir('tools')
subdir('grid')
subdir('intersection')


exe = executable(
                 meson.project_name(), 
                 'gridgenusg.cpp', 
                 include_directories : [
                                        inc_basic, 
                                        inc_shapefile,
                                        inc_shapelib,
                                        inc_vpass,
                                        inc_tools, 
                                        inc_grid,
                                        inc_intersection,
                                       ],
                 link_with: [basic_lib, tools_lib, grid_lib, intersection_lib],
                 link_language : 'cpp',  
                 install: true
                )

