unit_sources = files(
    'test_asciigrid.cpp',
    'test_box.cpp',
    'test_def_parser.cpp',
    'test_griddomain.cpp',
    'test_gridvertices.cpp',
    'test_ghostnode.cpp',
    'test_intersection.cpp',
    'test_modflowgrid.cpp',
    'test_quadtree.cpp',
    'test_quadtree_io.cpp',
    'test_shapereader.cpp',
    'test_vtk.cpp',
)

unit_exe = executable(
                      'unit_tests',
                      unit_sources,
                      dependencies : [gridgen_dep, gtest_dep],
                      # the test inputs, such as shapefiles
                      cpp_args : ['-DGRIDGEN_TEST_DIR="'
                                  + (meson.project_source_root() / 'test').replace('\\', '/')
                                  + '"'],
                     )

test('unit', unit_exe, protocol : 'gtest', timeout : 300)
