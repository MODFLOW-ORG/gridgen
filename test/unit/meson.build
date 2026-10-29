unit_sources = files(
    'test_box.cpp',
    'test_def_parser.cpp',
    'test_ghostnode.cpp',
    'test_modflowgrid.cpp',
    'test_quadtree.cpp',
    'test_quadtree_io.cpp',
    'test_vtk.cpp',
)

unit_exe = executable(
                      'unit_tests',
                      unit_sources,
                      dependencies : [gridgen_dep, gtest_dep],
                     )

test('unit', unit_exe, protocol : 'gtest', timeout : 300)
