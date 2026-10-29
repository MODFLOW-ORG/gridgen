project(
  'gridgen',
  ['cpp'],
  version: '1.0.02',
  license: 'CC0',
  meson_version: '>= 1.3.1',
  default_options : [
    'b_vscrt=static_from_buildtype', # Link runtime libraries statically on Windows
    'cpp_std=c++14',
    'buildtype=release',
])

cc = meson.get_compiler('cpp')
cc_id = cc.get_id()
message('Compiler ID:', cc_id)

compile_args = []
link_args = []
if cc_id == 'intel'
  compile_args += ['-frtti', '-diag-disable=10441']
  link_args += ['-diag-disable=10441']
endif

add_project_arguments(cc.get_supported_arguments(compile_args), language: 'cpp')
add_project_link_arguments(cc.get_supported_arguments(link_args), language: 'cpp')

subdir('src')

# unit tests (GoogleTest 1.17 requires C++17; 1.16 is not in WrapDB)
gtest_dep = dependency(
                       'gtest_main',
                       version: ['>=1.14', '<1.17'],
                       required: get_option('unit_tests'),
                       allow_fallback: true,
                      )
if gtest_dep.found()
  subdir('test/unit')
endif

# regression tests
python = find_program('python3', 'python', required: false)
if python.found()
  test('determinism',
       python,
       args : [meson.current_source_dir()/'test'/'compare_outputs.py',
               '--exe', exe],
       timeout: 600)
endif

# smoke tests share the test directory, so they run serially in order
test('initial', 
     exe,
     args : ['buildqtg', 'action01_buildqtg.dfn'],
     workdir: meson.current_source_dir()/'test',
     is_parallel: false,
     priority: 11)

test('write_data', 
     exe,
     args : ['grid02qtg-to-usgdata', 'action02_writeusgdata.dfn'],
     workdir: meson.current_source_dir()/'test',
     is_parallel: false,
     priority: 10)

test('shapefileA', 
     exe,
     args : ['grid01mfg-to-polyshapefile', 'action03_shapefile.dfn'],
     workdir: meson.current_source_dir()/'test',
     is_parallel: false,
     priority: 9)

test('shapefileB', 
     exe,
     args : ['grid02qtg-to-polyshapefile', 'action03_shapefile.dfn'],
     workdir: meson.current_source_dir()/'test',
     is_parallel: false,
     priority: 8)

test('shapefileC', 
     exe,
     args : ['grid01mfg-to-pointshapefile', 'action03_shapefile.dfn'],
     workdir: meson.current_source_dir()/'test',
     is_parallel: false,
     priority: 7)

test('shapefileD', 
     exe,
     args : ['grid02qtg-to-pointshapefile', 'action03_shapefile.dfn'],
     workdir: meson.current_source_dir()/'test',
     is_parallel: false,
     priority: 6)

test('intersectA', 
     exe,
     args : ['canal_grid02qtg_lay1_intersect', 'action04_intersect.dfn'],
     workdir: meson.current_source_dir()/'test',
     is_parallel: false,
     priority: 5)

test('intersectB', 
     exe,
     args : ['chd_grid02qtg_lay1_intersect', 'action04_intersect.dfn'],
     workdir: meson.current_source_dir()/'test',
     is_parallel: false,
     priority: 4)

test('vtkA', 
     exe,
     args : ['grid01mfg-to-vtkfile', 'action05_vtkfile.dfn'],
     workdir: meson.current_source_dir()/'test',
     is_parallel: false,
     priority: 3)

test('vtkB', 
     exe,
     args : ['grid02qtg-to-vtkfile', 'action05_vtkfile.dfn'],
     workdir: meson.current_source_dir()/'test',
     is_parallel: false,
     priority: 2)

test('vtkC', 
     exe,
     args : ['grid02qtg-to-vtkfilesv', 'action05_vtkfile.dfn'],
     workdir: meson.current_source_dir()/'test',
     is_parallel: false,
     priority: 1)
