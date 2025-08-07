project(
  'gridgen',
  ['cpp'],
  version: '1.0.02',
  license: 'CC0',
  meson_version: '>= 1.1.0',
  default_options : [
    'b_vscrt=static_from_buildtype', # Link runtime libraries statically on Windows
    'buildtype=release',
])

cc = meson.get_compiler('cpp')
cc_id = cc.get_id()
message('Compiler ID:', cc_id)

subdir('src')

# tests
test('build the quadtree grid grid02qtg', 
     find_program('gridgen'),
     args : ['buildqtg', 'action01_buildqtg.dfn'],
     workdir: meson.current_source_dir()/'test')

test('write modflow-usg data', 
     find_program('gridgen'),
     args : ['grid02qtg-to-usgdata', 'action02_writeusgdata.dfn'],
     workdir: meson.current_source_dir()/'test')

test('create shapefiles for both grids - 1', 
     find_program('gridgen'),
     args : ['grid01mfg-to-polyshapefile', 'action03_shapefile.dfn'],
     workdir: meson.current_source_dir()/'test')

test('create shapefiles for both grids - 2', 
     find_program('gridgen'),
     args : ['grid02qtg-to-polyshapefile', 'action03_shapefile.dfn'],
     workdir: meson.current_source_dir()/'test')

test('create shapefiles for both grids - 3', 
     find_program('gridgen'),
     args : ['grid01mfg-to-pointshapefile', 'action03_shapefile.dfn'],
     workdir: meson.current_source_dir()/'test')

test('create shapefiles for both grids - 4', 
     find_program('gridgen'),
     args : ['grid02qtg-to-pointshapefile', 'action03_shapefile.dfn'],
     workdir: meson.current_source_dir()/'test')

test('intersect features with quadtree grid - canal', 
     find_program('gridgen'),
     args : ['canal_grid02qtg_lay1_intersect', 'action04_intersect.dfn'],
     workdir: meson.current_source_dir()/'test')

test('intersect features with quadtree grid - chd', 
     find_program('gridgen'),
     args : ['chd_grid02qtg_lay1_intersect', 'action04_intersect.dfn'],
     workdir: meson.current_source_dir()/'test')

test('create a vtk file of the grid - 1', 
     find_program('gridgen'),
     args : ['grid01mfg-to-vtkfile', 'action05_vtkfile.dfn'],
     workdir: meson.current_source_dir()/'test')

test('create a vtk file of the grid - 2', 
     find_program('gridgen'),
     args : ['grid02qtg-to-vtkfile', 'action05_vtkfile.dfn'],
     workdir: meson.current_source_dir()/'test')

test('create a vtk file of the grid - 3', 
     find_program('gridgen'),
     args : ['grid02qtg-to-vtkfilesv', 'action05_vtkfile.dfn'],
     workdir: meson.current_source_dir()/'test')
