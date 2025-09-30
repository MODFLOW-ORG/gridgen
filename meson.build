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

subdir('src')

# tests
test('initial', 
     exe,
     args : ['buildqtg', 'action01_buildqtg.dfn'],
     workdir: meson.current_source_dir()/'test')

test('write_data', 
     exe,
     args : ['grid02qtg-to-usgdata', 'action02_writeusgdata.dfn'],
     workdir: meson.current_source_dir()/'test')

test('shapefileA', 
     exe,
     args : ['grid01mfg-to-polyshapefile', 'action03_shapefile.dfn'],
     workdir: meson.current_source_dir()/'test')

test('shapefileB', 
     exe,
     args : ['grid02qtg-to-polyshapefile', 'action03_shapefile.dfn'],
     workdir: meson.current_source_dir()/'test')

test('shapefileC', 
     exe,
     args : ['grid01mfg-to-pointshapefile', 'action03_shapefile.dfn'],
     workdir: meson.current_source_dir()/'test')

test('shapefileD', 
     exe,
     args : ['grid02qtg-to-pointshapefile', 'action03_shapefile.dfn'],
     workdir: meson.current_source_dir()/'test')

test('intersectA', 
     exe,
     args : ['canal_grid02qtg_lay1_intersect', 'action04_intersect.dfn'],
     workdir: meson.current_source_dir()/'test')

test('intersectB', 
     exe,
     args : ['chd_grid02qtg_lay1_intersect', 'action04_intersect.dfn'],
     workdir: meson.current_source_dir()/'test')

test('vtkA', 
     exe,
     args : ['grid01mfg-to-vtkfile', 'action05_vtkfile.dfn'],
     workdir: meson.current_source_dir()/'test')

test('vtkB', 
     exe,
     args : ['grid02qtg-to-vtkfile', 'action05_vtkfile.dfn'],
     workdir: meson.current_source_dir()/'test')

test('vtkC', 
     exe,
     args : ['grid02qtg-to-vtkfilesv', 'action05_vtkfile.dfn'],
     workdir: meson.current_source_dir()/'test')
