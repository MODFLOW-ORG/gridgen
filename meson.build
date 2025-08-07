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
