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

# gridgen libraries shared by the executable and the unit tests
gridgen_dep = declare_dependency(
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
                                )

# version header; vcs_tag adds the git commit to a development version
# (containing 'dev'), and a release version gets the release date
version = meson.project_version()
version_conf = configuration_data()
version_conf.set('VERSION', version)
version_conf.set('DATE', version.contains('dev') ? '' : release_date)
version_conf.set('VCS_TAG', '@VCS_TAG@')
version_in = configure_file(input: 'gridgen_version.h.in',
                            output: 'gridgen_version.h.in',
                            configuration: version_conf)
version_h = vcs_tag(command: [python, meson.project_source_root() / 'scripts' / 'vcs_tag_suffix.py', version],
                    input: version_in,
                    output: 'gridgen_version.h',
                    replace_string: '@VCS_TAG@',
                    fallback: '')

# on Linux the executable carries the C++ library rather than needing that of
# the computer it runs on; libgcc stays shared, since linking it statically
# raised the glibc requirement from 2.34 to 2.35 with GCC 13, and a sanitized
# build keeps the shared C++ library its runtime expects
gridgen_link_args = []
if host_machine.system() == 'linux' and get_option('b_sanitize') == 'none'
  gridgen_link_args = cc.get_supported_link_arguments('-static-libstdc++')
endif

exe = executable(
                 meson.project_name(),
                 ['gridgenusg.cpp', version_h],
                 dependencies: gridgen_dep,
                 link_args: gridgen_link_args,
                 install: true
                )

