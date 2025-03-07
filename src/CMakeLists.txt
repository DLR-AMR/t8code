if ( ${T8CODE_BUILD_AS_SHARED_LIBRARY} )
    add_library( T8 SHARED )
    set_target_properties( T8 PROPERTIES VERSION ${CMAKE_PROJECT_VERSION} SOVERSION ${T8CODE_VERSION_RAW})
    set_target_properties( T8 PROPERTIES POSITION_INDEPENDENT_CODE ON )
else()
    add_library( T8 STATIC )
endif()

add_library( T8CODE::T8 ALIAS T8 )

set_target_properties( T8 PROPERTIES OUTPUT_NAME t8 )

target_include_directories( T8 PUBLIC
  $<BUILD_INTERFACE:${CMAKE_CURRENT_LIST_DIR}>
  $<INSTALL_INTERFACE:${CMAKE_INSTALL_PREFIX}/include>
)

t8_add_cmake_properties (T8)

target_sources( T8 PRIVATE 
    t8.c 
    t8_eclass.c 
    t8_element_shape.c
    t8_element.cxx
    t8_refcount.c 
    t8_version.c 
    t8_vtk.c 
    t8_cmesh/t8_cmesh.cxx 
    t8_cmesh/t8_cmesh_cad.cxx 
    t8_cmesh/t8_cmesh_triangle.cxx 
    t8_cmesh/t8_cmesh_stash.c 
    t8_cmesh/t8_cmesh_vtk_reader.cxx 
    t8_cmesh/t8_cmesh_save.cxx 
    t8_cmesh/t8_cmesh_trees.cxx 
    t8_cmesh/t8_cmesh_commit.cxx 
    t8_cmesh/t8_cmesh_partition.cxx 
    t8_cmesh/t8_cmesh_copy.c 
    t8_cmesh/t8_cmesh_geometry.cxx 
    t8_cmesh/t8_cmesh_examples.cxx 
    t8_cmesh/t8_cmesh_helpers.cxx 
    t8_cmesh/t8_cmesh_offset.c 
    t8_cmesh/t8_cmesh_readmshfile.cxx 
    t8_data/t8_shmem.c 
    t8_data/t8_containers.cxx 
    t8_forest/t8_forest_adapt.cxx 
    t8_forest/t8_forest_partition.cxx 
    t8_forest/t8_forest.cxx 
    t8_forest/t8_forest_private.c 
    t8_forest/t8_forest_ghost.cxx 
    t8_forest/t8_forest_iterate.cxx 
    t8_forest/t8_forest_balance.cxx 
    t8_forest/t8_forest_search/t8_forest_search.cxx 
    t8_geometry/t8_geometry.cxx 
    t8_geometry/t8_geometry_helpers.c 
    t8_geometry/t8_geometry_base.cxx 
    t8_geometry/t8_geometry_handler.cxx 
    t8_geometry/t8_geometry_with_vertices.cxx 
    t8_geometry/t8_geometry_implementations/t8_geometry_analytic.cxx 
    t8_geometry/t8_geometry_implementations/t8_geometry_cad.cxx 
    t8_geometry/t8_geometry_implementations/t8_geometry_linear.cxx 
    t8_geometry/t8_geometry_implementations/t8_geometry_linear_axis_aligned.cxx 
    t8_geometry/t8_geometry_implementations/t8_geometry_lagrange.cxx 
    t8_geometry/t8_geometry_implementations/t8_geometry_zero.cxx 
    t8_geometry/t8_geometry_implementations/t8_geometry_examples.cxx 
    t8_schemes/t8_scheme.cxx
    t8_schemes/t8_default/t8_default.cxx
    t8_schemes/t8_default/t8_default_hex/t8_default_hex.cxx
    t8_schemes/t8_default/t8_default_hex/t8_dhex_bits.c
    t8_schemes/t8_default/t8_default_line/t8_default_line.cxx
    t8_schemes/t8_default/t8_default_line/t8_dline_bits.c
    t8_schemes/t8_default/t8_default_prism/t8_default_prism.cxx
    t8_schemes/t8_default/t8_default_prism/t8_dprism_bits.c
    t8_schemes/t8_default/t8_default_pyramid/t8_default_pyramid.cxx
    t8_schemes/t8_default/t8_default_pyramid/t8_dpyramid_bits.c
    t8_schemes/t8_default/t8_default_pyramid/t8_dpyramid_connectivity.c 
    t8_schemes/t8_default/t8_default_quad/t8_default_quad.cxx
    t8_schemes/t8_default/t8_default_quad/t8_default_quad_bits.cxx
    t8_schemes/t8_default/t8_default_tet/t8_default_tet.cxx
    t8_schemes/t8_default/t8_default_tet/t8_dtet_bits.c
    t8_schemes/t8_default/t8_default_tet/t8_dtet_connectivity.c
    t8_schemes/t8_default/t8_default_tri/t8_default_tri.cxx
    t8_schemes/t8_default/t8_default_tri/t8_dtri_bits.c
    t8_schemes/t8_default/t8_default_tri/t8_dtri_connectivity.c
    t8_schemes/t8_default/t8_default_vertex/t8_default_vertex.cxx
    t8_types/t8_vec.cxx
    t8_schemes/t8_standalone/t8_standalone.cxx
    t8_vtk/t8_vtk_polydata.cxx 
    t8_vtk/t8_vtk_unstructured.cxx 
    t8_vtk/t8_vtk_parallel.cxx 
    t8_vtk/t8_vtk_reader.cxx 
    t8_vtk/t8_vtk_writer.cxx
    t8_vtk/t8_vtk_write_ASCII.cxx
    t8_vtk/t8_vtk_writer_helper.cxx
)
if (T8CODE_ENABLE_NETCDF)
  target_sources( T8 PRIVATE
    t8_netcdf.c
    t8_cmesh/t8_cmesh_netcdf.c
    t8_forest/t8_forest_netcdf.cxx
  )
endif(T8CODE_ENABLE_NETCDF)

get_target_property( T8_LIBS_LIST T8 LINK_LIBRARIES )
string( REPLACE ";" " " T8_LIBS "${T8_LIBS_LIST}" )
target_compile_definitions( T8 PUBLIC T8_LIBS="${T8_LIBS}" )

install( FILES
    t8.h
    t8_cmesh.h
    t8_cmesh.hxx
    t8_cmesh_readmshfile.h
    t8_cmesh_tetgen.h
    t8_cmesh_triangle.h
    t8_cmesh_vtk_reader.hxx
    t8_eclass.h
    t8_element.h
    t8_element_shape.h
    t8_mat.h
    t8_refcount.h
    t8_types/t8_vec.hxx
    t8_version.h
    t8_vtk.h
    t8_windows.h DESTINATION include
)
if (T8CODE_ENABLE_NETCDF)
  install( FILES
    t8_cmesh_netcdf.h
    t8_forest_netcdf.h
    t8_netcdf.h
    DESTINATION include
  )
endif(T8CODE_ENABLE_NETCDF)

install( DIRECTORY t8_cmesh DESTINATION include FILES_MATCHING PATTERN "*.h" )
install( DIRECTORY t8_data DESTINATION include FILES_MATCHING PATTERN "*.h" )
install( DIRECTORY t8_forest DESTINATION include FILES_MATCHING
  PATTERN "*.h"
  PATTERN "*private.h" EXCLUDE )
install( DIRECTORY t8_geometry DESTINATION include FILES_MATCHING PATTERN "*.h" )
install( DIRECTORY t8_schemes DESTINATION include FILES_MATCHING PATTERN "*.h" )
install( DIRECTORY t8_vtk DESTINATION include FILES_MATCHING PATTERN "*.h" )

install( DIRECTORY t8_cmesh DESTINATION include FILES_MATCHING PATTERN "*.hxx" )
install( DIRECTORY t8_data DESTINATION include FILES_MATCHING PATTERN "*.hxx" )
install( DIRECTORY t8_forest DESTINATION include FILES_MATCHING PATTERN "*.hxx" )
install( DIRECTORY t8_geometry DESTINATION include FILES_MATCHING PATTERN "*.hxx" )
install( DIRECTORY t8_schemes DESTINATION include FILES_MATCHING PATTERN "*.hxx" )
install( DIRECTORY t8_vtk DESTINATION include FILES_MATCHING PATTERN "*.hxx" )
install( DIRECTORY t8_types DESTINATION include FILES_MATCHING PATTERN "*.hxx" )

install( TARGETS T8 DESTINATION lib )
install( TARGETS T8 EXPORT ${PROJECT_NAME}-targets )

include( CMakePackageConfigHelpers )

configure_package_config_file( ${CMAKE_CURRENT_LIST_DIR}/config.cmake.in
  ${CMAKE_CURRENT_BINARY_DIR}/CMakeFiles/${PROJECT_NAME}Config.cmake
INSTALL_DESTINATION cmake
)

write_basic_package_version_file(
${CMAKE_CURRENT_BINARY_DIR}/CMakeFiles/${PROJECT_NAME}ConfigVersion.cmake
COMPATIBILITY SameMajorVersion
VERSION "${T8CODE_VERSION_MAJOR}.${T8CODE_VERSION_MINOR}.${T8CODE_VERSION_PATCH}"
)

install( EXPORT ${PROJECT_NAME}-targets
NAMESPACE ${PROJECT_NAME}::
DESTINATION cmake
)

install( FILES
${CMAKE_CURRENT_BINARY_DIR}/CMakeFiles/${PROJECT_NAME}Config.cmake
${CMAKE_CURRENT_BINARY_DIR}/CMakeFiles/${PROJECT_NAME}ConfigVersion.cmake
DESTINATION cmake
)
