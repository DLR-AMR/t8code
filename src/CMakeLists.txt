if ( ${T8CODE_BUILD_AS_SHARED_LIBRARY} )
    add_library( T8 SHARED )
    set_target_properties( T8 PROPERTIES VERSION ${CMAKE_PROJECT_VERSION} SOVERSION ${T8CODE_VERSION_RAW})
    set_target_properties( T8 PROPERTIES POSITION_INDEPENDENT_CODE ON )
else()
    add_library( T8 STATIC )
endif()

add_library( T8CODE::T8 ALIAS T8 )

if( CMAKE_BUILD_TYPE STREQUAL "Release" )
  set (T8_CXXFLAGS "${CMAKE_CXX_FLAGS} ${CMAKE_CXX_FLAGS_RELEASE}")
  set (T8_CFLAGS "${CMAKE_C_FLAGS} ${CMAKE_C_FLAGS_RELEASE}")
endif()

if( CMAKE_BUILD_TYPE STREQUAL "RelWithDebInfo" )
  set (T8_CXXFLAGS "${CMAKE_CXX_FLAGS} ${CMAKE_CXX_FLAGS_RELWITHDEBINFO}")
  set (T8_CFLAGS "${CMAKE_C_FLAGS} ${CMAKE_C_FLAGS_RELWITHDEBINFO}")
endif()

if( CMAKE_BUILD_TYPE STREQUAL "Debug" )
  set (T8_CXXFLAGS "${CMAKE_CXX_FLAGS} ${CMAKE_CXX_FLAGS_DEBUG}")
  set (T8_CFLAGS "${CMAKE_C_FLAGS} ${CMAKE_C_FLAGS_DEBUG}")

    target_compile_definitions( T8 PUBLIC T8_ENABLE_DEBUG )
endif()

if( T8CODE_ENABLE_NETCDF )
  target_link_libraries( T8 PUBLIC netCDF::netcdf )
  target_compile_definitions(T8 PUBLIC
    T8_WITH_NETCDF
    $<$<AND:$<BOOL:${NETCDF_HAVE_NETCDF_PAR}>,$<BOOL:${T8CODE_ENABLE_MPI}>>:T8_WITH_NETCDF_PAR> )
endif()

set_target_properties( T8 PROPERTIES OUTPUT_NAME t8 )

target_include_directories( T8 PUBLIC
  $<BUILD_INTERFACE:${CMAKE_CURRENT_LIST_DIR}>
  $<INSTALL_INTERFACE:${CMAKE_INSTALL_PREFIX}/include>
)

target_link_libraries( T8 PUBLIC P4EST::P4EST SC::SC )

if ( T8CODE_ENABLE_MPI )
    target_compile_definitions( T8 PUBLIC T8_ENABLE_MPI )
    target_compile_definitions( T8 PUBLIC T8_ENABLE_MPIIO )
    target_link_libraries( T8 PUBLIC MPI::MPI_C )
endif()

if( T8CODE_ENABLE_VTK )
    target_compile_definitions( T8 PUBLIC T8_VTK_VERSION_USED="${VTK_MAJOR_VERSION}.${VTK_MINOR_VERSION}" )
    target_compile_definitions( T8 PUBLIC T8_WITH_VTK=1 )
    target_include_directories( T8 PUBLIC ${VTK_INCLUDE_DIR} )
    target_link_libraries( T8 PUBLIC ${VTK_LIBRARIES} )
endif()

if( T8CODE_ENABLE_OCC )
    target_compile_definitions( T8 PUBLIC T8_WITH_OCC=1 )
    target_include_directories( T8 PUBLIC ${OpenCASCADE_INCLUDE_DIR} )
    target_link_libraries( T8 PUBLIC ${OpenCASCADE_LIBRARIES} )
endif()

if( T8CODE_ENABLE_LESS_TESTS )
    target_compile_definitions( T8 PUBLIC T8_ENABLE_LESS_TESTS=1 )
endif()

if( T8CODE_BUILD_PEDANTIC )
    # TODO: Add `-Wextra -Wpedantic`. As of now, t8code would not compile with this.
    # See also: https://github.com/DLR-AMR/t8code/issues/1155.
    target_compile_options( T8 PUBLIC -Wall -pedantic -Werror )
    set (T8_CXXFLAGS "${T8_CXXFLAGS} -W -Wall -Wpedantic -Werror")
    set (T8_CFLAGS "${T8_CFLAGS} -W -Wall -Wpedantic -Werror")
endif()

target_sources( T8 PRIVATE 
    t8.c 
    t8_eclass.c 
    t8_element_shape.c
    t8_element.cxx 
    t8_element_c_interface.cxx 
    t8_mesh.c
    t8_netcdf.c 
    t8_refcount.c 
    t8_version.c 
    t8_vtk.c 
    t8_cmesh/t8_cmesh.cxx 
    t8_cmesh/t8_cmesh_cad.cxx 
    t8_cmesh/t8_cmesh_triangle.cxx 
    t8_cmesh/t8_cmesh_stash.c 
    t8_cmesh/t8_cmesh_vtk_reader.cxx 
    t8_cmesh/t8_cmesh_save.cxx 
    t8_cmesh/t8_cmesh_netcdf.c 
    t8_cmesh/t8_cmesh_trees.c 
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
    t8_forest/t8_forest_netcdf.cxx 
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
    t8_schemes/t8_default/t8_default.cxx
    t8_schemes/t8_default/t8_default_common/t8_default_common.cxx
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
    t8_schemes/t8_default/t8_default_quad/t8_dquad_bits.c
    t8_schemes/t8_default/t8_default_tet/t8_default_tet.cxx
    t8_schemes/t8_default/t8_default_tet/t8_dtet_bits.c
    t8_schemes/t8_default/t8_default_tet/t8_dtet_connectivity.c
    t8_schemes/t8_default/t8_default_tri/t8_default_tri.cxx
    t8_schemes/t8_default/t8_default_tri/t8_dtri_bits.c
    t8_schemes/t8_default/t8_default_tri/t8_dtri_connectivity.c
    t8_schemes/t8_default/t8_default_vertex/t8_default_vertex.cxx
    t8_schemes/t8_default/t8_default_vertex/t8_dvertex_bits.c
    t8_vtk/t8_vtk_polydata.cxx 
    t8_vtk/t8_vtk_unstructured.cxx 
    t8_vtk/t8_vtk_parallel.cxx 
    t8_vtk/t8_vtk_reader.cxx 
    t8_vtk/t8_vtk_writer.cxx
    t8_vtk/t8_vtk_write_ASCII.cxx
    t8_vtk/t8_vtk_writer_helper.cxx
)

target_compile_definitions( T8 PUBLIC T8_CMAKE_BUILD )
target_compile_definitions( T8 PUBLIC T8_CC="${CMAKE_C_COMPILER}" )
target_compile_definitions( T8 PUBLIC T8_CFLAGS="${T8_CFLAGS}" )
target_compile_definitions( T8 PUBLIC T8_CXX="${CMAKE_CXX_COMPILER}" )
target_compile_definitions( T8 PUBLIC T8_CXXFLAGS="${T8_CXXFLAGS}" )
target_compile_definitions( T8 PUBLIC T8_LDFLAGS="${CMAKE_SHARED_LINKER_FLAGS}" )

target_compile_definitions( T8 PUBLIC T8_PACKAGE_STRING="t8 ${T8CODE_VERSION}" )
target_compile_definitions( T8 PUBLIC T8_VERSION="${T8CODE_VERSION}" )
target_compile_definitions( T8 PUBLIC T8_VERSION_MAJOR=${T8CODE_VERSION_MAJOR} )
target_compile_definitions( T8 PUBLIC T8_VERSION_MINOR=${T8CODE_VERSION_MINOR} )
target_compile_definitions( T8 PUBLIC T8_VERSION_PATCH=${T8CODE_VERSION_PATCH} )
target_compile_definitions( T8 PUBLIC T8_VERSION_POINT=${T8CODE_VERSION_POINT} )

get_target_property( T8_LIBS_LIST T8 LINK_LIBRARIES )
string( REPLACE ";" " " T8_LIBS "${T8_LIBS_LIST}" )
target_compile_definitions( T8 PUBLIC T8_LIBS="${T8_LIBS}" )

install( FILES
    t8.h
    t8_cmesh.h
    t8_cmesh.hxx
    t8_cmesh_netcdf.h
    t8_cmesh_readmshfile.h
    t8_cmesh_tetgen.h
    t8_cmesh_triangle.h
    t8_cmesh_vtk_reader.hxx
    t8_eclass.h
    t8_element.h
    t8_element_c_interface.h
    t8_element.hxx
    t8_element_shape.h
    t8_forest_netcdf.h
    t8_mat.h
    t8_mesh.h
    t8_netcdf.h
    t8_refcount.h
    t8_vec.h
    t8_version.h
    t8_vtk.h
    t8_windows.h DESTINATION ${CMAKE_INSTALL_PREFIX}/include
)

install( DIRECTORY t8_cmesh DESTINATION ${CMAKE_INSTALL_PREFIX}/include FILES_MATCHING PATTERN "*.h" )
install( DIRECTORY t8_data DESTINATION ${CMAKE_INSTALL_PREFIX}/include FILES_MATCHING PATTERN "*.h" )
install( DIRECTORY t8_forest DESTINATION ${CMAKE_INSTALL_PREFIX}/include FILES_MATCHING
    PATTERN "*.h"
    PATTERN "*private.h" EXCLUDE )
install( DIRECTORY t8_geometry DESTINATION ${CMAKE_INSTALL_PREFIX}/include FILES_MATCHING PATTERN "*.h" )
install( DIRECTORY t8_schemes DESTINATION ${CMAKE_INSTALL_PREFIX}/include FILES_MATCHING PATTERN "*.h" )
install( DIRECTORY t8_vtk DESTINATION ${CMAKE_INSTALL_PREFIX}/include FILES_MATCHING PATTERN "*.h" )

install( DIRECTORY t8_cmesh DESTINATION ${CMAKE_INSTALL_PREFIX}/include FILES_MATCHING PATTERN "*.hxx" )
install( DIRECTORY t8_data DESTINATION ${CMAKE_INSTALL_PREFIX}/include FILES_MATCHING PATTERN "*.hxx" )
install( DIRECTORY t8_forest DESTINATION ${CMAKE_INSTALL_PREFIX}/include FILES_MATCHING PATTERN "*.hxx" )
install( DIRECTORY t8_geometry DESTINATION ${CMAKE_INSTALL_PREFIX}/include FILES_MATCHING PATTERN "*.hxx" )
install( DIRECTORY t8_schemes DESTINATION ${CMAKE_INSTALL_PREFIX}/include FILES_MATCHING PATTERN "*.hxx" )
install( DIRECTORY t8_vtk DESTINATION ${CMAKE_INSTALL_PREFIX}/include FILES_MATCHING PATTERN "*.hxx" )

install( TARGETS T8 DESTINATION ${CMAKE_INSTALL_PREFIX}/lib )
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
