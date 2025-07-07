#  This file is part of t8code.
#  t8code is a C library to manage a collection (a forest) of multiple
#  connected adaptive space-trees of general element types in parallel.
#
#  Copyright (C) 2025 the developers
#
#  t8code is free software; you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation; either version 2 of the License, or
#  (at your option) any later version.
#
#  t8code is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  You should have received a copy of the GNU General Public License
#  along with t8code; if not, write to the Free Software Foundation, Inc.,
#  51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.

if ( ${T8CODE_BUILD_AS_SHARED_LIBRARY} )
    add_library( T8 SHARED )
    set_target_properties( T8 PROPERTIES VERSION ${CMAKE_PROJECT_VERSION} SOVERSION ${T8CODE_VERSION_RAW})
    set_target_properties( T8 PROPERTIES POSITION_INDEPENDENT_CODE ON )
else()
    add_library( T8 STATIC )
endif()

target_compile_features( T8 PUBLIC c_std_11 cxx_std_20 )
add_library( T8CODE::T8 ALIAS T8 )
target_compile_definitions( T8 PUBLIC T8CODE_TEST_LEVEL=${T8CODE_TEST_LEVEL})

if( CMAKE_BUILD_TYPE STREQUAL "Release" )
  set (T8_CXXFLAGS "${CMAKE_CXX_FLAGS} ${CMAKE_CXX_FLAGS_RELEASE}")
  set (T8_CFLAGS "${CMAKE_C_FLAGS} ${CMAKE_C_FLAGS_RELEASE}")
endif()

if( CMAKE_BUILD_TYPE STREQUAL "RelWithDebInfo" )
  set (T8_CXXFLAGS "${CMAKE_CXX_FLAGS} ${CMAKE_CXX_FLAGS_RELWITHDEBINFO}")
  set (T8_CFLAGS "${CMAKE_C_FLAGS} ${CMAKE_C_FLAGS_RELWITHDEBINFO}")
endif()

if( CMAKE_BUILD_TYPE STREQUAL "Debug" )
  # Option to generate code coverage target.
  option( T8CODE_CODE_COVERAGE "Enable code coverage reporting" OFF)
  if(T8CODE_CODE_COVERAGE)
    include(CodeCoverage)
  endif()

  set (T8_CXXFLAGS "${CMAKE_CXX_FLAGS} ${CMAKE_CXX_FLAGS_DEBUG}")
  set (T8_CFLAGS "${CMAKE_C_FLAGS} ${CMAKE_C_FLAGS_DEBUG}")

  target_compile_definitions( T8 PUBLIC T8_ENABLE_DEBUG=1 )
endif()

if( T8CODE_EXPORT_COMPILE_COMMANDS )
  set_target_properties( T8 PROPERTIES EXPORT_COMPILE_COMMANDS ON )
endif( T8CODE_EXPORT_COMPILE_COMMANDS )

if( T8CODE_ENABLE_NETCDF )
  target_link_libraries( T8 PUBLIC netCDF::netcdf )
  target_compile_definitions(T8 PUBLIC
    T8_ENABLE_NETCDF=1
    $<$<AND:$<BOOL:${NETCDF_HAVE_NETCDF_PAR}>,$<BOOL:${T8CODE_ENABLE_MPI}>>:T8_ENABLE_NETCDF_PAR=1> )
  target_sources( T8 PRIVATE
  t8_netcdf.c
  t8_cmesh/t8_cmesh_netcdf.c
  t8_forest/t8_forest_netcdf.cxx
  )
  install( FILES
  t8_cmesh_netcdf.h
  t8_forest_netcdf.h
  t8_netcdf.h
  TYPE INCLUDE )
endif()

set_target_properties( T8 PROPERTIES OUTPUT_NAME t8 )

target_include_directories( T8 PUBLIC
  $<BUILD_INTERFACE:${CMAKE_CURRENT_LIST_DIR}>
  $<INSTALL_INTERFACE:${CMAKE_INSTALL_PREFIX}/include>
)

target_link_libraries( T8 PUBLIC P4EST::P4EST SC::SC )

if ( T8CODE_ENABLE_MPI )
    target_compile_definitions( T8 PUBLIC T8_ENABLE_MPI=1 )
    target_compile_definitions( T8 PUBLIC T8_ENABLE_MPIIO=1 )
    target_link_libraries( T8 PUBLIC MPI::MPI_C )
endif()

if( T8CODE_ENABLE_VTK )
    target_compile_definitions( T8 PUBLIC T8_VTK_VERSION_USED="${VTK_MAJOR_VERSION}.${VTK_MINOR_VERSION}" )
    target_compile_definitions( T8 PUBLIC T8_ENABLE_VTK=1 )
    target_include_directories( T8 PUBLIC ${VTK_INCLUDE_DIR} )
    target_link_libraries( T8 PUBLIC ${VTK_LIBRARIES} )
    target_sources(T8 PRIVATE
      t8_vtk/t8_with_vtk/t8_vtk_parallel.cxx
      t8_vtk/t8_with_vtk/t8_vtk_polydata.cxx
      t8_vtk/t8_with_vtk/t8_vtk_unstructured.cxx
      t8_vtk/t8_with_vtk/t8_vtk_reader.cxx
    )
    install( FILES
      t8_vtk/t8_with_vtk/t8_vtk_parallel.hxx
      t8_vtk/t8_with_vtk/t8_vtk_polydata.hxx
      t8_vtk/t8_with_vtk/t8_vtk_unstructured.hxx
      t8_vtk/t8_with_vtk/t8_vtk_reader.hxx
    DESTINATION include)
endif()

if( T8CODE_ENABLE_OCC )
    target_compile_definitions( T8 PUBLIC T8_ENABLE_OCC=1 )
    target_include_directories( T8 SYSTEM PUBLIC ${OpenCASCADE_INCLUDE_DIR} )
    target_link_libraries( T8 PUBLIC ${OpenCASCADE_LIBRARIES} )
    target_sources(T8 PRIVATE 
      t8_geometry/t8_geometry_implementations/t8_geometry_cad.cxx
    )
    install( FILES 
      t8_geometry/t8_geometry_implementations/t8_geometry_cad.hxx
      t8_geometry/t8_geometry_implementations/t8_geometry_cad.h
      DESTINATION include
    )
endif()

if( T8CODE_ENABLE_MRA )
  target_compile_definitions( T8 PUBLIC T8_ENABLE_MRA )

  # Add GSL library 
  target_link_libraries( T8 PUBLIC GSL::gsl GSL::gslcblas )

  target_sources( T8 PRIVATE 
    t8_mra/t8_gsl.cxx
  )
  install( FILES 
    t8_mra/t8_gsl.hxx
  DESTINATION include )

  # Add unordered_dense library https://github.com/martinus/unordered_dense
  set( UNORDERED_DENSE_DIR "${CMAKE_SOURCE_DIR}/thirdparty/unordered_dense/include" )

  target_include_directories(T8 PUBLIC
    $<BUILD_INTERFACE:${CMAKE_SOURCE_DIR}/src>
    $<BUILD_INTERFACE:${UNORDERED_DENSE_DIR}>
    $<INSTALL_INTERFACE:include>
  )

  target_sources( T8 PRIVATE
    t8_mra/t8_unordered_dense.cxx
    t8_mra/t8_mra.cxx
    t8_mra/vecmat.cxx
    t8_mra/basis_functions.cxx
    t8_mra/num/dunavant.cxx
    t8_mra/mask_coefficients.cxx
  )

  install( FILES
    t8_mra/data/levelindex_map.hpp
    t8_mra/data/cell_data.hpp
    t8_mra/data/levelmultiindex.hpp
    t8_mra/data/triangle_order.hpp
    t8_mra/t8_unordered_dense.hxx
    t8_mra/t8_mra.hxx
    t8_mra/vecmat.hxx
    t8_mra/basis_functions.hxx
    t8_mra/num/dunavant.hxx
    t8_mra/mask_coefficients.hxx
    # t8_mra/levelindex_map.hpp
    # t8_mra/cell_data.hpp
    DESTINATION include)
endif()

if( T8CODE_BUILD_PEDANTIC )
  target_compile_options( T8 PUBLIC -pedantic )
  set (T8_CXXFLAGS "${T8_CXXFLAGS} -Wpedantic")
  set (T8_CFLAGS "${T8_CFLAGS} -Wpedantic")
endif()

if( T8CODE_BUILD_WALL)
  target_compile_options( T8 PUBLIC -Wall )
  set (T8_CXXFLAGS "${T8_CXXFLAGS} -Wall")
  set (T8_CFLAGS "${T8_CFLAGS} -Wall")
endif()

if( T8CODE_BUILD_WERROR )
  target_compile_options( T8 PUBLIC -Werror )
  set (T8_CXXFLAGS "${T8_CXXFLAGS} -Werror")
  set (T8_CFLAGS "${T8_CFLAGS} -Werror")
endif()

if( T8CODE_BUILD_WEXTRA )
  target_compile_options( T8 PUBLIC -Wextra )
  set (T8_CXXFLAGS "${T8_CXXFLAGS} -Wextra")
  set (T8_CFLAGS "${T8_CFLAGS} -Wextra")
endif()


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
    t8_cmesh/t8_cmesh_save.cxx 
    t8_cmesh/t8_cmesh_trees.cxx 
    t8_cmesh/t8_cmesh_commit.cxx 
    t8_cmesh/t8_cmesh_partition.cxx 
    t8_cmesh/t8_cmesh_copy.c 
    t8_cmesh/t8_cmesh_geometry.cxx 
    t8_cmesh/t8_cmesh_examples.cxx 
    t8_cmesh/t8_cmesh_helpers.cxx 
    t8_cmesh/t8_cmesh_offset.c 
    t8_cmesh/t8_cmesh_vertex_connectivity/t8_cmesh_vertex_conn_tree_to_vertex.cxx 
    t8_cmesh/t8_cmesh_vertex_connectivity/t8_cmesh_vertex_conn_vertex_to_tree.cxx 
    t8_cmesh/t8_cmesh_vertex_connectivity/t8_cmesh_vertex_connectivity.cxx 
    t8_cmesh/t8_cmesh_readmshfile.cxx 
    t8_data/t8_shmem.c 
    t8_data/t8_containers.cxx 
    t8_forest/t8_forest_adapt.cxx 
    t8_forest/t8_forest_partition.cxx 
    t8_forest/t8_forest.cxx 
    t8_forest/t8_forest_private.cxx
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
    t8_cmesh_readmshfile.h
    t8_cmesh_tetgen.h
    t8_cmesh_triangle.h
    t8_eclass.h
    t8_element.h
    t8_element_shape.h
    t8_mat.h
    t8_refcount.h
    t8_types/t8_vec.hxx
    t8_version.h
    t8_vtk.h
    TYPE INCLUDE
)

install ( DIRECTORY t8_vector_helper DESTINATION PRIVATE_HEADER FILES_MATCHING PATTERN "*.hxx" )

install( DIRECTORY t8_cmesh TYPE INCLUDE FILES_MATCHING PATTERN "*.h" 
  PATTERN "t8_cmesh/t8_cmesh_copy.h"      EXCLUDE
  PATTERN "t8_cmesh/t8_cmesh_offset.h"    EXCLUDE 
  PATTERN "t8_cmesh/t8_cmesh_partition.h" EXCLUDE
  PATTERN "t8_cmesh/t8_cmesh_trees.h"     EXCLUDE 
  PATTERN "t8_cmesh/t8_cmesh_vertex_connectivity/*.h"     EXCLUDE 
)

install( DIRECTORY t8_data TYPE INCLUDE FILES_MATCHING PATTERN "*.h" )
install( DIRECTORY t8_forest TYPE INCLUDE FILES_MATCHING PATTERN "*.h"
  PATTERN "t8_forest/t8_forest_balance.h" EXCLUDE
  PATTERN "t8_forest/t8_forest_ghost.h "  EXCLUDE
  PATTERN "t8_forest/t8_forest_private.h" EXCLUDE
  PATTERN "t8_forest/t8_forest_types.h"   EXCLUDE
  )
install( DIRECTORY t8_geometry DESTINATION include FILES_MATCHING PATTERN "*.h" 
  PATTERN "t8_geometry/t8_geometry_implementations/t8_geometry_cad.h" EXCLUDE)
install( DIRECTORY t8_schemes DESTINATION include FILES_MATCHING PATTERN "*.h" )
install( DIRECTORY t8_vtk DESTINATION include FILES_MATCHING PATTERN "*.h" )

install( DIRECTORY t8_cmesh DESTINATION include FILES_MATCHING PATTERN "*.hxx"
  PATTERN "t8_cmesh/t8_cmesh_vertex_connectivity/t8_cmesh_vertex_connectivity.hxx" EXCLUDE 
  PATTERN "t8_cmesh/t8_cmesh_vertex_connectivity/t8_cmesh_vertex_conn_vertex_to_tree.hxx" EXCLUDE 
  PATTERN "t8_cmesh/t8_cmesh_vertex_connectivity/t8_cmesh_vertex_conn_tree_to_vertex.hxx" EXCLUDE
)

install( DIRECTORY t8_data DESTINATION include FILES_MATCHING PATTERN "*.hxx" )
install( DIRECTORY t8_forest DESTINATION include FILES_MATCHING PATTERN "*.hxx" )
install( DIRECTORY t8_geometry DESTINATION include FILES_MATCHING PATTERN "*.hxx" 
  PATTERN "t8_geometry/t8_geometry_implementations/t8_geometry_cad.hxx" EXCLUDE)
install( DIRECTORY t8_schemes DESTINATION include FILES_MATCHING PATTERN "*.hxx" )
install( DIRECTORY t8_vtk DESTINATION include FILES_MATCHING PATTERN "*.hxx" 
  PATTERN "t8_vtk/t8_vtk_parallel.hxx"      EXCLUDE
  PATTERN "t8_vtk/t8_vtk_polydata.hxx"      EXCLUDE
  PATTERN "t8_vtk/t8_vtk_unstructured.hxx"  EXCLUDE
  PATTERN "t8_vtk/t8_vtk_write_ASCII.hxx"   EXCLUDE
  PATTERN "t8_vtk/t8_vtk_writer_helper.hxx" EXCLUDE
  )
install( DIRECTORY t8_types TYPE INCLUDE FILES_MATCHING PATTERN "*.hxx" )

install( TARGETS T8 LIBRARY )
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
