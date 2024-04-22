if ( ${T8CODE_BUILD_AS_SHARED_LIBRARY} )
    add_library( T8 SHARED )
    set_target_properties( T8 PROPERTIES POSITION_INDEPENDENT_CODE ON )
else()
    add_library( T8 STATIC )
endif()

set_target_properties( T8 PROPERTIES OUTPUT_NAME t8 )

target_compile_definitions( T8 PUBLIC T8_CMAKE_BUILD )
target_compile_definitions( T8 PUBLIC T8_CC="${CMAKE_C_COMPILER}" )
target_compile_definitions( T8 PUBLIC T8_CFLAGS="${CMAKE_C_FLAGS}" )
target_compile_definitions( T8 PUBLIC T8_CPP="${CMAKE_CXX_COMPILER}" )
target_compile_definitions( T8 PUBLIC T8_CPPFLAGS="${CMAKE_CXX_FLAGS}" )
target_compile_definitions( T8 PUBLIC T8_LDFLAGS="${CMAKE_SHARED_LINKER_FLAGS}" )
target_compile_definitions( T8 PUBLIC T8_LIBS="Not available with CMake builds" )


find_package( Git REQUIRED )
execute_process( COMMAND ${GIT_EXECUTABLE} describe --tags 
                 COMMAND cut -c 2- 
                 WORKING_DIRECTORY ${CMAKE_CURRENT_LIST_DIR}
                 OUTPUT_VARIABLE T8CODE_VERSION 
                 OUTPUT_STRIP_TRAILING_WHITESPACE )
execute_process( COMMAND echo ${T8CODE_VERSION} 
                 COMMAND cut -d. -f1 
                 OUTPUT_VARIABLE T8CODE_VERSION_MAJOR 
                 OUTPUT_STRIP_TRAILING_WHITESPACE )
execute_process( COMMAND echo ${T8CODE_VERSION} 
                 COMMAND cut -d. -f2 
                 OUTPUT_VARIABLE T8CODE_VERSION_MINOR 
                 OUTPUT_STRIP_TRAILING_WHITESPACE )
execute_process( COMMAND echo ${T8CODE_VERSION} 
                 COMMAND cut -d. -f3 
                 OUTPUT_VARIABLE T8CODE_VERSION_POINT 
                 OUTPUT_STRIP_TRAILING_WHITESPACE )

target_compile_definitions( T8 PUBLIC T8_PACKAGE_STRING="t8 ${T8CODE_VERSION}" )
target_compile_definitions( T8 PUBLIC T8_VERSION="${T8CODE_VERSION}" )
target_compile_definitions( T8 PUBLIC T8_VERSION_MAJOR=${T8CODE_VERSION_MAJOR} )
target_compile_definitions( T8 PUBLIC T8_VERSION_MINOR=${T8CODE_VERSION_MINOR} )
target_compile_definitions( T8 PUBLIC T8_VERSION_POINT=${T8CODE_VERSION_POINT} )


target_include_directories( T8 PUBLIC ${CMAKE_CURRENT_LIST_DIR} )
target_link_libraries( T8 PUBLIC P4EST::P4EST SC::SC )

if ( CMAKE_BUILD_TYPE STREQUAL "Debug" )
    target_compile_definitions( T8 PUBLIC T8_ENABLE_DEBUG )
endif()

if ( T8CODE_ENABLE_MPI )
    target_compile_definitions( T8 PUBLIC T8_ENABLE_MPI )
    target_compile_definitions( T8 PUBLIC T8_ENABLE_MPIIO )
    target_link_libraries( T8 PUBLIC MPI::MPI_C )
endif()

if( T8CODE_ENABLE_VTK )
    target_compile_definitions( T8 PUBLIC T8_VTK_VERSION_USED="${VTK_MAJOR_VERSION}.${VTK_MINOR_VERSION}" )
    target_compile_definitions( T8 PUBLIC T8_WITH_VTK=1 )
    target_include_directories( T8 PUBLIC ${VTK_INCLUDE_DIRS} )
    target_link_libraries( T8 PUBLIC ${VTK_LIBRARIES} )
endif()

target_sources( T8 PRIVATE 
    t8.c 
    t8_eclass.c 
    t8_element.c 
    t8_element_shape.c
    t8_element_cxx.cxx 
    t8_element_c_interface.cxx 
    t8_mesh.c
    t8_netcdf.c 
    t8_refcount.c 
    t8_version.c 
    t8_vtk.c 
    t8_cmesh/t8_cmesh.cxx
    t8_cmesh/t8_cmesh_cad.cxx 
    t8_cmesh/t8_cmesh_cxx.cxx 
    t8_cmesh/t8_cmesh_triangle.cxx 
    t8_cmesh/t8_cmesh_vtk_writer.c 
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
    t8_cmesh/t8_cmesh_helpers.c 
    t8_cmesh/t8_cmesh_offset.c 
    t8_cmesh/t8_cmesh_readmshfile.cxx 
    t8_data/t8_shmem.c 
    t8_data/t8_containers.cxx 
    t8_forest/t8_forest.c 
    t8_forest/t8_forest_adapt.cxx 
    t8_forest/t8_forest_partition.cxx 
    t8_forest/t8_forest_cxx.cxx 
    t8_forest/t8_forest_private.c 
    t8_forest/t8_forest_vtk.cxx 
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
    t8_schemes/t8_default/t8_default_cxx.cxx
    t8_schemes/t8_default/t8_default_common/t8_default_common_cxx.cxx
    t8_schemes/t8_default/t8_default_hex/t8_default_hex_cxx.cxx
    t8_schemes/t8_default/t8_default_hex/t8_dhex_bits.c
    t8_schemes/t8_default/t8_default_line/t8_default_line_cxx.cxx
    t8_schemes/t8_default/t8_default_line/t8_dline_bits.c
    t8_schemes/t8_default/t8_default_prism/t8_default_prism_cxx.cxx
    t8_schemes/t8_default/t8_default_prism/t8_dprism_bits.c
    t8_schemes/t8_default/t8_default_pyramid/t8_default_pyramid_cxx.cxx
    t8_schemes/t8_default/t8_default_pyramid/t8_dpyramid_bits.c
    t8_schemes/t8_default/t8_default_pyramid/t8_dpyramid_connectivity.c 
    t8_schemes/t8_default/t8_default_quad/t8_default_quad_cxx.cxx
    t8_schemes/t8_default/t8_default_quad/t8_dquad_bits.c
    t8_schemes/t8_default/t8_default_tet/t8_default_tet_cxx.cxx
    t8_schemes/t8_default/t8_default_tet/t8_dtet_bits.c
    t8_schemes/t8_default/t8_default_tet/t8_dtet_connectivity.c
    t8_schemes/t8_default/t8_default_tri/t8_default_tri_cxx.cxx
    t8_schemes/t8_default/t8_default_tri/t8_dtri_bits.c
    t8_schemes/t8_default/t8_default_tri/t8_dtri_connectivity.c
    t8_schemes/t8_default/t8_default_vertex/t8_default_vertex_cxx.cxx
    t8_schemes/t8_default/t8_default_vertex/t8_dvertex_bits.c
    t8_vtk/t8_vtk_polydata.cxx 
    t8_vtk/t8_vtk_unstructured.cxx 
    t8_vtk/t8_vtk_parallel.cxx 
    t8_vtk/t8_vtk_reader.cxx 
)

set( T8_PUBLIC_HEADERS 
    t8.h 
    t8_eclass.h 
    t8_mesh.h 
    t8_element_cxx.hxx 
    t8_element.h
    t8_element_c_interface.h 
    t8_refcount.h 
    t8_cmesh.h t8_cmesh_triangle.h 
    t8_cmesh_tetgen.h t8_cmesh_readmshfile.h
    t8_cmesh_vtk_writer.h
    t8_cmesh_vtk_reader.hxx 
    t8_vec.h 
    t8_version.h 
    t8_vtk.h 
    t8_cmesh_netcdf.h 
    t8_forest_netcdf.h 
    t8_element_shape.h 
    t8_netcdf.h
    t8_cmesh/t8_cmesh_save.h
    t8_cmesh/t8_cmesh_examples.h 
    t8_cmesh/t8_cmesh_geometry.h 
    t8_cmesh/t8_cmesh_helpers.h 
    t8_cmesh/t8_cmesh_cad.hxx
    t8_data/t8_shmem.h 
    t8_data/t8_containers.h
    t8_forest/t8_forest.h
    t8_forest/t8_forest_general.h
    t8_forest/t8_forest_geometrical.h
    t8_forest/t8_forest_profiling.h
    t8_forest/t8_forest_io.h
    t8_forest/t8_forest_adapt.h 
    t8_forest/t8_forest_vtk.h 
    t8_forest/t8_forest_to_vtkUnstructured.hxx
    t8_forest/t8_forest_iterate.h 
    t8_forest/t8_forest_partition.h
    t8_geometry/t8_geometry.h
    t8_geometry/t8_geometry_base.hxx 
    t8_geometry/t8_geometry_base.h 
    t8_geometry/t8_geometry_with_vertices.hxx 
    t8_geometry/t8_geometry_with_vertices.h 
    t8_geometry/t8_geometry_helpers.h
    t8_geometry/t8_geometry_implementations/t8_geometry_linear.h
    t8_geometry/t8_geometry_implementations/t8_geometry_linear_axis_aligned.h
    t8_geometry/t8_geometry_implementations/t8_geometry_analytic.hxx
    t8_geometry/t8_geometry_implementations/t8_geometry_examples.h
    t8_geometry/t8_geometry_implementations/t8_geometry_cad.h
    t8_geometry/t8_geometry_implementations/t8_geometry_cad.hxx
    t8_geometry/t8_geometry_implementations/t8_geometry_linear.hxx
    t8_geometry/t8_geometry_implementations/t8_geometry_linear_axis_aligned.hxx
    t8_geometry/t8_geometry_implementations/t8_geometry_examples.hxx
    t8_geometry/t8_geometry_implementations/t8_geometry_zero.hxx 
    t8_vtk/t8_vtk_reader.hxx 
    t8_vtk/t8_vtk_types.h
)

install( FILES ${T8_PUBLIC_HEADERS} DESTINATION ${CMAKE_INSTALL_PREFIX}/include )
install( TARGETS T8 DESTINATION ${CMAKE_INSTALL_PREFIX}/lib )
