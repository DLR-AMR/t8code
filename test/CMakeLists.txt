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

# Function to add a test executable and register it with CTest.
# The test executable is built from the provided source files.
function( add_t8_test )
    set( options "" )
    set( oneValueArgs "NAME" )
    set( multiValueArgs "SOURCES" )
    cmake_parse_arguments( ADD_T8_TEST "${options}" "${oneValueArgs}" "${multiValueArgs}" ${ARGN} )

    # Get the path of the first file listed in the SOURCES list and use it to determine the build directory.
    list (LENGTH ADD_T8_TEST_SOURCES TEST_SOURCES_LENGTH)
    list(GET ADD_T8_TEST_SOURCES 0 TEST_SOURCE)
    
    get_filename_component(TEST_SOURCE_DIR "${CMAKE_CURRENT_LIST_DIR}/${TEST_SOURCE}" DIRECTORY)
    file(RELATIVE_PATH TEST_RELATIVE_DIR "${CMAKE_SOURCE_DIR}" "${TEST_SOURCE_DIR}")
    set(TEST_BUILD_DIR "${CMAKE_BINARY_DIR}/${TEST_RELATIVE_DIR}")

    add_executable( ${ADD_T8_TEST_NAME} ${ADD_T8_TEST_SOURCES} )

    # Link base libraries.
    target_include_directories( ${ADD_T8_TEST_NAME} PRIVATE ${CMAKE_SOURCE_DIR} )

    target_link_libraries( ${ADD_T8_TEST_NAME} PRIVATE GTest::gtest_main T8 pthread )

    # Check if test is a Fortran or mesh handle file.
    string( FIND ${ADD_T8_TEST_NAME} "fortran" is_fortran_file )
    string( FIND ${TEST_BUILD_DIR} "mesh_handle" is_mesh_handle_file )

    if ( (${is_fortran_file} GREATER_EQUAL 0) AND T8CODE_ENABLE_MPI ) 
        target_link_libraries( ${ADD_T8_TEST_NAME} PRIVATE MPI::MPI_Fortran )
    endif ()
    if ( (${is_mesh_handle_file} GREATER_EQUAL 0) )
        target_link_libraries( ${ADD_T8_TEST_NAME} PRIVATE T8_MESH_HANDLE )
    endif ()

    set_target_properties(${ADD_T8_TEST_NAME} PROPERTIES
        RUNTIME_OUTPUT_DIRECTORY "${TEST_BUILD_DIR}"
        LIBRARY_OUTPUT_DIRECTORY "${TEST_BUILD_DIR}"
        ARCHIVE_OUTPUT_DIRECTORY "${TEST_BUILD_DIR}"
    )

    if( T8CODE_EXPORT_COMPILE_COMMANDS )
        set_target_properties( ${ADD_T8_TEST_NAME} PROPERTIES EXPORT_COMPILE_COMMANDS ON )
    endif( T8CODE_EXPORT_COMPILE_COMMANDS )

    # Split custom test command into a list by whitespace.
    # If MPI is enabled, and no custom test command is provided, use mpirun -np 4 as the default command.
    # If MPI is not enabled, use the custom test command as is.
    if ( T8CODE_ENABLE_MPI )
        if (${ADD_T8_TEST_NAME} MATCHES "_serial")
            separate_arguments(T8CODE_TEST_COMMAND_LIST NATIVE_COMMAND ${T8CODE_CUSTOM_SERIAL_TEST_COMMAND})
        elseif ( T8CODE_CUSTOM_PARALLEL_TEST_COMMAND STREQUAL "" )
            separate_arguments (T8CODE_TEST_COMMAND_LIST NATIVE_COMMAND "mpirun -np 4")
        else ( T8CODE_CUSTOM_PARALLEL_TEST_COMMAND STREQUAL "" )
            separate_arguments (T8CODE_TEST_COMMAND_LIST NATIVE_COMMAND ${T8CODE_CUSTOM_PARALLEL_TEST_COMMAND})
        endif ()
    else( T8CODE_ENABLE_MPI )
        separate_arguments (T8CODE_TEST_COMMAND_LIST NATIVE_COMMAND ${T8CODE_CUSTOM_SERIAL_TEST_COMMAND})
    endif ( T8CODE_ENABLE_MPI )
    add_test( NAME ${ADD_T8_TEST_NAME} COMMAND ${T8CODE_TEST_COMMAND_LIST} ${TEST_BUILD_DIR}/${ADD_T8_TEST_NAME} )
    if ( T8CODE_CODE_COVERAGE )
        # Extent time that a single test is allowed to run to collect coverage information (as this naturally takes more time).
        # The default timeout of 1500 seconds is actually exceeded for some tests.
        # Without the collection of coverage information, an extension of the default timeout is not needed.
        set_tests_properties(${ADD_T8_TEST_NAME} PROPERTIES TIMEOUT 3000)
    endif ( T8CODE_CODE_COVERAGE )
endfunction()

# Function to add a test executable for C++ tests.
# This function wraps the `add_t8_test` function to include common test sources.
# Common test sources are t8_gtest_main.cpp and t8_gtest_memory_macros.cxx
function ( add_t8_cpp_test)
    set( options "" )
    set( oneValueArgs "NAME" )
    set( multiValueArgs "SOURCES" )
    cmake_parse_arguments( ADD_T8_CPP_TEST "${options}" "${oneValueArgs}" "${multiValueArgs}" ${ARGN} )

    # add cpp sources to SOURCES list and call add_t8_test function
    list ( APPEND ADD_T8_CPP_TEST_SOURCES t8_gtest_main.cxx t8_gtest_memory_macros.cxx )
    add_t8_test( NAME ${ADD_T8_CPP_TEST_NAME} SOURCES ${ADD_T8_CPP_TEST_SOURCES} )
endfunction()

# Copy test files to build folder so that the t8_test programs can find them.
function( copy_test_file TEST_FILE_NAME )
    configure_file(${CMAKE_CURRENT_LIST_DIR}/testfiles/${TEST_FILE_NAME} ${CMAKE_CURRENT_BINARY_DIR}/testfiles/${TEST_FILE_NAME} COPYONLY)
endfunction()

add_t8_cpp_test( NAME t8_gtest_cmesh_bcast_parallel     SOURCES t8_cmesh/t8_gtest_bcast.cxx )
add_t8_cpp_test( NAME t8_gtest_eclass_serial            SOURCES t8_gtest_eclass.cxx )
add_t8_cpp_test( NAME t8_gtest_mat_serial               SOURCES t8_gtest_mat.cxx )
add_t8_cpp_test( NAME t8_gtest_refcount_serial          SOURCES t8_gtest_refcount.cxx )
add_t8_cpp_test( NAME t8_gtest_occ_linkage_serial       SOURCES t8_gtest_occ_linkage.cxx )
add_t8_cpp_test( NAME t8_gtest_version_serial           SOURCES t8_gtest_version.cxx )
add_t8_cpp_test( NAME t8_gtest_basics_serial            SOURCES t8_gtest_basics.cxx )
add_t8_cpp_test( NAME t8_gtest_vtk_linkage_serial       SOURCES t8_gtest_vtk_linkage.cxx )

add_t8_cpp_test( NAME t8_gtest_type_serial                  SOURCES t8_types/t8_gtest_type.cxx )
add_t8_cpp_test( NAME t8_gtest_random_accessible_serial     SOURCES t8_types/t8_gtest_random_accessible.cxx )
add_t8_cpp_test( NAME t8_gtest_vec_serial                   SOURCES t8_types/t8_gtest_vec.cxx )

add_t8_cpp_test( NAME t8_gtest_hypercube_parallel                           SOURCES t8_cmesh/t8_gtest_hypercube.cxx )
add_t8_cpp_test( NAME t8_gtest_cmesh_readmshfile_serial                     SOURCES t8_cmesh/t8_gtest_cmesh_readmshfile.cxx )
add_t8_cpp_test( NAME t8_gtest_cmesh_copy_serial                            SOURCES t8_cmesh/t8_gtest_cmesh_copy.cxx )
add_t8_cpp_test( NAME t8_gtest_cmesh_face_is_boundary_parallel              SOURCES t8_cmesh/t8_gtest_cmesh_face_is_boundary.cxx )
add_t8_cpp_test( NAME t8_gtest_cmesh_partition_parallel                     SOURCES t8_cmesh/t8_gtest_cmesh_partition.cxx )
add_t8_cpp_test( NAME t8_gtest_cmesh_set_partition_offsets_parallel         SOURCES t8_cmesh/t8_gtest_cmesh_set_partition_offsets.cxx )
add_t8_cpp_test( NAME t8_gtest_cmesh_set_join_by_vertices_serial            SOURCES t8_cmesh/t8_gtest_cmesh_set_join_by_vertices.cxx )
add_t8_cpp_test( NAME t8_gtest_cmesh_add_attributes_when_derive_parallel    SOURCES t8_cmesh/t8_gtest_cmesh_add_attributes_when_derive.cxx )
add_t8_cpp_test( NAME t8_gtest_cmesh_vertex_conn_tree_to_vertex_parallel    SOURCES t8_cmesh/t8_gtest_cmesh_vertex_conn_tree_to_vertex.cxx )
add_t8_cpp_test( NAME t8_gtest_cmesh_vertex_conn_vertex_to_tree_parallel    SOURCES t8_cmesh/t8_gtest_cmesh_vertex_conn_vertex_to_tree.cxx )
add_t8_cpp_test( NAME t8_gtest_cmesh_vertex_conn_serial                     SOURCES t8_cmesh/t8_gtest_cmesh_vertex_conn.cxx )
add_t8_cpp_test( NAME t8_gtest_compute_first_element_serial                 SOURCES t8_cmesh/t8_gtest_compute_first_element.cxx )
add_t8_cpp_test( NAME t8_gtest_multiple_attributes_parallel                 SOURCES t8_cmesh/t8_gtest_multiple_attributes.cxx )
add_t8_cpp_test( NAME t8_gtest_attribute_gloidx_array_serial                SOURCES t8_cmesh/t8_gtest_attribute_gloidx_array.cxx )
add_t8_cpp_test( NAME t8_gtest_cmesh_bounding_box_serial                    SOURCES t8_cmesh/t8_gtest_cmesh_bounding_box.cxx )

add_t8_cpp_test( NAME t8_gtest_shmem_parallel  SOURCES t8_data/t8_gtest_shmem.cxx )
add_t8_cpp_test( NAME t8_gtest_data_pack_parallel       SOURCES t8_data/t8_gtest_data_handler.cxx t8_data/t8_data_handler_specs.cxx)

add_t8_cpp_test( NAME t8_gtest_element_volume_serial              SOURCES t8_forest/t8_gtest_element_volume.cxx )
add_t8_cpp_test( NAME t8_gtest_search_parallel                    SOURCES t8_forest/t8_gtest_search.cxx )
add_t8_cpp_test( NAME t8_gtest_half_neighbors_parallel            SOURCES t8_forest/t8_gtest_half_neighbors.cxx )
add_t8_cpp_test( NAME t8_gtest_find_owner_parallel                SOURCES t8_forest/t8_gtest_find_owner.cxx )
add_t8_cpp_test( NAME t8_gtest_user_data_parallel                 SOURCES t8_forest/t8_gtest_user_data.cxx )
add_t8_cpp_test( NAME t8_gtest_transform_serial                   SOURCES t8_forest/t8_gtest_transform.cxx )
add_t8_cpp_test( NAME t8_gtest_ghost_exchange_parallel            SOURCES t8_forest/t8_gtest_ghost_exchange.cxx )
add_t8_cpp_test( NAME t8_gtest_ghost_delete_parallel              SOURCES t8_forest/t8_gtest_ghost_delete.cxx )
add_t8_cpp_test( NAME t8_gtest_ghost_and_owner_parallel           SOURCES t8_forest/t8_gtest_ghost_and_owner.cxx )
add_t8_cpp_test( NAME t8_gtest_balance_parallel                   SOURCES t8_forest/t8_gtest_balance.cxx )
add_t8_cpp_test( NAME t8_gtest_forest_commit_parallel             SOURCES t8_forest/t8_gtest_forest_commit.cxx )
add_t8_cpp_test( NAME t8_gtest_forest_face_normal_serial          SOURCES t8_forest/t8_gtest_forest_face_normal.cxx )
add_t8_cpp_test( NAME t8_gtest_element_is_leaf_serial             SOURCES t8_forest/t8_gtest_element_is_leaf.cxx )
add_t8_cpp_test( NAME t8_gtest_partition_data_parallel            SOURCES t8_forest/t8_gtest_partition_data.cxx )
add_t8_cpp_test( NAME t8_gtest_set_partition_offset_parallel      SOURCES t8_forest/t8_gtest_set_partition_offset.cxx )
add_t8_cpp_test( NAME t8_gtest_partition_for_coarsening_parallel  SOURCES t8_forest/t8_gtest_partition_for_coarsening.cxx )
add_t8_cpp_test( NAME t8_gtest_weighted_partitioning_parallel     SOURCES t8_forest/t8_gtest_weighted_partitioning.cxx )

add_t8_cpp_test( NAME t8_gtest_permute_hole_serial          SOURCES t8_forest_incomplete/t8_gtest_permute_hole.cxx )
add_t8_cpp_test( NAME t8_gtest_recursive_serial             SOURCES t8_forest_incomplete/t8_gtest_recursive.cxx )
add_t8_cpp_test( NAME t8_gtest_iterate_replace_parallel     SOURCES t8_forest_incomplete/t8_gtest_iterate_replace.cxx )
add_t8_cpp_test( NAME t8_gtest_empty_local_tree_parallel    SOURCES t8_forest_incomplete/t8_gtest_empty_local_tree.cxx )
add_t8_cpp_test( NAME t8_gtest_empty_global_tree_parallel   SOURCES t8_forest_incomplete/t8_gtest_empty_global_tree.cxx )

if( CMAKE_BUILD_TYPE STREQUAL "Debug" )
    add_t8_cpp_test( NAME t8_gtest_geometry_negative_volume_serial      SOURCES t8_geometry/t8_gtest_geometry_negative_volume.cxx )
endif()

if( T8_ENABLE_OCC )
    add_t8_cpp_test( NAME t8_gtest_geometry_cad_serial                  SOURCES t8_geometry/t8_geometry_implementations/t8_gtest_geometry_cad.cxx )
endif()
add_t8_cpp_test( NAME t8_gtest_geometry_linear_serial                   SOURCES t8_geometry/t8_geometry_implementations/t8_gtest_geometry_linear.cxx )
add_t8_cpp_test( NAME t8_gtest_geometry_lagrange_serial                 SOURCES t8_geometry/t8_geometry_implementations/t8_gtest_geometry_lagrange.cxx )
add_t8_cpp_test( NAME t8_gtest_geometry_triangular_interpolation_serial SOURCES t8_geometry/t8_gtest_geometry_triangular_interpolation.cxx )
add_t8_cpp_test( NAME t8_gtest_geometry_handling_serial                 SOURCES t8_geometry/t8_gtest_geometry_handling.cxx )
add_t8_cpp_test( NAME t8_gtest_point_inside_serial                      SOURCES t8_geometry/t8_gtest_point_inside.cxx )

add_t8_cpp_test( NAME t8_gtest_vtk_writer_parallel          SOURCES t8_IO/t8_gtest_vtk_writer.cxx )
if( T8CODE_ENABLE_VTK )
    add_t8_cpp_test( NAME t8_gtest_vtk_reader_parallel      SOURCES t8_IO/t8_gtest_vtk_reader.cxx )
endif()

add_t8_cpp_test( NAME t8_gtest_nca_serial                   SOURCES t8_schemes/t8_gtest_nca.cxx )
add_t8_cpp_test( NAME t8_gtest_pyra_connectivity_serial     SOURCES t8_schemes/t8_gtest_pyra_connectivity.cxx )
add_t8_cpp_test( NAME t8_gtest_face_neigh_serial            SOURCES t8_schemes/t8_gtest_face_neigh.cxx )
add_t8_cpp_test( NAME t8_gtest_get_linear_id_serial         SOURCES t8_schemes/t8_gtest_get_linear_id.cxx )
add_t8_cpp_test( NAME t8_gtest_ancestor_serial              SOURCES t8_schemes/t8_gtest_ancestor.cxx )
add_t8_cpp_test( NAME t8_gtest_ancestor_id_serial           SOURCES t8_schemes/t8_gtest_ancestor_id.cxx )
add_t8_cpp_test( NAME t8_gtest_element_count_leaves_serial  SOURCES t8_schemes/t8_gtest_element_count_leaves.cxx )
add_t8_cpp_test( NAME t8_gtest_element_ref_coords_serial    SOURCES t8_schemes/t8_gtest_element_ref_coords.cxx )
add_t8_cpp_test( NAME t8_gtest_descendant_serial            SOURCES t8_schemes/t8_gtest_descendant.cxx )
add_t8_cpp_test( NAME t8_gtest_find_parent_serial           SOURCES t8_schemes/t8_gtest_find_parent.cxx )
add_t8_cpp_test( NAME t8_gtest_equal_serial                 SOURCES t8_schemes/t8_gtest_equal.cxx )
add_t8_cpp_test( NAME t8_gtest_successor_serial             SOURCES t8_schemes/t8_gtest_successor.cxx )
add_t8_cpp_test( NAME t8_gtest_boundary_extrude_serial      SOURCES t8_schemes/t8_gtest_boundary_extrude.cxx )
add_t8_cpp_test( NAME t8_gtest_face_descendant_serial       SOURCES t8_schemes/t8_gtest_face_descendant.cxx )
add_t8_cpp_test( NAME t8_gtest_default_serial               SOURCES t8_schemes/t8_gtest_default.cxx )
add_t8_cpp_test( NAME t8_gtest_child_parent_face_serial     SOURCES t8_schemes/t8_gtest_child_parent_face.cxx )
add_t8_cpp_test( NAME t8_gtest_pack_unpack_serial           SOURCES t8_schemes/t8_gtest_pack_unpack.cxx )
add_t8_cpp_test( NAME t8_gtest_root_serial                  SOURCES t8_schemes/t8_gtest_root.cxx )
add_t8_cpp_test( NAME t8_gtest_scheme_consistency_serial    SOURCES t8_schemes/t8_gtest_scheme_consistency.cxx )
add_t8_cpp_test( NAME t8_gtest_input_equal_output_serial    SOURCES t8_schemes/t8_gtest_input_equal_output.cxx )
add_t8_cpp_test( NAME t8_gtest_face_corner_serial           SOURCES t8_schemes/t8_gtest_face_corner.cxx )
add_t8_cpp_test( NAME t8_gtest_set_linear_id_serial         SOURCES t8_schemes/t8_gtest_set_linear_id.cxx )
add_t8_cpp_test( NAME t8_gtest_elements_are_family_serial   SOURCES t8_schemes/t8_gtest_elements_are_family.cxx )
if( T8CODE_BUILD_FORTRAN_INTERFACE AND T8CODE_ENABLE_MPI )
  add_t8_test( NAME t8_test_fortran_mpi_interface_init_parallel    SOURCES api/t8_fortran_interface/t8_test_mpi_init.f90 )
endif()
 
add_t8_cpp_test( NAME t8_gtest_vector_split_serial          SOURCES t8_helper_functions/t8_gtest_vector_split.cxx )

if( T8CODE_BUILD_MESH_HANDLE )
  add_t8_cpp_test( NAME t8_gtest_mesh_handle_parallel             SOURCES mesh_handle/t8_gtest_mesh_handle.cxx )
  add_t8_cpp_test( NAME t8_gtest_compare_handle_to_forest_serial  SOURCES mesh_handle/t8_gtest_compare_handle_to_forest.cxx )
  add_t8_cpp_test( NAME t8_gtest_handle_ghost_parallel            SOURCES mesh_handle/t8_gtest_ghost.cxx )
  add_t8_cpp_test( NAME t8_gtest_custom_competence_serial         SOURCES mesh_handle/t8_gtest_custom_competence.cxx )
  add_t8_cpp_test( NAME t8_gtest_cache_competence_serial          SOURCES mesh_handle/t8_gtest_cache_competence.cxx )
  add_t8_cpp_test( NAME t8_gtest_handle_data_parallel             SOURCES mesh_handle/t8_gtest_handle_data.cxx )
endif()

copy_test_file( test_cube_unstructured_1.inp )
copy_test_file( test_cube_unstructured_2.inp )
copy_test_file( test_vtk_tri.vtu )
copy_test_file( test_vtk_cube.vtp )
copy_test_file( test_parallel_file.pvtu )
copy_test_file( test_parallel_file_0.vtu )
copy_test_file( test_parallel_file_1.vtu )
copy_test_file( test_polydata.pvtp )
copy_test_file( test_polydata_0.vtp )
copy_test_file( test_polydata_1.vtp )
copy_test_file( test_msh_file_vers4_ascii.msh )
copy_test_file( test_msh_file_vers4_bin.msh )
