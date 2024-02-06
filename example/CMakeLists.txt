# Conditionally adding new sources to the main T8 target as in the original BS looked wrong...
add_library( T8ex )
target_sources( T8ex PRIVATE common/t8_example_common.cxx common/t8_example_common_functions.cxx )
target_include_directories( T8ex PRIVATE ${CMAKE_CURRENT_LIST_DIR}/.. )
target_link_libraries( T8ex PRIVATE T8 )

function( add_t8_example TEST_NAME TEST_SOURCES )
    add_executable( ${TEST_NAME} ${TEST_SOURCES} )
    target_link_libraries( ${TEST_NAME} PRIVATE T8 T8ex SC::SC )
    target_include_directories( ${TEST_NAME} PRIVATE ${CMAKE_CURRENT_LIST_DIR}/.. )
    if( T8CODE_ENABLE_MPI )
        target_compile_definitions( ${TEST_NAME} PRIVATE T8_ENABLE_MPI )
        target_compile_definitions( ${TEST_NAME} PRIVATE T8_ENABLE_MPIIO )
    endif()
endfunction()

add_t8_example( t8_advection advect/t8_advection.cxx )

add_t8_example( t8_cmesh_partition                      cmesh/t8_cmesh_partition.cxx )
add_t8_example( t8_cmesh_set_join_by_vertices           cmesh/t8_cmesh_set_join_by_vertices.cxx )
add_t8_example( t8_cmesh_geometry_examples              cmesh/t8_cmesh_geometry_examples.cxx )
add_t8_example( t8_cmesh_create_partitioned             cmesh/t8_cmesh_create_partitioned.cxx )
add_t8_example( t8_cmesh_hypercube_pad                  cmesh/t8_cmesh_hypercube_pad.cxx )

add_t8_example( t8_face_neighbor                        forest/t8_face_neighbor.cxx )
add_t8_example( t8_test_ghost                           forest/t8_test_ghost.cxx )
add_t8_example( t8_test_face_iterate                    forest/t8_test_face_iterate.cxx )
add_t8_example( t8_test_ghost_large_level_diff          forest/t8_test_ghost_large_level_diff.cxx )

add_t8_example( t8_example_geometries                   geometry/t8_example_geometries.cxx )