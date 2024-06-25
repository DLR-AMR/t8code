# Conditionally adding new sources to the main T8 target as in the original BS looked wrong...
add_library( t8example )
target_sources( t8example PRIVATE common/t8_example_common.cxx common/t8_example_common_functions.cxx )
target_include_directories( t8example PRIVATE ${CMAKE_CURRENT_LIST_DIR}/.. )
target_link_libraries( t8example PRIVATE T8 )
install( TARGETS t8example DESTINATION ${CMAKE_INSTALL_PREFIX}/lib )

function( add_t8_example )
    set( options "" )
    set( oneValueArgs "NAME" )  
    set( multiValueArgs "SOURCES" )
    cmake_parse_arguments( ADD_T8_EXAMPLE "${options}" "${oneValueArgs}" "${multiValueArgs}" ${ARGN} )

    add_executable( ${ADD_T8_EXAMPLE_NAME} ${ADD_T8_EXAMPLE_SOURCES} )
    target_link_libraries( ${ADD_T8_EXAMPLE_NAME} PRIVATE T8 t8example SC::SC )
    target_include_directories( ${ADD_T8_EXAMPLE_NAME} PRIVATE ${CMAKE_CURRENT_LIST_DIR}/.. )

    install( TARGETS ${ADD_T8_EXAMPLE_NAME} DESTINATION ${CMAKE_INSTALL_PREFIX}/example )
endfunction()

add_t8_example( NAME t8_advection                       SOURCES advect/t8_advection.cxx )

add_t8_example( NAME t8_cmesh_partition                 SOURCES cmesh/t8_cmesh_partition.cxx )
add_t8_example( NAME t8_cmesh_set_join_by_vertices      SOURCES cmesh/t8_cmesh_set_join_by_vertices.cxx )
add_t8_example( NAME t8_cmesh_geometry_examples         SOURCES cmesh/t8_cmesh_geometry_examples.cxx )
add_t8_example( NAME t8_cmesh_create_partitioned        SOURCES cmesh/t8_cmesh_create_partitioned.cxx )
add_t8_example( NAME t8_cmesh_hypercube_pad             SOURCES cmesh/t8_cmesh_hypercube_pad.cxx )

add_t8_example( NAME t8_face_neighbor                   SOURCES forest/t8_face_neighbor.cxx )
add_t8_example( NAME t8_test_ghost                      SOURCES forest/t8_test_ghost.cxx )
add_t8_example( NAME t8_test_face_iterate               SOURCES forest/t8_test_face_iterate.cxx )
add_t8_example( NAME t8_test_ghost_large_level_diff     SOURCES forest/t8_test_ghost_large_level_diff.cxx )

add_t8_example( NAME t8_example_geometries              SOURCES geometry/t8_example_geometries.cxx )

add_t8_example( NAME t8_cmesh_load_save                 SOURCES IO/cmesh/t8_cmesh_load_save.cxx )
add_t8_example( NAME t8_read_msh_file                   SOURCES IO/cmesh/gmsh/t8_read_msh_file.cxx )
add_t8_example( NAME t8_load_and_refine_square_w_hole   SOURCES IO/cmesh/gmsh/t8_load_and_refine_square_w_hole.cxx )
add_t8_example( NAME t8_write_cmesh_netcdf              SOURCES IO/cmesh/netcdf/t8_write_cmesh_netcdf.cxx )
add_t8_example( NAME t8_read_tetgen                     SOURCES IO/cmesh/tetgen/t8_read_tetgen_file.cxx )
add_t8_example( NAME t8_time_tetgen                     SOURCES IO/cmesh/tetgen/t8_time_tetgen_file.cxx )
add_t8_example( NAME t8_forest_tetgen                   SOURCES IO/cmesh/tetgen/t8_forest_from_tetgen.cxx )
add_t8_example( NAME t8_read_triangle                   SOURCES IO/cmesh/triangle/t8_read_triangle_file.cxx )
add_t8_example( NAME t8_cmesh_read_from_vtk             SOURCES IO/cmesh/vtk/t8_cmesh_read_from_vtk.cxx )
add_t8_example( NAME t8_netcdf_compilation_status       SOURCES IO/forest/netcdf/t8_netcdf_status.c )
add_t8_example( NAME t8_write_forest_netcdf             SOURCES IO/forest/netcdf/t8_write_forest_netcdf.cxx )

add_t8_example( NAME t8_example_spheres                 SOURCES remove/t8_example_spheres.cxx )
add_t8_example( NAME t8_example_gauss_blob              SOURCES remove/t8_example_gauss_blob.cxx )
add_t8_example( NAME t8_example_empty_trees             SOURCES remove/t8_example_empty_trees.cxx )

add_t8_example( NAME t8_version                         SOURCES version/t8_version.cxx )