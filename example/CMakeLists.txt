# Conditionally adding new sources to the main T8 target as in the original BS looked wrong...
if ( ${T8CODE_BUILD_AS_SHARED_LIBRARY} )
    add_library( t8example SHARED )
else()
    add_library( t8example STATIC )
endif()
target_sources( t8example PRIVATE common/t8_example_common.cxx common/t8_example_common_functions.cxx )
target_include_directories( t8example PRIVATE ${CMAKE_CURRENT_LIST_DIR}/.. ${SC_INCLUDE_DIR} )
target_link_libraries( t8example PRIVATE T8 ${SC_LIBRARIES} m )
install( TARGETS t8example DESTINATION lib )

function( add_t8_example )
    set( options "" )
    set( oneValueArgs "NAME" )  
    set( multiValueArgs "SOURCES" )
    cmake_parse_arguments( ADD_T8_EXAMPLE "${options}" "${oneValueArgs}" "${multiValueArgs}" ${ARGN} )

    # Get the path of the first file listed in the SOURCES list and use it to determine the build directory.
    # The executable will be build in the same directory as the first source file.
    list(GET ADD_T8_EXAMPLE_SOURCES 0 FIRST_SOURCE)
    get_filename_component(EXAMPLE_SOURCE_DIR "${CMAKE_CURRENT_LIST_DIR}/${FIRST_SOURCE}" DIRECTORY)
    file(RELATIVE_PATH EXAMPLE_RELATIVE_DIR "${CMAKE_SOURCE_DIR}" "${EXAMPLE_SOURCE_DIR}")
    set(EXAMPLE_BUILD_DIR "${CMAKE_BINARY_DIR}/${EXAMPLE_RELATIVE_DIR}")

    add_executable( ${ADD_T8_EXAMPLE_NAME} ${ADD_T8_EXAMPLE_SOURCES} )
    target_link_libraries( ${ADD_T8_EXAMPLE_NAME} PRIVATE T8 t8example SC::SC )
    target_include_directories( ${ADD_T8_EXAMPLE_NAME} PRIVATE ${CMAKE_CURRENT_LIST_DIR}/.. )

    set_target_properties(${ADD_T8_EXAMPLE_NAME} PROPERTIES
        RUNTIME_OUTPUT_DIRECTORY "${EXAMPLE_BUILD_DIR}"
        LIBRARY_OUTPUT_DIRECTORY "${EXAMPLE_BUILD_DIR}"
        ARCHIVE_OUTPUT_DIRECTORY "${EXAMPLE_BUILD_DIR}"
    )

    if( T8CODE_EXPORT_COMPILE_COMMANDS )
        set_target_properties( ${ADD_T8_EXAMPLE_NAME} PROPERTIES EXPORT_COMPILE_COMMANDS ON )
    endif( T8CODE_EXPORT_COMPILE_COMMANDS )

    install( TARGETS ${ADD_T8_EXAMPLE_NAME} DESTINATION bin )
endfunction()

add_t8_example( NAME t8_advection                       SOURCES advect/t8_advection.cxx )

add_t8_example( NAME t8_cmesh_partition                 SOURCES cmesh/t8_cmesh_partition.cxx )
add_t8_example( NAME t8_cmesh_set_join_by_vertices      SOURCES cmesh/t8_cmesh_set_join_by_vertices.cxx )
add_t8_example( NAME t8_cmesh_geometry_examples         SOURCES cmesh/t8_cmesh_geometry_examples.cxx )
add_t8_example( NAME t8_cmesh_create_partitioned        SOURCES cmesh/t8_cmesh_create_partitioned.cxx )
add_t8_example( NAME t8_cmesh_hypercube_pad             SOURCES cmesh/t8_cmesh_hypercube_pad.cxx )

add_t8_example( NAME t8_test_ghost                      SOURCES forest/t8_test_ghost.cxx )
add_t8_example( NAME t8_test_face_iterate               SOURCES forest/t8_test_face_iterate.cxx )
add_t8_example( NAME t8_test_ghost_large_level_diff     SOURCES forest/t8_test_ghost_large_level_diff.cxx )

add_t8_example( NAME t8_example_geometries              SOURCES geometry/t8_example_geometries.cxx )

add_t8_example( NAME t8_cmesh_load_save                 SOURCES IO/cmesh/t8_cmesh_load_save.cxx )
add_t8_example( NAME t8_read_msh_file                   SOURCES IO/cmesh/gmsh/t8_read_msh_file.cxx )
add_t8_example( NAME t8_load_and_refine_square_w_hole   SOURCES IO/cmesh/gmsh/t8_load_and_refine_square_w_hole.cxx )
add_t8_example( NAME t8_read_tetgen                     SOURCES IO/cmesh/tetgen/t8_read_tetgen_file.cxx )
add_t8_example( NAME t8_time_tetgen                     SOURCES IO/cmesh/tetgen/t8_time_tetgen_file.cxx )
add_t8_example( NAME t8_forest_tetgen                   SOURCES IO/cmesh/tetgen/t8_forest_from_tetgen.cxx )
add_t8_example( NAME t8_read_triangle                   SOURCES IO/cmesh/triangle/t8_read_triangle_file.cxx )
add_t8_example( NAME t8_cmesh_read_from_vtk             SOURCES IO/cmesh/vtk/t8_cmesh_read_from_vtk.cxx )

if(T8CODE_ENABLE_NETCDF)
    add_t8_example( NAME t8_write_cmesh_netcdf              SOURCES IO/cmesh/netcdf/t8_write_cmesh_netcdf.cxx )
    add_t8_example( NAME t8_write_forest_netcdf             SOURCES IO/forest/netcdf/t8_write_forest_netcdf.cxx )
endif()

add_t8_example( NAME t8_gmsh_to_vtk                     SOURCES IO/forest/gmsh/t8_gmsh_to_vtk.cxx )

add_t8_example( NAME t8_example_spheres                 SOURCES remove/t8_example_spheres.cxx )
add_t8_example( NAME t8_example_gauss_blob              SOURCES remove/t8_example_gauss_blob.cxx )
add_t8_example( NAME t8_example_empty_trees             SOURCES remove/t8_example_empty_trees.cxx )

add_t8_example( NAME t8_version                         SOURCES version/t8_version.cxx )
