function( add_t8_benchmark )
    set( options "" )
    set( oneValueArgs "NAME" )
    set( multiValueArgs "SOURCES" )
    cmake_parse_arguments( ADD_T8_BENCHMARK "${options}" "${oneValueArgs}" "${multiValueArgs}" ${ARGN} )

    # Get the path of the first file listed in the SOURCES list and use it to determine the build directory.
    # The executable will be build in the same directory as the first source file.
    list(GET ADD_T8_BENCHMARK_SOURCES 0 FIRST_SOURCE)
    get_filename_component(BENCHMARK_SOURCE_DIR "${CMAKE_CURRENT_LIST_DIR}/${FIRST_SOURCE}" DIRECTORY)
    file(RELATIVE_PATH BENCHMARK_RELATIVE_DIR "${CMAKE_SOURCE_DIR}" "${BENCHMARK_SOURCE_DIR}")
    set(BENCHMARK_BUILD_DIR "${CMAKE_BINARY_DIR}/${BENCHMARK_RELATIVE_DIR}")

    add_executable( ${ADD_T8_BENCHMARK_NAME} ${ADD_T8_BENCHMARK_SOURCES} )
    target_include_directories( ${ADD_T8_BENCHMARK_NAME} PRIVATE ${PROJECT_SOURCE_DIR} )
    target_link_libraries( ${ADD_T8_BENCHMARK_NAME} PRIVATE T8 SC::SC )

    set_target_properties(${ADD_T8_BENCHMARK_NAME} PROPERTIES
        RUNTIME_OUTPUT_DIRECTORY "${BENCHMARK_BUILD_DIR}"
        LIBRARY_OUTPUT_DIRECTORY "${BENCHMARK_BUILD_DIR}"
        ARCHIVE_OUTPUT_DIRECTORY "${BENCHMARK_BUILD_DIR}"
    )

    if( T8CODE_EXPORT_COMPILE_COMMANDS )
        set_target_properties( ${ADD_T8_BENCHMARK_NAME} PROPERTIES EXPORT_COMPILE_COMMANDS ON )
    endif( T8CODE_EXPORT_COMPILE_COMMANDS )

    install( TARGETS ${ADD_T8_BENCHMARK_NAME} DESTINATION bin )
endfunction()

add_t8_benchmark( NAME t8_time_partition SOURCES t8_time_partition.cxx )
add_t8_benchmark( NAME t8_time_forest_partition SOURCES t8_time_forest_partition.cxx )
add_t8_benchmark( NAME t8_time_prism_adapt SOURCES t8_time_prism_adapt.cxx )
add_t8_benchmark( NAME t8_time_fractal SOURCES t8_time_fractal.cxx )
add_t8_benchmark( NAME t8_time_set_join_by_vertices SOURCES t8_time_set_join_by_vertices.cxx )
add_t8_benchmark( NAME t8_time_new_refine SOURCES t8_time_new_refine.c )
