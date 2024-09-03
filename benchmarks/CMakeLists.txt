function( add_t8_benchmark )
    set( options "" )
    set( oneValueArgs "NAME" )
    set( multiValueArgs "SOURCES" )
    cmake_parse_arguments( ADD_T8_BENCHMARK "${options}" "${oneValueArgs}" "${multiValueArgs}" ${ARGN} )

    add_executable( ${ADD_T8_BENCHMARK_NAME} ${ADD_T8_BENCHMARK_SOURCES} )
    target_link_libraries( ${ADD_T8_BENCHMARK_NAME} PRIVATE T8 ${SC_LIBRARIES} m )
    target_include_directories( ${ADD_T8_BENCHMARK_NAME} PRIVATE ${PROJECT_SOURCE_DIR} ${SC_INCLUDE_DIR})

    if( T8CODE_EXPORT_COMPILE_COMMANDS )
        set_target_properties( ${ADD_T8_BENCHMARK_NAME} PROPERTIES EXPORT_COMPILE_COMMANDS ON )
    endif( T8CODE_EXPORT_COMPILE_COMMANDS )

    install( TARGETS ${ADD_T8_BENCHMARK_NAME} DESTINATION ${CMAKE_INSTALL_PREFIX}/bin )
endfunction()

add_t8_benchmark( NAME t8_time_partition SOURCES time_partition.cxx )
add_t8_benchmark( NAME t8_time_forest_partition SOURCES time_forest_partition.cxx )
add_t8_benchmark( NAME t8_time_prism_adapt SOURCES t8_time_prism_adapt.cxx )
add_t8_benchmark( NAME t8_time_fractal SOURCES t8_time_fractal.cxx )
add_t8_benchmark( NAME t8_time_set_join_by_vertices SOURCES t8_time_set_join_by_vertices.cxx )
add_t8_benchmark( NAME t8_time_new_refine SOURCES time_new_refine.c )
add_t8_benchmark( NAME t8_bunny SOURCES ExtremeScaling/bunny.cxx )
