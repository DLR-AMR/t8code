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
    if( T8CODE_ENABLE_PROFILING_BARRIERS )
        target_compile_definitions( ${ADD_T8_BENCHMARK_NAME} PUBLIC T8_ENABLE_PROFILE_BARRIER=1 )
    endif( T8CODE_ENABLE_PROFILING_BARRIERS )

    install( TARGETS ${ADD_T8_BENCHMARK_NAME} DESTINATION ${CMAKE_INSTALL_BINDIR} )
endfunction()

add_t8_benchmark( NAME t8_time_partition SOURCES t8_time_partition.cxx )
add_t8_benchmark( NAME t8_time_forest_partition SOURCES t8_time_forest_partition.cxx )
add_t8_benchmark( NAME t8_time_prism_adapt SOURCES t8_time_prism_adapt.cxx )
add_t8_benchmark( NAME t8_time_fractal SOURCES t8_time_fractal.cxx )
add_t8_benchmark( NAME t8_time_set_join_by_vertices SOURCES t8_time_set_join_by_vertices.cxx )
add_t8_benchmark( NAME t8_time_new_refine SOURCES t8_time_new_refine.c )
add_t8_benchmark( NAME t8_time_search SOURCES t8_time_search.cxx )
