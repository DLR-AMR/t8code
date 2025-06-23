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

function( add_t8_tutorial )
    set( options "" )
    set( oneValueArgs "NAME" )  
    set( multiValueArgs "SOURCES" )
    cmake_parse_arguments( ADD_T8_TUTORIAL "${options}" "${oneValueArgs}" "${multiValueArgs}" ${ARGN} )

    # Get the path of the first file listed in the SOURCES list and use it to determine the build directory.
    # The executable will be build in the same directory as the first source file.
    list(GET ADD_T8_TUTORIAL_SOURCES 0 FIRST_SOURCE)
    get_filename_component(TUTORIAL_SOURCE_DIR "${CMAKE_CURRENT_LIST_DIR}/${FIRST_SOURCE}" DIRECTORY)
    file(RELATIVE_PATH TUTORIAL_RELATIVE_DIR "${CMAKE_SOURCE_DIR}" "${TUTORIAL_SOURCE_DIR}")
    set(TUTORIAL_BUILD_DIR "${CMAKE_BINARY_DIR}/${TUTORIAL_RELATIVE_DIR}")

    add_executable( ${ADD_T8_TUTORIAL_NAME} ${ADD_T8_TUTORIAL_SOURCES} )
    target_link_libraries( ${ADD_T8_TUTORIAL_NAME} PRIVATE T8 SC::SC )
    target_include_directories( ${ADD_T8_TUTORIAL_NAME} PRIVATE ${CMAKE_CURRENT_LIST_DIR}/.. )
    
    set_target_properties(${ADD_T8_TUTORIAL_NAME} PROPERTIES
        RUNTIME_OUTPUT_DIRECTORY "${TUTORIAL_BUILD_DIR}"
        LIBRARY_OUTPUT_DIRECTORY "${TUTORIAL_BUILD_DIR}"
        ARCHIVE_OUTPUT_DIRECTORY "${TUTORIAL_BUILD_DIR}"
    )

    if( T8CODE_EXPORT_COMPILE_COMMANDS )
        set_target_properties( ${ADD_T8_TUTORIAL_NAME} PROPERTIES EXPORT_COMPILE_COMMANDS ON )
    endif( T8CODE_EXPORT_COMPILE_COMMANDS )

    install( TARGETS ${ADD_T8_TUTORIAL_NAME} RUNTIME )
endfunction()

#copy tutorial files to the exact same location in the builddir as they are in the source dir
function( copy_tutorial_file TUTORIAL_FILE_NAME_SUBPATH )
    configure_file(${CMAKE_CURRENT_LIST_DIR}/${TUTORIAL_FILE_NAME_SUBPATH} ${CMAKE_CURRENT_BINARY_DIR}/${TUTORIAL_FILE_NAME_SUBPATH} COPYONLY)
endfunction()

add_t8_tutorial( NAME t8_step0_helloworld               SOURCES general/t8_step0_helloworld.cxx )
add_t8_tutorial( NAME t8_step1_coarsemesh               SOURCES general/t8_step1_coarsemesh.cxx )
add_t8_tutorial( NAME t8_step2_uniform_forest           SOURCES general/t8_step2_uniform_forest.cxx )
add_t8_tutorial( NAME t8_step3_adapt_forest             SOURCES general/t8_step3_main.cxx general/t8_step3_adapt_forest.cxx )
add_t8_tutorial( NAME t8_step4_partition_balance_ghost  SOURCES general/t8_step4_main.cxx general/t8_step3_adapt_forest.cxx general/t8_step4_partition_balance_ghost.cxx )
add_t8_tutorial( NAME t8_step5_element_data             SOURCES general/t8_step5_main.cxx general/t8_step3_adapt_forest.cxx general/t8_step5_element_data.cxx )
add_t8_tutorial( NAME t8_step6_stencil                  SOURCES general/t8_step6_main.cxx general/t8_step3_adapt_forest.cxx general/t8_step6_stencil.cxx )
add_t8_tutorial( NAME t8_step7_interpolation            SOURCES general/t8_step7_main.cxx general/t8_step7_interpolation.cxx )
add_t8_tutorial( NAME t8_tutorial_build_cmesh           SOURCES general/t8_tutorial_build_cmesh.cxx general/t8_tutorial_build_cmesh_main.cxx)
add_t8_tutorial( NAME t8_tutorial_search                SOURCES general/t8_tutorial_search.cxx general/t8_step3_adapt_forest.cxx )
add_t8_tutorial( NAME t8_features_curved_meshes         SOURCES features/t8_features_curved_meshes.cxx)
add_t8_tutorial( NAME t8_tutorial_build_scheme          SOURCES general/t8_tutorial_build_scheme.cxx)   
copy_tutorial_file (features/t8_features_curved_meshes_generate_cmesh_hex.geo)
copy_tutorial_file (features/t8_features_curved_meshes_generate_cmesh_quad.geo)
copy_tutorial_file (features/t8_features_curved_meshes_generate_cmesh_tet.geo)
copy_tutorial_file (features/t8_features_curved_meshes_generate_cmesh_tri.geo)
