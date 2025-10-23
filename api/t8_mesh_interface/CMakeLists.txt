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

# Define own library for mesh interface.
if ( ${T8CODE_BUILD_AS_SHARED_LIBRARY} )
    add_library( T8_MESH_INTERFACE SHARED )
    set_target_properties( T8_MESH_INTERFACE PROPERTIES POSITION_INDEPENDENT_CODE ON )
else()
    add_library( T8_MESH_INTERFACE STATIC )
endif()

add_library( T8CODE::T8_MESH_INTERFACE ALIAS T8_MESH_INTERFACE )
set_target_properties(T8_MESH_INTERFACE PROPERTIES OUTPUT_NAME t8_mesh_interface)

target_include_directories( T8_MESH_INTERFACE PUBLIC
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/..>
  $<INSTALL_INTERFACE:include/t8_mesh_interface>
) 

target_sources(T8_MESH_INTERFACE PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/t8_interface_mesh.cxx)
install( DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR} TYPE INCLUDE FILES_MATCHING PATTERN "*.hxx" )

target_link_libraries(T8_MESH_INTERFACE PUBLIC T8)

install( TARGETS T8_MESH_INTERFACE LIBRARY )
install( TARGETS T8_MESH_INTERFACE EXPORT ${PROJECT_NAME}-targets )

if( T8CODE_EXPORT_COMPILE_COMMANDS )
    set_target_properties( T8_MESH_INTERFACE PROPERTIES EXPORT_COMPILE_COMMANDS ON )
endif( T8CODE_EXPORT_COMPILE_COMMANDS )