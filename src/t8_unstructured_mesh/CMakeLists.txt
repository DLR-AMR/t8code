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
    # Library versioning consistent to t8code.
    set_target_properties(T8_UNSTRUCTURED_MESH PROPERTIES
        VERSION ${T8CODE_VERSION_MAJOR}.${T8CODE_VERSION_MINOR}.${T8CODE_VERSION_PATCH}
        SOVERSION ${T8CODE_VERSION_RAW}
        POSITION_INDEPENDENT_CODE ON
    )
else()
    add_library( T8_MESH_INTERFACE STATIC )
endif()

target_compile_features( T8_MESH_INTERFACE PUBLIC c_std_11 cxx_std_20 )
add_library( T8CODE::T8_MESH_INTERFACE ALIAS T8_MESH_INTERFACE )

install( DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR} TYPE INCLUDE FILES_MATCHING PATTERN "*.hxx" )

target_include_directories( T8 PUBLIC
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/..>
  $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>
) 

target_link_libraries(T8_UNSTRUCTURED_MESH PUBLIC T8)

target_compile_options(T8_UNSTRUCTURED_MESH PRIVATE ${T8_CXXFLAGS})

install( TARGETS T8_MESH_INTERFACE LIBRARY )
install( TARGETS T8_MESH_INTERFACE EXPORT ${PROJECT_NAME}-mesh-interface )