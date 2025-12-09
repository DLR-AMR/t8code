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

set(CPACK_PACKAGE_VENDOR "DLR-SC AMR")
set(CPACK_PACKAGE_NAME "T8CODE")
set(CPACK_PACKAGE_DESCRIPTION_SUMMARY "Parallel algorithms and data structures for tree-based AMR with arbitrary element shapes.")
set(CPACK_PACKAGE_VERSION_MAJOR ${T8CODE_VERSION_MAJOR})
set(CPACK_PACKAGE_VERSION_MINOR ${T8CODE_VERSION_MINOR})
set(CPACK_PACKAGE_VERSION_PATCH ${T8CODE_VERSION_PATCH})
set(CPACK_RESOURCE_FILE_README "${CMAKE_CURRENT_SOURCE_DIR}/README.md")
set(CPACK_RESOURCE_FILE_LICENSE "${CMAKE_CURRENT_SOURCE_DIR}/COPYING")
set(CPACK_PACKAGE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/package)
set(CPACK_PACKAGE_ICON  ${CMAKE_CURRENT_SOURCE_DIR}/t8code_logo.png)

# Define a variable for the version file
set(VERSION_FILE "${CMAKE_CURRENT_SOURCE_DIR}/version.txt")


# Custom command to generate the version file
add_custom_command(
    OUTPUT ${VERSION_FILE}
    COMMAND ${CMAKE_COMMAND} -E echo "Version ${T8CODE_VERSION_MAJOR}.${T8CODE_VERSION_MINOR}.${T8CODE_VERSION_PATCH}" > ${VERSION_FILE}
    DEPENDS ${CMAKE_CURRENT_SOURCE_DIR}/CMakeLists.txt  # Change as needed
    )

# Create a custom target to ensure the version file is generated
add_custom_target(GenerateVersionFile ALL DEPENDS ${VERSION_FILE})

set(CPACK_SOURCE_GENERATOR "TGZ;ZIP")
set(CPACK_SOURCE_INCLUDE_FILES ${VERSION_FILE})
set(CPACK_SOURCE_IGNORE_FILES .git/ .github/ .vscode/ _CPack_Packages/
.gitmodules .gitignore
${PROJECT_BINARY_DIR}/
bin/
DartConfiguration.tcl
CMakeCache.txt
build/
compile_commands.json
)

set(CPACK_PACKAGE_NAME "T8CODE")
set(CPACK_VERBATIM_VARIABLES TRUE)


include(CPack)

