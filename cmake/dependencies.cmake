#  This file is part of t8code.
#  t8code is a C library to manage a collection (a forest) of multiple
#  connected adaptive space-trees of general element types in parallel.
#
#  Copyright (C) 2026 the developers
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

include(FetchContent)
set(FETCHCONTENT_QUIET FALSE)


# 1. Read and parse the JSON
file(READ "cmake/dependencies.json" DEPS_JSON)
# 2. Get the number of elements in the "dependencies" array
string(JSON DEPS_COUNT LENGTH "${DEPS_JSON}" "dependencies")
# 3. Loop through the array (subtract 1 because it's 0-indexed)
math(EXPR DEPS_RANGE "${DEPS_COUNT} - 1")

foreach(INDEX RANGE ${DEPS_RANGE})
    # Extract fields
    string(JSON DEP_NAME GET "${DEPS_JSON}" "dependencies" ${INDEX} "name")
    string(JSON DEP_CMAKE_OPTION GET "${DEPS_JSON}" "dependencies" ${INDEX} "depends_on_cmake_option")
    string(JSON DEP_TYPE GET "${DEPS_JSON}" "dependencies" ${INDEX} "source" "type")
    string(JSON DEP_URL  GET "${DEPS_JSON}" "dependencies" ${INDEX} "source" "url")
    string(JSON DEP_REF  GET "${DEPS_JSON}" "dependencies" ${INDEX} "source" "ref")
    string(JSON DEP_SHALLOW  GET "${DEPS_JSON}" "dependencies" ${INDEX} "source" "shallow")

    # If the DEP_CMAKE_OPTION field is non-empty, check the CMake option.
    if(NOT DEP_CMAKE_OPTION STREQUAL "")
        # If the named option variable does not exist in the CMake cache, warn and proceed.
        if(NOT DEFINED ${DEP_CMAKE_OPTION})
            message(FATAL_ERROR "Loading dependency ${DEP_NAME} at index ${INDEX} references unknown CMake option '${DEP_CMAKE_OPTION}'. Aborting.")
        else()
            # If the named option is defined but set to ON, skip this dependency.
            if(${${DEP_CMAKE_OPTION}})
                message(STATUS "Skipping FetchContent-step for dependency ${DEP_NAME} because CMake option '${DEP_CMAKE_OPTION}' is '${${DEP_CMAKE_OPTION}}'")
                continue()
            endif()
        endif()
    endif()

    message(STATUS "Configuring dependency: ${DEP_NAME} (Type: ${DEP_TYPE}, SHALLOW: ${DEP_SHALLOW})")

    #If Dep_SHALLOW is not a valid boolean, set it to FALSE and print a warning.
    if(NOT DEP_SHALLOW STREQUAL "TRUE" AND NOT DEP_SHALLOW STREQUAL "FALSE")
        message(STATUS "Invalid value for 'shallow' in dependency ${DEP_NAME}: '${DEP_SHALLOW}'. Expected 'true' or 'false'. Defaulting to 'false'.")
        set(DEP_SHALLOW FALSE)
    endif()

    # 2. Branch logic based on "type"
    if(DEP_TYPE STREQUAL "git")
        FetchContent_Declare(
            ${DEP_NAME}
            GIT_REPOSITORY ${DEP_URL}
            GIT_TAG        ${DEP_REF}
            GIT_PROGRESS   TRUE
            GIT_SHALLOW    ${DEP_SHALLOW}
        )
    else()
        message(WARNING "Unknown dependency type '${DEP_TYPE}' for ${DEP_NAME}")
    endif()

    # 3. Populate the dependency
    FetchContent_MakeAvailable(${DEP_NAME})
endforeach()
