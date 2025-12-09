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

find_package( Git REQUIRED )

    if(EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/version.txt")
        file(READ "${CMAKE_CURRENT_SOURCE_DIR}/version.txt" VERSION_CONTENT)
        # Extract the version number 
        string(REGEX MATCH "Version ([0-9])\.([0-9]+)\.([0-9]+)" VERSION_MATCH "${VERSION_CONTENT}" )
        if (VERSION_MATCH)
            # The version number will be in ${CMAKE_MATCH_1}
            set(T8CODE_VERSION_NUMBERS "${CMAKE_MATCH_1}.${CMAKE_MATCH_2}.${CMAKE_MATCH_3}")
            set(T8CODE_VERSION_RAW "${CMAKE_MATCH_1}.${CMAKE_MATCH_2}.${CMAKE_MATCH_3}")
            message(STATUS "Extracted Version: ${T8CODE_VERSION_NUMBERS}")
        else()
            message(WARNING "Version number not found in version.txt")
        endif()
    else()
        message(WARNING "Version information not found")
    endif()


execute_process( COMMAND echo ${T8CODE_VERSION_NUMBERS}
                COMMAND cut -d. -f1
                OUTPUT_VARIABLE T8CODE_VERSION_MAJOR
                OUTPUT_STRIP_TRAILING_WHITESPACE )
execute_process( COMMAND echo ${T8CODE_VERSION_NUMBERS}
                COMMAND cut -d. -f2
                OUTPUT_VARIABLE T8CODE_VERSION_MINOR
                OUTPUT_STRIP_TRAILING_WHITESPACE )
execute_process( COMMAND echo ${T8CODE_VERSION_NUMBERS}
                COMMAND cut -d. -f3
                OUTPUT_VARIABLE T8CODE_VERSION_PATCH
                OUTPUT_STRIP_TRAILING_WHITESPACE )
set(T8_VERSION ${T8CODE_VERSION_NUMBERS} CACHE INTERNAL "")


