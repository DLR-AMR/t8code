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

include(FindPackageHandleStandardArgs)

find_package(PythonInterp)
if(PYTHONINTERP_FOUND)
    get_filename_component(_PYTHON_DIR "${PYTHON_EXECUTABLE}" DIRECTORY)
    set(
        _PYTHON_PATHS
        "${_PYTHON_DIR}"
        "{_PYTHON_DIR}/bin"
        "{_PYTHON_DIR}/Scripts"
    )
endif()

#Look for an executable called sphinx-build
find_program(SPHINX_EXECUTABLE
             NAMES sphinx-build sphinx-build.exe sphinx-build.py
             HINTS ${_PYTHON_PATHS}
             DOC "Path to sphinx-build executable")


#Handle standard arguments to find_package like REQUIRED and QUIET
find_package_handle_standard_args(Sphinx
                                  "Failed to find sphinx-build executable"
                                  SPHINX_EXECUTABLE)