#!/bin/bash

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

cd "$(git rev-parse --show-toplevel)"

git submodule init
git submodule update

# Remove build directory if it exists
if [ -d build ]; then
    rm -rf build
fi

# Create the build directory
mkdir build

# Navigate into the build directory
cd build

if [ "$READTHEDOCS" = "True" ]; then
    DOXYFILE_PATH="../doc/Doxyfile.in"
    if [ ! -f "$DOXYFILE_PATH" ]; then
        echo "Error: $DOXYFILE_PATH does not exist or is not a regular file."
        exit 1
    fi
    echo "Configuring Doxygen for ReadTheDocs: Excluding source files."
    # Exclude source files from documentation
    echo "EXCLUDE_PATTERNS += *.c *.cc *.cpp *.cxx" >> "$DOXYFILE_PATH"
fi

cmake .. -DT8CODE_BUILD_DOCUMENTATION=ON -DT8CODE_BUILD_DOCUMENTATION_SPHINX=ON -DT8CODE_ENABLE_MPI=OFF

# Return to the parent directory
cd ..