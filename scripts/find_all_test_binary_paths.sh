#!/bin/bash

#  This file is part of t8code.
#  t8code is a C library to manage a collection (a forest) of multiple
#  connected adaptive space-trees of general element classes in parallel.
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

#
# This file lists all paths to test binaries that exist in the build/test directory.
# The script can be used to run the check_valgrind script for every test binary.
# The paths are relative paths assuming an execution from the test/ folder in the build directory.
#

build_test_directory="../build/test/"

# Check that path to test folder in build directory is correct.
if [ ! -d "$build_test_directory" ]; then
  echo "Directory build/test/ not found!"
  exit 1
fi

# Find all executables in the build/test/ directory (that do not have a .so file ending) 
# and store the relative paths to the build_test_directory with leading "./".
test_bin_paths=$(find "$build_test_directory" -type f -executable -not -name "*.so" -exec realpath --relative-to="$build_test_directory" {} \; | sed 's|^|./|')

echo $test_bin_paths
