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
# This file lists all binary paths to test files except the tests for the api.
# The script can be used to run the check_valgrind script for every test binary. 
# The CMakeLists.txt file in the test folder is used to generate the list of binary paths.
# The paths are relative paths assuming an execution from the current folder to a build folder in the main t8code folder.
#

cmake_file="../test/CMakeLists.txt"

# Check that path to cmake_file is correct.
if [ ! -f "$cmake_file" ]; then
  echo "CMakeLists.txt not found!"
  exit 1
fi

# List with test binary paths.
test_bin_paths=""
name=""

# Iterate over the lines of the cmake file. 
while IFS= read -r line; do
  # Consider line if it contains “add_t8_test(”.
  if [[ "$line" == *"add_t8_test("* ]]; then
    # Extract the name of the executable.
    [[ "$line" =~ NAME[[:space:]]+([a-zA-Z0-9_]+) ]]
    name="${BASH_REMATCH[1]}"

    # Extract the subfolder of the executable.
    [[ "$line" =~ SOURCES[[:space:]]+([a-zA-Z0-9_/.-]+) ]]
    related_source="${BASH_REMATCH[1]}"
    # Check if the first argument is 't8_gtest_main.cxx'.
    if [[ "$related_source" == "t8_gtest_main.cxx" ]]; then
      # Extract the second argument if the first one is 't8_gtest_main.cxx'.
      if [[ "$line" =~ SOURCES[[:space:]]+[a-zA-Z0-9_/.-]+[[:space:]]+([a-zA-Z0-9_/.-]+) ]]; then
        related_source="${BASH_REMATCH[1]}"
      fi
    fi
    # Add folder location if file is in a subfolder.
    if [[ "$related_source" =~ ^(.*)/ ]]; then
      name="${BASH_REMATCH[1]}/${name}"
    fi
    # Set correct path and add it to the list.
    name="../build/test/${name}"
    # Do not include api in the binary paths as we do not want to check them with valgrind.
    if ! [[ "$name" =~ "/api/" ]]; then
      test_bin_paths="${test_bin_paths} ${name}"
    fi
  fi
done < "$cmake_file"

echo $test_bin_paths