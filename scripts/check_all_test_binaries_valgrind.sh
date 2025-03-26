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
# This script performs a valgrind check on each test binary given by find_all_test_binary_paths.sh.
# The valgrind check is done by the check_valgrind.sh script.
# The script returns 1 if an error is found and 0 otherwise. 
# This script must be executed from the scripts/ folder.
# It is assumed that the build folder ../build/test/ with the correct test binaries exists.
#

# Script must be executed from the scripts/ folder.
if [ `basename $PWD` != scripts ]; then
  if [ -d scripts ]; then
    # The directory stack is automatically reset on script exit.
    pushd scripts/ > /dev/null
  else
    echo "ERROR: scripts/ directory not found."
    exit 1
  fi
fi

# Find all test binary paths.
test_bin_paths=`bash ./find_all_test_binary_paths.sh`
num_paths=$(echo $test_bin_paths | wc -w)

# This is necessary because some tests use test files specified by relative paths. 
# These tests only work when run from the build/test/ directory. 
if [ -d ../build/test ]; then
  # The directory stack is automatically reset on script exit.
  pushd ../build/test/ > /dev/null
else
  echo "ERROR: Couldn't find a the directory ../build/test/."
  exit 1
fi

status=0
counter=0
valgrind_suppressions_file=../../scripts/valgrind_suppressions_file.supp

for bin_path in $test_bin_paths; do
  echo "[$counter/$num_paths] Valgrind check of $bin_path..."
  # Run check_valgrind script for each test binary.
  bash ../../scripts/check_valgrind.sh $bin_path $valgrind_suppressions_file 2>&1
  status=$?
  # If status is not 0, an error occurred.
  if test $status -ne 0; then
    echo "Error occurred during the valgrind check of $bin_path."
    exit 1
  fi
  counter=$(( $counter + 1 ))
done

echo "Valgrind found no errors in the test executables."
exit 0