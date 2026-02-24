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
# With "--ntasks=[NUMBER]", you can provide the number of processes to use with MPI for parallel tests (default is 1).
#

USAGE="\nUSAGE: This script executes valgrind in parallel on each test binary available. Use the syntax \n
$0 [TEST_BINARY_PATH] --ntasks=[NUM_TASKS]\n
Providing the number of parallel processes to use with MPI for parallel tests is optional.\n"

# Check directory exists
if [ -z "$1" ]; then
  echo "ERROR: Need to provide a directory as first argument."
  echo -e "$USAGE"
  exit 1
fi

# Check if it is a directory
if [ -d "$1" ]; then
  TEST_BINARY_PATH="$1"
else
  echo "ERROR: Directory does not exist: $1"
  echo -e "$USAGE"
  exit 1
fi

#convert to abspath
TEST_BINARY_ABSPATH=$(realpath $TEST_BINARY_PATH)

# Check if a number of processes is provided. If not, set to 1.
num_procs=1
for arg in "$@"; do
  if [[ "$arg" == --ntasks=* ]]; then
    ntasks_val="${arg#--ntasks=}"
    if [[ "$ntasks_val" =~ ^[0-9]+$ ]]; then
      num_procs="$ntasks_val"
    else
      echo "ERROR: --ntasks value '$ntasks_val' is not a valid number."
      echo -e "$USAGE"
      exit 1
    fi
  fi
done

# Script must be executed from the scripts/ folder.
if [ `basename $PWD` != scripts ]; then
  if [ -d scripts ]; then
    # The directory stack is automatically reset on script exit.
    pushd scripts/ > /dev/null
  else
    echo "ERROR: scripts/ directory not found."
    echo -e "$USAGE"
    exit 1
  fi
fi

# Find all test binary paths.
test_bin_paths=$(bash ./find_all_test_binary_paths.sh "$TEST_BINARY_ABSPATH")
status=$?

if [ $status -ne 0 ]; then
  echo "$test_bin_paths"
  echo "Failed to collect test binaries."
  exit $status
fi

num_paths=$(echo $test_bin_paths | wc -w)

# This is necessary because some tests use test files specified by relative paths.
# These tests only work when run from the build/test/ directory.
if [ -d ../build/test ]; then
  # The directory stack is automatically reset on script exit.
  pushd ../build/test/ > /dev/null
else
  echo "ERROR: Couldn't find a the directory ../build/test/."
  echo -e "$USAGE"
  exit 1
fi

status=0
counter=0
valgrind_suppressions_file=../../scripts/valgrind_suppressions_file.supp

# First run all serial tests in parallel.
serial_count=$(echo $test_bin_paths | tr ' ' '\n' | grep -c '_serial')
echo "Running valgrind checks on $serial_count serial test binaries in parallel..."
for bin_path in $test_bin_paths; do
  if [[ $bin_path == *"_serial"* ]]; then
    counter=$(( $counter + 1 ))
    # Run check_valgrind script for each test binary. The & at the end allows parallel execution.
    bash ../../scripts/check_valgrind.sh $bin_path --supp=$valgrind_suppressions_file --ntasks=1 2>&1 &
    status=$?
    # If status is not 0, an error occurred.
    if test $status -ne 0; then
      echo "Error occurred during the valgrind check of $bin_path."
      kill $(jobs -p)
      exit 1
    fi
  fi
done
# Wait until all serial valgrind checks are done.
wait
echo "Valgrind found no errors in the ${counter} serial test executables."
echo "Check remaining parallel tests."
for bin_path in $test_bin_paths; do
  if [[ $bin_path == *"_parallel"* ]]; then
    counter=$(( $counter + 1 ))
    echo -n "[$counter/$num_paths] "
    # Run check_valgrind script for each test binary.
    bash ../../scripts/check_valgrind.sh $bin_path --supp=$valgrind_suppressions_file --ntasks=$num_procs 2>&1
    status=$?
    # If status is not 0, an error occurred.
    if test $status -ne 0; then
      echo "Error occurred during the valgrind check of $bin_path."
      exit 1
    fi
  fi
done

echo "Valgrind found no errors in the test executables."
exit 0