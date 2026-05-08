#!/bin/bash

#  This file is part of t8code.
#  t8code is a C library to manage a collection (a forest) of multiple
#  connected adaptive space-trees of general element classes in parallel.
#
#  Copyright (C) 2023 the developers
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


# This script checks whether all .c .h .cxx and .hxx
# files in t8code are properly indented.
#
# Returns 0 if yes and not 0 if not.
#

#
# This script must be executed from within the repository.
#
repo_main_dir=$(git rev-parse --show-toplevel 2>/dev/null)

if [ $? -ne 0 ]; then
  echo "ERROR: check_if_all_files_indented.sh was not called from inside the git repository."
  exit 1
fi

# Find all files with the appropriate suffix.
files=$($repo_main_dir/scripts/internal/find_all_source_files.sh) || {
  echo $files # return error message of find_all_source_files.sh
  echo "ERROR: find_all_source_files.sh returned exit code 1"
  exit 1
}

if [ -z "$files" ]; then
  echo "ERROR: find_all_source_files.sh returned nothing."
  exit 1
fi

notallindented=0
file_found=0
for file in $files
do
  if [ -f "$file" ]
  then
    file_found=1
    $repo_main_dir/scripts/check_if_file_indented.sh "$file"
    status=$?
    if test $status -ne 0
    then
      notallindented=1
    else
      echo "REMOVE THIS: File $file is indented."
    fi
  else
    echo "ERROR: find_all_source_files.sh returned a file which does not exist."
    exit 1
  fi
done

if test $file_found -eq 0
then
  echo Error: Could not find any source files.
  exit 1
fi

if test $notallindented -eq 0
then
  echo All files are indented.
  exit 0
fi

exit 1
