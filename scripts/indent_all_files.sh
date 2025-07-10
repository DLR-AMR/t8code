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

# This script indents all .c .h .cxx and .hxx
# files in t8code.
#
# Returns 0 if yes and not 0 if not.
#

#
# This script must be executed from the scripts/ folder.
#
if [ `basename $PWD` != scripts ]
then
  if [ -d scripts ]
  then
    # The directory stack is automatically reset on script exit.
    pushd scripts/ > /dev/null
  else
    echo ERROR: scripts/ directory not found.
    exit 1
  fi
fi

# Find all files with the appropriate suffix.
# Excluding the sc/ and p4est/ subfolders.
files=`./find_all_source_files.sh`

notallindented=0

INDENT=./t8indent.sh

echo This script will change all source files found in 
echo	$PWD/../src/
echo	$PWD/../example/
echo	$PWD/../test/
echo	$PWD/../tutorials/
echo	$PWD/../benchmarks/
echo	$PWD/../api/
echo
read -p "Are you sure? ('Y' or 'y' to continue)" -n 1 -r
echo
if [[ $REPLY =~ ^[Yy]$ ]]
then
  echo Indenting all files...
  for file in $files
  do
    # Find also give as directories,
    # so we ensure that $file is a proper
    # file before checking for indentation.
    if [ -f $file ]
    then
      $INDENT $file
      status=$?
      if test $status -ne 0 
      then
        echo "File $file is not indented."
        notallindented=1
      fi
    fi
  done
  
  if test $notallindented -eq 0
  then
    echo All files are indented.
  fi
  echo done.
else
  echo Aborted.
fi

exit $notallindented
