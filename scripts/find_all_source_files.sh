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

# This script lists all .c .h .cxx and .hxx
# files in t8code's src/ example/ and test/ subfolders.
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

# All valid file suffixes.
# Separated by '|' in order to be directly used
# as a regex in the find command.

suffixes="c|cxx|h|hxx"

if find --version >/dev/null 2>&1; then
    FIND=find
else
    echo "GNU find not found, trying gfind..."
    if gfind --version >/dev/null 2>&1; then
        FIND=gfind
    else
        echo "Error: GNU find not found."
        exit 1
    fi
fi

# Find all files with the appropriate suffix in the
# src/, example/, test/, tutorials/, benchmark/, and /api subfolders.
files=`$FIND ../src ../example ../test ../tutorials ../benchmarks ../api -regextype egrep -iregex ".*\.($suffixes)"`

echo $files
