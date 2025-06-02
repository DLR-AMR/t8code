#! /bin/bash

# This file is part of t8code.
# t8code is a C library to manage a collection (a forest) of multiple
# connected adaptive space-trees of general element classes in parallel.
#
# Copyright (C) 2023 the developers
#
# t8code is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 2 of the License, or
# (at your option) any later version.
#
# t8code is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with t8code; if not, write to the Free Software Foundation, Inc.,
# 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.

GIT_REPO_PATH=$(git rev-parse --show-toplevel)

INDENT_SCRIPT=${GIT_REPO_PATH}/scripts/t8indent.sh

usage="$0 [FILE_TO_INDENT]"

# Check if first argument given
if [ ${1-x} = x ]
then
  echo ERROR: Need to provide a file as first argument.
  echo $usage
  exit 1
fi

# Check if first argument is a file and store it in variable
if [ -f "$1" ]
  then
  file="$1"
else
  # Try from folder above.
  if [ -f "../$1" ]
  then
    file="../$1"
  else
    echo "ERROR: Non existing file: $1"
    echo $usage
    exit 1
  fi
fi

#
# Check if the file is indented
#
$INDENT_SCRIPT NO_CHANGE $file
status=$?
if [ $status != 0 ]
then
  echo $file is not indented.
  echo
fi
exit $status

