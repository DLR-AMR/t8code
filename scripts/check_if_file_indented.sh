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

# We use clang-format to indent our code. There are different base styles
# available, but we use a modification. The options for this modification
# are located in the .clang-format file.
# We have to use the --dry-run argument, so that clang-format does not change
# anything. -Werror produces an error if the file is not indented cerrctly.
# The --style=file arguments tells clang-format to look for a *.clang-format file.
FORMAT_OPTIONS="--dry-run --Werror --style=file"

# Required version of the clang format program.
REQUIRED_VERSION_MAJOR="17"
REQUIRED_VERSION_MINOR="0"
REQUIRED_VERSION_STRING="${REQUIRED_VERSION_MAJOR}.${REQUIRED_VERSION_MINOR}"

FORMAT=`which clang-format 2> /dev/null`

if [ -z "$FORMAT" ]
then
  # Exit if the spell checking script was not found
  echo "ERROR: clang-format not found."
  echo "Please install clang-format version ${REQUIRED_VERSION_STRING}."
  echo "See https://github.com/ssciwr/clang-format-wheel"
  exit 1
fi

CLANG_VERSION_STRING=`$FORMAT --version`

VERSION=`echo $CLANG_VERSION_STRING | cut -d " " -f 3`
MAJOR=`echo $VERSION | cut -d. -f1`
MINOR=`echo $VERSION | cut -d. -f2`
PATCH=`echo $VERSION | cut -d. -f3`

if [[ "$MAJOR" != "$REQUIRED_VERSION_MAJOR" || $MINOR != "$REQUIRED_VERSION_MINOR" ]]; then
  echo "Please install clang-format version $REQUIRED_VERSION_STRING"
  exit 1
fi

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

# Here is the actual checking part
#
# We only check files with .c, .h, .cxx, .hxx suffix
FILE_SUFFIX="${file##*.}"
if [ $FILE_SUFFIX = "c" -o $FILE_SUFFIX = "h" -o $FILE_SUFFIX = "cxx" -o $FILE_SUFFIX = "hxx" ]
  then
  $FORMAT $FORMAT_OPTIONS $file 2>&1
  status=$?
  if [ $status != 0 ]
  then
    echo $file is not indented.
    echo
  fi
  exit $status
else
  echo "ERROR: File does not have valid suffix (.c .h .cxx .hxx)."
fi
