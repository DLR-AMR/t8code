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
# We have to use the -i argument, so that clang-format directly alters the
# files instead of printing the changes to stdout. The --style=file
# arguments tells clang-format to look for a *.clang-format file.
FORMAT_OPTIONS="-i --style=file"

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

#
#  Parsing of input files and throwing out files to be ignored
#
# Read all lines from the IGNORE_FILE 
# that are not empty and are not comments (i.e. start with '#').
# Determine base directory of git repo
GIT_REPO_PATH=$(git rev-parse --show-toplevel)

IGNORE_FILE=${GIT_REPO_PATH}/scripts/t8indent_ignore.sh
files_to_ignore=()
while read line; do
    if [[ ${line:0:1} != "#" ]] && [[ $line != "" ]]
    then
        files_to_ignore+=("$line")
    fi
done <$IGNORE_FILE

#
# Check if first argument is "NO_CHANGE", if so
# the file content is not changed.
#
OUTFILE_OPTION=
NO_CHANGE=FALSE
if [[ $1 == "NO_CHANGE" ]]
then
  shift # Removes first argument from $@ list
  NO_CHANGE=TRUE
fi


# Iterate over all arguments and throw
# aways those filenames that we should ignore.
# Also check if suffix is ".c" ".cxx" ".h" or ".hxx"
for arg in "$@"
do
  FILE_SUFFIX="${arg##*.}"
  if ! [ $FILE_SUFFIX = "c" -o $FILE_SUFFIX = "h" -o $FILE_SUFFIX = "cxx" -o $FILE_SUFFIX = "hxx" ]
  then
    echo "ERROR: File "$arg" does not have valid suffix (.c .h .cxx .hxx)."
    exit 1
  fi

  ignore_arg=0
  # Iterate over each ignore filename
  for ignore_file in "${files_to_ignore[@]}"
    do
    if [[ "$arg" -ef "${GIT_REPO_PATH}/$ignore_file" ]]
    then 
      # arg matches and will be ignored
      echo The file \"$arg\" will be ignored by indentation as specified in \"$IGNORE_FILE\".
      ignore_arg=1
    fi
  done
  # Now add all non-ignored files to a new argument array
  if [[ $ignore_arg == 0 ]]
  then
    newargs+=("$arg")
  fi
done


for arg in "$@" ; do
  if [ "x$arg" == "x-o" ]; then
    WANTSOUT=1
  fi
done
if [ -z "$WANTSOUT" ]; then
  for NAME in "${newargs[@]}" ; do
    if [[ $NO_CHANGE == "TRUE" ]]
    then
      $FORMAT $FORMAT_OPTIONS "$NAME" 2>&1
    else
      $FORMAT $FORMAT_OPTIONS "$NAME"
    fi
    status=$?
  done
else
  if [[ $NO_CHANGE == "TRUE" ]]
  then
    $FORMAT $FORMAT_OPTIONS ${newargs[@]} 2>&1
  else
    $FORMAT $FORMAT_OPTIONS ${newargs[@]}
  fi  
  status=$?
fi

# If the file content was not change, the return
# value determines whether or not the file was
# indented.
if [[ $NO_CHANGE == "TRUE" ]]
then
  exit $status
fi