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

TYPOS=`which typos 2> /dev/null`
if [ -z "$TYPOS" ]
then
  echo "ERROR: typos not found."
  echo "Please install typos."
  echo "See https://github.com/crate-ci/typos#install"
  exit 1
fi

repo_main_dir=`git rev-parse --show-toplevel`
(cd $repo_main_dir && typos)

echo "This script will correct all previously listed typos."
echo
read -p "Are you sure? ('Y' or 'y' to continue) " -n 1 -r
echo
if [[ $REPLY =~ ^[Yy]$ ]]
then
  echo Correcting all typos...
  (cd $repo_main_dir && typos -w)
  echo done.

  echo
  echo "Changed files according to git diff:"
  changed_files=``
  for file in `(cd $repo_main_dir && git diff --name-only)`
  do
    # only check existing files, this is necessary since if we rename or delete
    # a file it is added to the committed files and we thus would try to indent a
    # nonexisting file.
    if [ ! -e $repo_main_dir/$file ]
    then
      continue
    fi
    # We only indent .c, .cxx, .h and .hxx files
    #-a ${file: -2} != ".h" ]
    FILE_ENDING="${file##*.}"
    if [ $FILE_ENDING = "c" -o $FILE_ENDING = "h" -o $FILE_ENDING = "cxx" -o $FILE_ENDING = "hxx" ]
    then
      echo $file
      changed_files="$changed_files $file"
    fi
  done
  echo
  echo "Should also all changed (according to git diff) files be indented?"
  read -p "('Y' or 'y' to continue) " -n 1 -r
  echo
  if [[ $REPLY =~ ^[Yy]$ ]]
  then
    echo Indenting all files...
    (cd $repo_main_dir && ./scripts/t8indent.sh $changed_files)
    echo done.
  else
    echo Aborted.
  fi
else
  echo Aborted.
fi
