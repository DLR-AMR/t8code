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

# This shell script will parse all files given in the variable SEARCH_PATH
# for C-style TODO comments. That is comments starting with "TODO" and
# ending with "*/"
#
# The script then prints out all these comments together with file and line.

PROJECT_DIR=".."

for folder in src example test
do
  SEARCH_PATH="$SEARCH_PATH $PROJECT_DIR/$folder/*.h"
  SEARCH_PATH="$SEARCH_PATH $PROJECT_DIR/$folder/*.c"
  SEARCH_PATH="$SEARCH_PATH $PROJECT_DIR/$folder/*.hxx"
  SEARCH_PATH="$SEARCH_PATH $PROJECT_DIR/$folder/*.cxx"
  SEARCH_PATH="$SEARCH_PATH $PROJECT_DIR/$folder/*/*.h"
  SEARCH_PATH="$SEARCH_PATH $PROJECT_DIR/$folder/*/*.c"
  SEARCH_PATH="$SEARCH_PATH $PROJECT_DIR/$folder/*/*.hxx"
  SEARCH_PATH="$SEARCH_PATH $PROJECT_DIR/$folder/*/*.cxx"
done
grep TODO $SEARCH_PATH -n -A20 | awk '/TODO/,/*\//'
