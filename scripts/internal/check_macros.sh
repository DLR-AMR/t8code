#!/bin/bash

# This file is part of t8code.
# t8code is a C library to manage a collection (a forest) of multiple
# connected adaptive space-trees of general element classes in parallel.
#
# Copyright (C) 2025 the developers
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

if [ "$#" -ne 1 ]; then
    echo "Usage: $0 <file_path>"
    exit 1
fi

file_path=$1

echo "$file_path"

if [[ "$file_path" -ef "src/t8_misc/t8_with_macro_error.h" ]]
then
  echo The file \"src/t8_misc/t8_with_macro_error.h\" will be ignored by the check_macros.sh script.
  exit 0
fi

#
# This script searches for lines containing a macro definition in the style of '#ifdef T8_ENABLE_'
# in the specified file and processes each matching line.
# It uses 'grep' to find all occurrences of '#ifdef T8_ENABLE_' in the file located
# at the path stored in the variable 'file_path'. The '-n' option with 'grep'
# ensures that the line numbers of the matching lines are included in the output.
# The output of 'grep' is then piped into a 'while' loop, which reads each line
# and splits it into the line number and the line content using ':' as the delimiter.
# Variables:
# - file_path: The path to the file to be searched.
# - line_number: The line number where the macro definition is found.
# - line: The content of the line where the macro definition is found.
#

found_macros=FALSE

# Check for #ifdef T8_ENABLE
while IFS=: read -r line_number line; do
    macro_name=$(echo "$line" | grep -o 'T8_ENABLE_[^ ]*')
    echo "Incorrect macro found in $file_path on line $line_number: $macro_name. Please use '#if T8_ENABLE_' instead."
    found_macros=TRUE
done < <(grep -n '#ifdef T8_ENABLE_' "$file_path")

# Check for #ifdef T8_WITH or #if T8_WITH
while IFS=: read -r line_number line; do
    macro_name=$(echo "$line" | grep -o 'T8_WITH_[^ ]*')
    echo "Incorrect macro found in $file_path on line $line_number: $macro_name. Please use '#if T8_ENABLE_' instead."
    found_macros=TRUE
done < <(grep -E -n '#if T8_WITH_|#ifdef T8_WITH' "$file_path")

if [ "$found_macros" = "TRUE" ]; then
    exit 1
else
    exit 0
fi
