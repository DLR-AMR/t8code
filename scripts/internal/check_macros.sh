#!/bin/bash

# This file is part of t8code.
# t8code is a C library to manage a collection (a forest) of multiple
# connected adaptive space-trees of general element classes in parallel.
#
# Copyright (C) 2026 the developers
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

if [[ "$file_path" -ef "src/t8_misc/t8_with_macro_error.h" ]]
then
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

found_macros=false

# Check for #ifdef and T8_WITH
while IFS=: read -r line_number line; do
    printf '\033[38;5;196m%s\n%s\n\033[0m' \
           "Incorrect macro found in $file_path" \
           "Line $line_number: $line"
    # Print message for #ifdef
    if echo $line | grep -q '#ifdef'; then
        printf '\033[38;5;196m%s\n\033[0m' \
               "Please do not use #ifdef with T8_ENABLE macros. Use #if instead."
    fi
    # Print message for T8_WITH
    if echo $line | grep -q 'T8_WITH'; then
        printf '\033[38;5;196m%s\n\033[0m' \
               "T8_WITH macros are not longer supported. Use T8_ENABLE instead."
    fi
    printf "\n"
    found_macros=true
done < <(grep -nE '(#ifdef T8_ENABLE|T8_WITH_)' "$file_path")

if [ "$found_macros" == true ]; then
    exit 1
else
    exit 0
fi
