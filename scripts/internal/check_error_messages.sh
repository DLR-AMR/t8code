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
found_faulty_message=false

# Check if a warning or error message is not capitalized
while IFS=: read -r line_number line; do
    if ! echo "$line" | grep -qE "[\"'](WARNING|ERROR)"; then
        printf '\033[38;5;196m%s\n%s\n%s\n\033[0m\n' \
               "Incorrect error/warning message found in $file_path" \
               "Line $line_number: $line" \
               "Please use 'ERROR' or 'WARNING' instead."
        found_faulty_message=true
    fi
done < <(grep -nP '["`](?i:(warning|error))' "$file_path")

if [ "$found_faulty_message" == true ]; then
    exit 1
else
    exit 0
fi
