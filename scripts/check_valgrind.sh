#!/bin/bash

#  This file is part of t8code.
#  t8code is a C library to manage a collection (a forest) of multiple
#  connected adaptive space-trees of general element classes in parallel.
#
#  Copyright (C) 2025 the developers
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

# Check that an argument is given and that the argument is a file
# Check if argument given
if [ ${1-x} = x ]
then
  echo ERROR: Need to provide a file as first argument.
  echo $usage
  exit 1
fi

# Check if first argument is a file and store it in variable
if [ -f "$1" ]
  then
  FILE="$1"
else
  # Try from folder above.
  if [ -f "../$1" ]
  then
    FILE="../$1"
  else
    echo "ERROR: Non existing file: $1"
    echo $usage
    exit 1
  fi
fi

# Write valgrind output to variable OUTPUT_FILE
OUTPUT_FILE="valgrind-output.log"
# Set valgrind flags.
VALGRIND_FLAGS="--leak-check=full --track-origins=yes --read-var-info=yes --trace-children=yes"
VALGRIND_FLAGS="$VALGRIND_FLAGS --show-leak-kinds=all --read-inline-info=yes --errors-for-leak-kinds=all"
VALGRIND_FLAGS="$VALGRIND_FLAGS --expensive-definedness-checks=yes --gen-suppressions=all --redzone-size=16"
VALGRIND_FLAGS="$VALGRIND_FLAGS --track-fds=yes"

# Run valgrind on given file with flags and write output to OUTPUT_FILE.
valgrind $VALGRIND_FLAGS "${FILE}" 2>"${OUTPUT_FILE}"

# Parse valgrind output
declare -a VALGRIND_RULES=(
        "^==.*== .* bytes in .* blocks are definitely lost in loss record .* of .*$"
        "^==.*== .* bytes in .* blocks are still reachable in loss record .* of .*$"
        "^==.*== Invalid .* of size .*$"
        "^==.*== Open file descriptor .*: .*$"
        "^==.*== Invalid free() / delete / delete\[\] / realloc()$"
        "^==.*== Mismatched free() / delete / delete \[\].*$"
        "^==.*== Syscall param .* points to uninitialised byte(s).*$"
        "^==.*== Source and destination overlap in .*$"
        "^==.*== Argument .* of function .* has a fishy (possibly negative) value: .*$"
        "^==.*== .*alloc() with size 0$"
        "^==.*== Invalid alignment value: .* (should be power of 2)$"
    )
report_id=1
status=0
error=""
kind="error"

while IFS= read -r line; do
    if [[ "${error}" != "" ]]; then
        if [[ $(echo "${line}" | grep '^==.*== $') ]]; then
            echo "::${kind} title=Valgrind Report '${FILE}' (${report_id})::${error}"
            report_id=$(( $report_id + 1 ))
            error=""
            status=1
        else
            error="${error}%0A${line}"
        fi
    fi
    for rule in "${VALGRIND_RULES[@]}"; do
        if [[ $(echo "${line}" | grep "${rule}") ]]; then
            error="${line}"
            break
        fi
    done
done < "${OUTPUT_FILE}"
rm -f "${OUTPUT_FILE}"
exit "${status}"
