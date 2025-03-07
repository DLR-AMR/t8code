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

#
# This script runs Valgrind on an input binary paths with specified memory leak detection flags. 
# The Valgrind output is parsed. If any errors are found, they are printed and the script exits with a status of 1.
# As a second argument, you can provide a path to a suppression file that is used by Valgrind to suppress certain errors.
#

# Check that an argument is given and that the argument is a file
# Check if argument given
if [ ${1-x} = x ]; then
  echo ERROR: Need to provide a file as first argument.
  exit 1
fi

# Check if first argument is a file and store it in variable
if [ -f "$1" ]; then
  FILE="$1"
else
  # Try from folder above.
  if [ -f "../$1" ]; then
    FILE="../$1"
  else
    echo "ERROR: Non existing file: $1"
    exit 1
  fi
fi

# Write valgrind output to variable OUTPUT_FILE.
OUTPUT_FILE="valgrind-output.log"
# Set valgrind flags.
VALGRIND_FLAGS="--leak-check=full --track-origins=yes \
    --trace-children=yes --show-leak-kinds=definite,indirect,possible \
    --errors-for-leak-kinds=definite,indirect,possible --gen-suppressions=all"
# There are some more flags that can be reasonable to use, e.g., for debugging reasons if you found an error.
# We used minimal flags for performance reasons.
# Further flags include (but of course are not limited to): --expensive-definedness-checks=yes --track-fds=yes
# For more detailed outputs: -read-var-info=yes --read-inline-info=yes
# Warning: --show-leak-kinds=all will find a lot of still reachable leaks. This is not necessarily a problem.

# Check if a second argument is provided. If yes, add the flag to incorporate the Valgrind suppression file.
if ! [ ${2-x} = x ]; then
  if [ -f "$2" ]; then
    VALGRIND_FLAGS="${VALGRIND_FLAGS} --suppressions=${2}"
  else
    echo "ERROR: If a second argument is provided, this must be a valid valgrind suppression file."
    exit 1
  fi
fi

# Run valgrind on given file with flags and write output to OUTPUT_FILE.
valgrind $VALGRIND_FLAGS "${FILE}" > /dev/null 2>"${OUTPUT_FILE}"

# Parse valgrind output.
declare -a VALGRIND_RULES=(
        "^==.*== .* bytes in .* blocks are definitely lost in loss record .* of .*$"
        "^==.*== .* bytes in .* blocks are indirectly lost in loss record .* of .*$"
        "^==.*== .* bytes in .* blocks are possibly lost in loss record .* of .*$"
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

while IFS= read -r line; do
  if [[ "${error}" != "" ]]; then
    # Error message of valgrind always end with a line ==.*== without any further information.
    # Only print if we collected every line of the error message.
    if [[ $(echo "${line}" | grep '^==.*== $') ]]; then
      echo "::Error found in valgrind report '${FILE}' (${report_id})::"
      echo -e "${error}"
      echo ""
      report_id=$(( $report_id + 1 ))
      error=""
      status=1
    else
      # Add to error message that is printed with the last line of the error.
      error="${error}\n${line}"
    fi
  fi
  for rule in "${VALGRIND_RULES[@]}"; do
    # Check if we found one of the errors defined in VALGRIND_RULES.
    if [[ $(echo "${line}" | grep "${rule}") ]]; then
      error="${line}"
      break
    fi
  done
  if [[ $(echo "${line}" | grep '^==.*== ERROR SUMMARY:') ]]; then
    if ! [[ $(echo "${line}" | grep '^==.*== ERROR SUMMARY: 0 ') ]]; then
      # Set status to 1 if an error was found that is not included in VALGRIND_RULES.
      status=1
    fi
    echo "${line}"
  elif [[ $(echo "${line}" | grep 'valgrind:.*: command not found') ]]; then
    echo "${line}"
    status=1
  fi
done < "${OUTPUT_FILE}"

rm -f "${OUTPUT_FILE}"
exit "${status}"