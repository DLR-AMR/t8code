#!/bin/bash

#  This file is part of t8code.
#  t8code is a C library to manage a collection (a forest) of multiple
#  connected adaptive space-trees of general element classes in parallel.
#
#  Copyright (C) 2023 Johannes Markert <johannes.markert@dlr.de>
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


# Search and remove all '#if 0' pragmas from C/C++ source files in
# 'src', 'tutorials', and 'example' directories.
#
# Uses 'unidef', a non-standard commandline tool: https://dotat.at/prog/unifdef/
#
# On Ubuntu: sudo apt-get install unidef

#
# Usage: ./scripts/remove-all-if-zero-pragmas.sh
#

grep -r -i -l -E '^\s*#\s*if\s+0' src/ tutorials/ example/ test/ | xargs unifdef -k -m
