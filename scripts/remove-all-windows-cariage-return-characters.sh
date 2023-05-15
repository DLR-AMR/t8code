#!/bin/bash

# This file is part of t8code.  t8code is a C library to manage a collection (a
# forest) of multiple connected adaptive space-trees of general element classes
# in parallel.
#
# Copyright (C) 2023 Johannes Markert <johannes.markert@dlr.de>
#
# t8code is free software; you can redistribute it and/or modify it under the
# terms of the GNU General Public License as published by the Free Software
# Foundation; either version 2 of the License, or (at your option) any later
# version.
#
# t8code is distributed in the hope that it will be useful, but WITHOUT ANY
# WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR
# A PARTICULAR PURPOSE.  See the GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License along with
# t8code; if not, write to the Free Software Foundation, Inc., 51 Franklin
# Street, Fifth Floor, Boston, MA 02110-1301, USA.

# Usage
#
# Search and remove all Windows(TM) cariage return characters '\r' from C/C++
# source files in source code directories. Usually, these
# characters are visible as '^M' characters and are introduced by misconfigured
# source code editors on Windows(TM).

grep -R -l $'\r' src tutorials example test | xargs -r -l1 perl -i -p -e 's/\r//g'
