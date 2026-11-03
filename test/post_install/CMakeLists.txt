#  This file is part of t8code.
#  t8code is a C library to manage a collection (a forest) of multiple
#  connected adaptive space-trees of general element types in parallel.
#
#  Copyright (C) 2026 the developers
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

# Post-install tests for t8code.
#
# This is a small downstream project that consumes an *installed* t8code via
# find_package(), the same way an external project would. Each test case lives
# in its own subdirectory and builds its own executable.
#
# Usually run via the "post_install_tests" target of the main t8code build.
# Manual use against an install prefix:
#   cmake -S <this dir> -B build-post-install -DCMAKE_PREFIX_PATH=<install prefix>
#   cmake --build build-post-install
#   ctest --test-dir build-post-install --output-on-failure

cmake_minimum_required( VERSION 3.20 )

project( T8CODE_POST_INSTALL_TESTS LANGUAGES C CXX )

# Must be called in the top-level directory, otherwise ctest finds no tests.
enable_testing()

find_package( T8CODE CONFIG REQUIRED )

# Test cases
add_subdirectory( sanity_check )
add_subdirectory( header_checks )