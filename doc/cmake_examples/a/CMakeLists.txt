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

# A CMakeLists.txt example file showing a way of detecting an installed instance of T8CODE and using it in program test.cxx.
# This example relies on what is being declared in files installed by T8CODE.

cmake_minimum_required(VERSION 3.22.0)
project(t8test DESCRIPTION "t8code test program (requires CMAKE_PREFIX_PATH to point to T8CODE's package files are)" LANGUAGES CXX)
find_package(T8CODE)
add_executable(test_a test.cxx)
target_link_libraries(test_a PRIVATE T8CODE::T8)
add_custom_target(tests_a COMMAND ./test_a)
