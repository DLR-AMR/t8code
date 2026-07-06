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

# A CMakeLists.txt file showing a way of detecting an installed instance of T8CODE for use in program test.cxx.
# This example does not rely on what is being declared in files installed by T8CODE.

cmake_minimum_required(VERSION 3.22)
include(CheckIncludeFileCXX)
include(CheckLibraryExists)
project(t8test DESCRIPTION "t8code test program (needs T8TEST_ROOT to contain relevant include and lib directories and CMAKE_PREFIX_PATH to point to T8CODE's package files are.)" LANGUAGES CXX)
option(T8TEST_ROOT "Custom directory containing t8code's include, lib, ..." "")
find_library(T8CODE_LIBRARY t8 REQUIRED) # one may have done this in the if conditional construct after check_library_exists
CHECK_LIBRARY_EXISTS(t8 t8_init "${T8TEST_ROOT}/lib" HAVE_T8) # check_library_exists looks for a library and sets one variable accordingly
if(NOT HAVE_T8)
    message(FATAL_ERROR "no t8 library found!")
else()
    message("t8 library found in ${T8TEST_ROOT}/lib")
endif()
find_path(SC_INCLUDE_DIR sc.h REQUIRED)
find_path(T8CODE_INCLUDE_DIR t8.h REQUIRED)
add_executable(test_b test.cxx)
target_include_directories(test_b PUBLIC "${SC_INCLUDE_DIR}" )
target_include_directories(test_b PUBLIC "${T8CODE_INCLUDE_DIR}" )
target_link_directories(test_b PUBLIC "${T8TEST_ROOT}/lib" )
target_link_libraries(test_b PRIVATE t8 sc)
add_custom_target(tests_b COMMAND ./test_b)
