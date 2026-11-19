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

set(T8CODE_POST_INSTALL_TEST_SOURCE_DIR ${PROJECT_SOURCE_DIR}/test/post_install)
set(T8CODE_POST_INSTALL_TEST_BINARY_DIR ${PROJECT_BINARY_DIR}/post_install_tests)
set(T8CODE_POST_INSTALL_TEST_CMAKE_ARGS
  -DCMAKE_PREFIX_PATH=${CMAKE_INSTALL_PREFIX}
)

if(CMAKE_BUILD_TYPE)
  list(APPEND T8CODE_POST_INSTALL_TEST_CMAKE_ARGS -DCMAKE_BUILD_TYPE=${CMAKE_BUILD_TYPE})
endif()

if(DEFINED SC_DIR AND NOT "${SC_DIR}" STREQUAL "")
  list(APPEND T8CODE_POST_INSTALL_TEST_CMAKE_ARGS -DSC_DIR=${SC_DIR})
endif()

if(DEFINED P4EST_DIR AND NOT "${P4EST_DIR}" STREQUAL "")
  list(APPEND T8CODE_POST_INSTALL_TEST_CMAKE_ARGS -DP4EST_DIR=${P4EST_DIR})
endif()

add_custom_target( post_install_tests
  COMMAND ${CMAKE_COMMAND} -E rm -rf ${T8CODE_POST_INSTALL_TEST_BINARY_DIR}
  COMMAND ${CMAKE_COMMAND} -E make_directory ${T8CODE_POST_INSTALL_TEST_BINARY_DIR}
  COMMAND ${CMAKE_COMMAND}
          -G ${CMAKE_GENERATOR}
          -S ${T8CODE_POST_INSTALL_TEST_SOURCE_DIR}
          -B ${T8CODE_POST_INSTALL_TEST_BINARY_DIR}
          ${T8CODE_POST_INSTALL_TEST_CMAKE_ARGS}
  COMMAND ${CMAKE_COMMAND} --build ${T8CODE_POST_INSTALL_TEST_BINARY_DIR}
  COMMAND ${CMAKE_CTEST_COMMAND} --test-dir ${T8CODE_POST_INSTALL_TEST_BINARY_DIR} --output-on-failure
  COMMENT "Build and run the installed t8code post-installation sanity test"
  VERBATIM
)
