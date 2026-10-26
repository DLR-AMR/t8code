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

# Post-install header check for the installed t8code headers.
#
# For every public header, this generates:
#   - one translation unit that includes the header twice
#     (self-containment + include guard check)
# plus
#   - one translation unit that includes all headers together
#     (headers coexist: no conflicting declarations, macro clashes, order issues)
#   - a trivial main()
# and links everything into one executable, so non-inline definitions in
# headers show up as "multiple definition" link errors.
#
# The check is the build itself: if header_check compiles and links, it passes.
# T8CODE is found by the parent CMakeLists.txt.

# ---------------------------------------------------------------------------
# Configuration
# ---------------------------------------------------------------------------
set(HEADER_CHECK_PACKAGE       "T8CODE"      CACHE STRING "Package name passed to find_package()")
set(HEADER_CHECK_TARGET        "T8CODE::T8"  CACHE STRING "Imported target whose headers are checked")
# Matched against the header path relative to the include dir, e.g.
# "t8_foo.h", "t8_cmesh/t8_cmesh.h", "t8_schemes/sub/bar.hxx".
# "^t8(_|\\.h$)" = everything starting with "t8_" (files and folders) plus t8.h
set(HEADER_CHECK_INCLUDE_REGEX "^t8(_|\\.h$)" CACHE STRING "Regex selecting this library's headers")
set(HEADER_CHECK_EXCLUDE_REGEX "/detail/"     CACHE STRING "Regex of headers to skip (empty = check all)")
set(HEADER_CHECK_CXX_STANDARD  "20"           CACHE STRING "Minimum C++ standard the headers must support")

# Optional allowlist for an incremental rollout: if this file exists, ONLY the
# headers listed in it are checked (one path per line, relative to the include
# dir, '#' starts a comment line). Delete the file to check all headers.
# It is a normal variable, not a cache entry, so editing it takes effect on
# the next configure without clearing the cache.
set(HEADER_CHECK_LIST_FILE "${CMAKE_CURRENT_SOURCE_DIR}/headers_to_check.txt")

# ---------------------------------------------------------------------------
# Find the installed package (normally already done by the parent project)
# ---------------------------------------------------------------------------
if(NOT TARGET ${HEADER_CHECK_TARGET})
  find_package(${HEADER_CHECK_PACKAGE} CONFIG REQUIRED)
endif()

if(NOT TARGET ${HEADER_CHECK_TARGET})
  message(FATAL_ERROR "Package ${HEADER_CHECK_PACKAGE} does not provide target ${HEADER_CHECK_TARGET}")
endif()

get_target_property(inc_dirs ${HEADER_CHECK_TARGET} INTERFACE_INCLUDE_DIRECTORIES)
if(NOT inc_dirs)
  message(FATAL_ERROR "${HEADER_CHECK_TARGET} has no INTERFACE_INCLUDE_DIRECTORIES")
endif()

# ---------------------------------------------------------------------------
# Collect installed public headers
# ---------------------------------------------------------------------------
set(all_headers)
foreach(dir IN LISTS inc_dirs)
  if(NOT IS_DIRECTORY "${dir}")
    continue()
  endif()
  # Glob all headers, then keep only this library's (the include dir may also
  # contain headers of other libraries installed into the same prefix)
  file(GLOB_RECURSE headers RELATIVE "${dir}" CONFIGURE_DEPENDS
       "${dir}/*.h"
       "${dir}/*.hh"
       "${dir}/*.hpp"
       "${dir}/*.hxx")
  list(FILTER headers INCLUDE REGEX "${HEADER_CHECK_INCLUDE_REGEX}")
  if(HEADER_CHECK_EXCLUDE_REGEX)
    list(FILTER headers EXCLUDE REGEX "${HEADER_CHECK_EXCLUDE_REGEX}")
  endif()
  list(APPEND all_headers ${headers})
endforeach()

list(REMOVE_DUPLICATES all_headers)
list(SORT all_headers)

# TODO: Remove this later again once the workflow works.
# Restrict to the allowlist, if present
if(EXISTS "${HEADER_CHECK_LIST_FILE}")
  # Re-run configure automatically when the list changes
  set_property(DIRECTORY APPEND PROPERTY CMAKE_CONFIGURE_DEPENDS "${HEADER_CHECK_LIST_FILE}")

  file(STRINGS "${HEADER_CHECK_LIST_FILE}" lines)
  set(selected)
  set(missing)
  foreach(line IN LISTS lines)
    string(STRIP "${line}" line)
    if(line STREQUAL "" OR line MATCHES "^#")
      continue()
    endif()
    list(APPEND selected "${line}")
    if(NOT line IN_LIST all_headers)
      list(APPEND missing "${line}")
    endif()
  endforeach()

  if(missing)
    list(JOIN missing "\n  " missing_str)
    message(FATAL_ERROR
      "Headers listed in ${HEADER_CHECK_LIST_FILE} are not installed "
      "(or do not match HEADER_CHECK_INCLUDE_REGEX):\n  ${missing_str}")
  endif()
  if(NOT selected)
    message(FATAL_ERROR "${HEADER_CHECK_LIST_FILE} lists no headers. "
                        "Add some, or delete the file to check all headers.")
  endif()

  list(LENGTH all_headers num_installed)
  list(REMOVE_DUPLICATES selected)
  set(all_headers ${selected})
  list(LENGTH all_headers num_selected)
  message(STATUS "Header check: allowlist active, checking ${num_selected} of ${num_installed} headers")
endif()

list(LENGTH all_headers num_headers)
if(num_headers EQUAL 0)
  message(FATAL_ERROR
    "No headers matching '${HEADER_CHECK_INCLUDE_REGEX}' found in: ${inc_dirs}\n"
    "Either the install is broken or HEADER_CHECK_INCLUDE_REGEX is wrong.")
endif()
message(STATUS "Header check: ${num_headers} public headers")

# ---------------------------------------------------------------------------
# Generate sources
# file(CONFIGURE) only rewrites files whose content changed, so
# re-configuring does not trigger needless rebuilds.
# ---------------------------------------------------------------------------
# Must not be named like the executable, which is written to the same dir
set(gen_dir "${CMAKE_CURRENT_BINARY_DIR}/generated_header_tests")
set(sources)
set(all_in_one_content "// All public headers in one translation unit\n")

foreach(h IN LISTS all_headers)
  string(MAKE_C_IDENTIFIER "${h}" id)
  set(src "${gen_dir}/${id}.cpp")
  file(CONFIGURE OUTPUT "${src}" @ONLY CONTENT
"#include <@h@>
#include <@h@> // second include: include guard / #pragma once check
")
  list(APPEND sources "${src}")
  string(APPEND all_in_one_content "#include <${h}>\n")
endforeach()

file(CONFIGURE OUTPUT "${gen_dir}/all_headers.cpp" @ONLY CONTENT "${all_in_one_content}")
file(CONFIGURE OUTPUT "${gen_dir}/main.cpp"        @ONLY CONTENT "int main() { return 0; }\n")

list(APPEND sources "${gen_dir}/all_headers.cpp" "${gen_dir}/main.cpp")

# ---------------------------------------------------------------------------
# Build target: part of the normal build, so a header problem fails the
# build just as it would in a downstream project.
# ---------------------------------------------------------------------------
add_executable(header_check ${sources})
target_link_libraries(header_check PRIVATE ${HEADER_CHECK_TARGET})

set_target_properties(header_check PROPERTIES
  CXX_STANDARD            ${HEADER_CHECK_CXX_STANDARD}
  CXX_STANDARD_REQUIRED   ON
  CXX_EXTENSIONS          OFF
  UNITY_BUILD             OFF   # unity builds would merge the TUs again
  NO_SYSTEM_FROM_IMPORTED ON)   # don't hide warnings coming from the headers

target_compile_options(header_check PRIVATE
  $<$<CXX_COMPILER_ID:GNU,Clang,AppleClang>:-Wall -Wextra -Wpedantic -Werror>
  $<$<CXX_COMPILER_ID:MSVC>:/W4 /WX /permissive- /Zc:__cplusplus>)

# ---------------------------------------------------------------------------
# CTest entry: the real check already happened when header_check was built.
# Running the (empty) executable only lists the check in the CTest summary.
# ---------------------------------------------------------------------------
add_test(NAME t8code_post_install_header_check COMMAND header_check)
set_tests_properties(t8code_post_install_header_check PROPERTIES LABELS "post_install;headers")