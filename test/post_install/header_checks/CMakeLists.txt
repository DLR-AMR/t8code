# Post-install header check for an installed CMake package.
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
# Usage (against an install prefix, ideally on a machine/container without
# the source and build trees):
#   cmake -S <this dir> -B build-header-check -DCMAKE_PREFIX_PATH=<install prefix>
#   cmake --build build-header-check
#   ctest --test-dir build-header-check     # optional, see add_test below

cmake_minimum_required(VERSION 3.20)
project(MyLibHeaderCheck LANGUAGES C CXX)

# ---------------------------------------------------------------------------
# Configuration — adapt the defaults to your library
# ---------------------------------------------------------------------------
set(HEADER_CHECK_PACKAGE       "T8CODE"        CACHE STRING "Package name passed to find_package()")
set(HEADER_CHECK_TARGET        "T8CODE::T8" CACHE STRING "Imported target whose headers are checked")
set(HEADER_CHECK_SUBDIR        ""        CACHE STRING "Header subdirectory below the include dir")
set(HEADER_CHECK_EXCLUDE_REGEX "/detail/"     CACHE STRING "Regex of headers to skip (empty = check all)")
set(HEADER_CHECK_CXX_STANDARD  "20"           CACHE STRING "Minimum C++ standard the headers must support")

# ---------------------------------------------------------------------------
# Find the installed package
# ---------------------------------------------------------------------------
find_package(${HEADER_CHECK_PACKAGE} REQUIRED)

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
  if(NOT IS_DIRECTORY "${dir}/${HEADER_CHECK_SUBDIR}")
    continue()
  endif()
  file(GLOB_RECURSE headers RELATIVE "${dir}" CONFIGURE_DEPENDS
       "${dir}/${HEADER_CHECK_SUBDIR}/*.h"
       "${dir}/${HEADER_CHECK_SUBDIR}/*.hh"
       "${dir}/${HEADER_CHECK_SUBDIR}/*.hpp"
       "${dir}/${HEADER_CHECK_SUBDIR}/*.hxx")
  if(HEADER_CHECK_EXCLUDE_REGEX)
    list(FILTER headers EXCLUDE REGEX "${HEADER_CHECK_EXCLUDE_REGEX}")
  endif()
  list(APPEND all_headers ${headers})
endforeach()

list(REMOVE_DUPLICATES all_headers)
list(SORT all_headers)

list(LENGTH all_headers num_headers)
if(num_headers EQUAL 0)
  message(FATAL_ERROR
    "No headers found below '${HEADER_CHECK_SUBDIR}' in: ${inc_dirs}\n"
    "Either the install is broken or HEADER_CHECK_SUBDIR is wrong.")
endif()
message(STATUS "Header check: ${num_headers} public headers")

# ---------------------------------------------------------------------------
# Generate sources
# file(CONFIGURE) only rewrites files whose content changed, so
# re-configuring does not trigger needless rebuilds.
# ---------------------------------------------------------------------------
set(gen_dir "${CMAKE_CURRENT_BINARY_DIR}/header_check")
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
# Build target
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
# CTest integration: the test is that the target builds (compile + link).
# The executable itself does nothing, so it is not run.
# ---------------------------------------------------------------------------
enable_testing()
add_test(NAME header_check
         COMMAND ${CMAKE_COMMAND} --build ${CMAKE_BINARY_DIR}
                 --target header_check --config $<CONFIG>)