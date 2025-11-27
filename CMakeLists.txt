#  This file is part of t8code.
#  t8code is a C library to manage a collection (a forest) of multiple
#  connected adaptive space-trees of general element types in parallel.
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

cmake_minimum_required( VERSION 3.16 )

include( cmake/GitProjectVersion.cmake )

project(
    T8CODE
    DESCRIPTION "Parallel algorithms and data structures for tree-based AMR with arbitrary element shapes."
    LANGUAGES C CXX
    VERSION "${T8CODE_VERSION_MAJOR}.${T8CODE_VERSION_MINOR}.${T8CODE_VERSION_PATCH}"
    )

include( GNUInstallDirs)
include( CTest )
include( CMakeDependentOption )

option( T8CODE_BUILD_AS_SHARED_LIBRARY "Whether t8code should be built as a shared or a static library" ON )
option( T8CODE_BUILD_PEDANTIC "Compile t8code with `-pedantic` as done in the Github CI." OFF )
option( T8CODE_BUILD_WALL "Compile t8code with `-Wall` as done in the Github CI." OFF )
option( T8CODE_BUILD_WERROR "Compile t8code with `-Werror` as done in the Github CI." OFF )
option( T8CODE_BUILD_WEXTRA "Compile t8code with extra warnings as done in the Github CI." OFF )
option( T8CODE_EXPORT_COMPILE_COMMANDS "Export the compile commands as json. Can be used by IDEs for code completion (e.g. intellisense, clangd)" OFF )

option( T8CODE_BUILD_TESTS "Build t8code's automated tests" ON )
cmake_dependent_option( T8CODE_BUILD_TPL_TESTS "Build the tests from libsc and p4est" OFF "T8CODE_BUILD_TESTS" OFF )

option( T8CODE_BUILD_EXAMPLES "Build t8code's examples" ON )
cmake_dependent_option( T8CODE_BUILD_TPL_EXAMPLES "Build the examples from libsc and p4est" OFF "T8CODE_BUILD_EXAMPLES" OFF )

option( T8CODE_BUILD_TUTORIALS "Build t8code's tutorials" ON )
option( T8CODE_BUILD_BENCHMARKS "Build t8code's benchmarks" ON )
option( T8CODE_BUILD_MESH_HANDLE "Build t8code's mesh handle" OFF )
option( T8CODE_BUILD_FORTRAN_INTERFACE "Build t8code's Fortran interface" OFF )

option( T8CODE_ENABLE_MPI "Enable t8code's features which rely on MPI" ON )
option( T8CODE_ENABLE_VTK "Enable t8code's features which rely on VTK" OFF )
option( T8CODE_ENABLE_OCC "Enable t8code's features which rely on OpenCASCADE" OFF )
option( T8CODE_ENABLE_NETCDF "Enable t8code's features which rely on netCDF" OFF )

option( T8CODE_USE_SYSTEM_SC "Use system-installed sc library" OFF )
option( T8CODE_USE_SYSTEM_P4EST "Use system-installed p4est library" OFF )

option( T8CODE_BUILD_DOCUMENTATION "Build t8code's documentation" OFF )
cmake_dependent_option( T8CODE_BUILD_DOCUMENTATION_SPHINX "Build t8code's documentation using sphinx" OFF "T8CODE_BUILD_DOCUMENTATION" OFF )

set(T8CODE_CUSTOM_PARALLEL_TEST_COMMAND "" CACHE STRING "Define a custom command for parallel tests , e.g.: mpirun -np 8 (overwrites standard mpirun -np 4 if build with mpi)")
set(T8CODE_CUSTOM_SERIAL_TEST_COMMAND "" CACHE STRING "Define a custom command for serial tests.")

# Reading the option T8CODE_TEST_LEVEL: Possible choices are T8_TEST_LEVEL_FULL, T8_TEST_LEVEL_MEDIUM, or T8_TEST_LEVEL_BASIC. Default is T8_TEST_LEVEL_FULL.
set(T8CODE_TEST_LEVEL "T8_TEST_LEVEL_FULL" CACHE STRING "Test level: T8_TEST_LEVEL_FULL for full tests, T8_TEST_LEVEL_MEDIUM for less thorough tests, T8_TEST_LEVEL_BASIC for minimal tests. (WARNING: Use with care.)")
set_property(CACHE T8CODE_TEST_LEVEL PROPERTY STRINGS "T8_TEST_LEVEL_FULL" "T8_TEST_LEVEL_MEDIUM" "T8_TEST_LEVEL_BASIC")

# Parse string T8CODE_TEST_LEVEL into integer T8_TEST_LEVEL_INT for easier usage within the code.
if(T8CODE_TEST_LEVEL STREQUAL "T8_TEST_LEVEL_BASIC" )
    set(T8_TEST_LEVEL_INT 2)
elseif(T8CODE_TEST_LEVEL STREQUAL "T8_TEST_LEVEL_MEDIUM" )
    set(T8_TEST_LEVEL_INT 1)
elseif(T8CODE_TEST_LEVEL STREQUAL "T8_TEST_LEVEL_FULL" )
    set(T8_TEST_LEVEL_INT 0)
else()
    message( FATAL_ERROR "Invalid string for T8CODE_TEST_LEVEL: ${T8CODE_TEST_LEVEL}. Valid options are T8_TEST_LEVEL_FULL, T8_TEST_LEVEL_MEDIUM, or T8_TEST_LEVEL_BASIC.")
endif()

# Set a default build type if none was specified
set(default_build_type "Release")

if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
  message(STATUS "Setting build type to '${default_build_type}' as none was specified.")
  set(CMAKE_BUILD_TYPE "${default_build_type}" CACHE
      STRING "Choose the type of build. Build types available: Release Debug RelWithDebInfo" FORCE)
  # Set the possible values of build type for cmake-gui
  set_property(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS
    "Debug" "Release" "RelWithDebInfo")
endif()

if( CMAKE_BUILD_TYPE STREQUAL "Debug" )
  # Option to generate code coverage target.
  option( T8CODE_CODE_COVERAGE "Enable code coverage reporting" OFF)
  if(T8CODE_CODE_COVERAGE)
    include(./cmake/CodeCoverage.cmake)
  endif()
endif()

set(CMAKE_MODULE_PATH "${PROJECT_SOURCE_DIR}/cmake" ${CMAKE_MODULE_PATH})

if( T8CODE_BUILD_FORTRAN_INTERFACE )
    enable_language( Fortran )
endif()

if( T8CODE_ENABLE_MPI )
    if( T8CODE_BUILD_FORTRAN_INTERFACE )
      find_package( MPI 3.0 COMPONENTS C Fortran REQUIRED )
    else()
      find_package( MPI 3.0 COMPONENTS C REQUIRED )
    endif()

    if( NOT MPIEXEC_EXECUTABLE )
        message( FATAL_ERROR "MPIEXEC was not found" )
    endif()
    set( SC_ENABLE_MPI ON )
endif()

if( T8CODE_ENABLE_VTK )
    find_package( VTK REQUIRED COMPONENTS
        IOXML CommonExecutionModel CommonDataModel
        IOGeometry IOXMLParser IOParallelXML IOPLY
        ParallelMPI FiltersCore vtksys CommonCore zlib IOLegacy)
    if(VTK_FOUND)
        message("Found VTK")
    endif (VTK_FOUND)
endif( T8CODE_ENABLE_VTK )

if( T8CODE_ENABLE_OCC )
    find_package( OpenCASCADE REQUIRED COMPONENTS
    TKBO TKPrim TKTopAlgo
    TKGeomAlgo TKBRep
    TKG3d TKG2d TKMath TKernel )
    if(OpenCASCADE_FOUND)
        message("Found OpenCASCADE")
    endif (OpenCASCADE_FOUND)
endif( T8CODE_ENABLE_OCC )

if( T8CODE_ENABLE_NETCDF )
    find_package( NetCDF REQUIRED )
    if(NetCDF_FOUND)
        message("Found NetCDF")
        include(cmake/CheckNetCDFPar.cmake)
    endif (NetCDF_FOUND)
endif( T8CODE_ENABLE_NETCDF )

# Override default for this libsc option
set( BUILD_SHARED_LIBS ON CACHE BOOL "Build libsc as a shared library" )

# Prevent `libsc` and `p4est` from overwriting the default install prefix.
set(CMAKE_INSTALL_PREFIX_INITIALIZED_TO_DEFAULT FALSE)

# Rpath options necessary for shared library install to work correctly in user projects.
# Compute the relative Rpath for proper library installation
set(libPath ${CMAKE_CURRENT_BINARY_DIR}/${CMAKE_INSTALL_LIBDIR}) # Library directory)
cmake_path(RELATIVE_PATH libPath
    BASE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/${CMAKE_INSTALL_BINDIR} # Binary directory
    OUTPUT_VARIABLE relativeRpath
)
set(CMAKE_INSTALL_RPATH $ORIGIN $ORIGIN/${relativeRpath})
set(CMAKE_INSTALL_RPATH_USE_LINK_PATH true)

if ( T8CODE_USE_SYSTEM_SC )
    find_package( SC REQUIRED )
else()
    set( SC_BUILD_EXAMPLES ${T8CODE_BUILD_TPL_EXAMPLES} )
    set( SC_BUILD_TESTING ${T8CODE_BUILD_TPL_TESTS} )

    # Capture the list of variables before adding the subdirectory to mark the added ones as advanced.
    get_cmake_property(_vars_before_sc VARIABLES )
    list( APPEND CMAKE_MESSAGE_CONTEXT "sc" )
    add_subdirectory( ${CMAKE_CURRENT_LIST_DIR}/sc )
    list( POP_BACK CMAKE_MESSAGE_CONTEXT )
    add_library(SC INTERFACE)
    # Capture the list of variables after adding the subdirectory
    get_cmake_property(_vars_after_sc VARIABLES )
    # Compute the difference (new variables added by the subdirectory)
    set( _new_sc_vars)
    foreach( _var IN LISTS _vars_after_sc )
        if( NOT _var IN_LIST _vars_before_sc )
            list( APPEND _new_sc_vars ${_var} )
        endif()
    endforeach()

    # Mark the new variables as advanced
    mark_as_advanced( FORCE ${_new_sc_vars} )

    ########################## Fix for OpenMPI ############################################
    # OpenMPI sometimes has a mismatch between their C and CXX interface which produces a warning in SC.
    # The warning makes t8codes build fail when Werror is activated. This fix checks if t8code
    # is linked against OpenMPI and then suppresses that Werror treats this warning as error.
    # The warning is still displayed while compiling.
    # We compile a test program mich checks if a OpenMPI specific macro is defined and then apply the fix.

    include(CheckIncludeFile)

    # Check if mpi.h is available
    check_include_file("mpi.h" MPI_HEADER_FOUND)

    if(MPI_HEADER_FOUND)
        # Create a temporary C++ source file to check for OpenMPI-specific macros
        include(CheckCXXSourceCompiles)
        check_cxx_source_compiles("
        #include <mpi.h>
        #ifndef OPEN_MPI
        #error 'Not OpenMPI'
        #endif
        int main() { return 0; }" HAS_OPENMPI)

        if(HAS_OPENMPI)
            message(STATUS "OpenMPI detected via mpi.h")
            # Suppress warning
            target_compile_options(SC::SC INTERFACE -Wno-error=cast-function-type)
        endif()
    else()
        message(WARNING "mpi.h not found!")
    endif()
    ########################## End of fix for OpenMPI ############################################
endif()

if ( T8CODE_USE_SYSTEM_P4EST )
    find_package( P4EST REQUIRED )
else()
    set( P4EST_BUILD_EXAMPLES ${T8CODE_BUILD_TPL_EXAMPLES} )
    set( P4EST_BUILD_TESTING ${T8CODE_BUILD_TPL_TESTS} )

    # Capture the list of variables before adding the subdirectory to mark the added ones as advanced.
    get_cmake_property( _vars_before_p4est VARIABLES )
    list( APPEND CMAKE_MESSAGE_CONTEXT "p4est" )
    add_subdirectory( ${CMAKE_CURRENT_LIST_DIR}/p4est )
    list( POP_BACK CMAKE_MESSAGE_CONTEXT )
    # Capture the list of variables after adding the subdirectory
    get_cmake_property( _vars_after_p4est VARIABLES )
    # Compute the difference (new variables added by the subdirectory)
    set( _new_p4est_vars )
    foreach( _var IN LISTS _vars_after_p4est )
        if( NOT _var IN_LIST _vars_before_p4est )
            list( APPEND _new_p4est_vars ${_var} )
        endif()
    endforeach()

    # Mark the new variables as advanced
    mark_as_advanced( FORCE ${_new_p4est_vars} )
endif()

add_subdirectory( ${CMAKE_CURRENT_LIST_DIR}/src )

if( T8CODE_BUILD_MESH_HANDLE )
  add_subdirectory( ${CMAKE_CURRENT_LIST_DIR}/mesh_handle )
endif()

if ( T8CODE_BUILD_TESTS )
    add_subdirectory( ${CMAKE_CURRENT_LIST_DIR}/test )
endif()

if ( T8CODE_BUILD_TUTORIALS )
    add_subdirectory( ${CMAKE_CURRENT_LIST_DIR}/tutorials )
endif()

if ( T8CODE_BUILD_EXAMPLES )
    add_subdirectory( ${CMAKE_CURRENT_LIST_DIR}/example )
endif()

if ( T8CODE_BUILD_BENCHMARKS )
    add_subdirectory( ${CMAKE_CURRENT_LIST_DIR}/benchmarks )
endif()

if ( T8CODE_BUILD_DOCUMENTATION )
    add_subdirectory( ${CMAKE_CURRENT_LIST_DIR}/doc )
endif()

if( T8CODE_BUILD_FORTRAN_INTERFACE )
    add_subdirectory( ${CMAKE_CURRENT_LIST_DIR}/api/t8_fortran_interface )

    if( NOT T8CODE_ENABLE_MPI )
        message( FATAL_ERROR "Fortran API only available when MPI is enabled." )
    endif()
endif()

include (cmake/CPackConfig.cmake)
