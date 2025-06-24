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
option( T8CODE_BUILD_FORTRAN_INTERFACE "Build t8code's Fortran interface" OFF )

option( T8CODE_ENABLE_MPI "Enable t8code's features which rely on MPI" ON )
option( T8CODE_ENABLE_VTK "Enable t8code's features which rely on VTK" OFF )
option( T8CODE_ENABLE_OCC "Enable t8code's features which rely on OpenCASCADE" OFF )
option( T8CODE_ENABLE_NETCDF "Enable t8code's features which rely on netCDF" OFF )

option( T8CODE_USE_SYSTEM_SC "Use system-installed sc library" OFF )
option( T8CODE_USE_SYSTEM_P4EST "Use system-installed p4est library" OFF )

option( T8CODE_BUILD_DOCUMENTATION "Build t8code's documentation" OFF )
cmake_dependent_option( T8CODE_BUILD_DOCUMENTATION_SPHINX "Build t8code's documentation using sphinx" OFF "T8CODE_BUILD_DOCUMENTATION" OFF )

set(T8CODE_TEST_LEVEL 0 CACHE STRING "Test level: 0 for full tests, 1 for less thorough tests, 2 for minimal tests. (WARNING: Use with care.)")
set(T8CODE_CUSTOM_PARALLEL_TEST_COMMAND "" CACHE STRING "Define a custom command for parallel tests , e.g.: mpirun -np 8 (overwrites standard mpirun -np 4 if build with mpi)")
set(T8CODE_CUSTOM_SERIAL_TEST_COMMAND "" CACHE STRING "Define a custom command for serial tests.")

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

set(CMAKE_MODULE_PATH "${PROJECT_SOURCE_DIR}/cmake" ${CMAKE_MODULE_PATH})

if( T8CODE_BUILD_FORTRAN_INTERFACE )
    enable_language( Fortran )
endif()

if( T8CODE_ENABLE_MPI )
    if( T8CODE_BUILD_FORTRAN_INTERFACE )
      find_package( MPI COMPONENTS C Fortran REQUIRED )
    else()
      find_package( MPI COMPONENTS C REQUIRED )
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
    find_package( netCDF REQUIRED )
    if(netCDF_FOUND)
        message("Found netCDF")
        include(cmake/CheckNetCDFPar.cmake)
    endif (netCDF_FOUND)
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
    find_package( SC REQUIRED PATHS /path/to/system/sc )
else()
    include(FetchContent)
    FetchContent_Declare(
        sc
        GIT_REPOSITORY https://github.com/cburstedde/libsc.git
        GIT_TAG master
    )
    set( SC_BUILD_EXAMPLES ${T8CODE_BUILD_TPL_EXAMPLES} )
    set( SC_BUILD_TESTING ${T8CODE_BUILD_TPL_TESTS} )
    
    # Capture the list of variables before adding the subdirectory to mark the added ones as advanced.
    get_cmake_property(_vars_before_sc VARIABLES )
    add_subdirectory( ${CMAKE_CURRENT_LIST_DIR}/sc )
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
endif()

if (T8CODE_USE_SYSTEM_P4EST)
    find_package(P4EST REQUIRED PATHS /path/to/system/p4est)
else()
    include(FetchContent)
    FetchContent_Declare(
        p4est
        GIT_REPOSITORY https://github.com/cburstedde/p4est.git
        GIT_TAG master
    )
    FetchContent_GetProperties(p4est)

    if(NOT p4est_POPULATED)
        FetchContent_Populate(p4est)
        add_subdirectory(${p4est_SOURCE_DIR} ${p4est_BINARY_DIR})
    endif()

    set(P4EST_BUILD_EXAMPLES ${T8CODE_BUILD_TPL_EXAMPLES})
    set(P4EST_BUILD_TESTING ${T8CODE_BUILD_TPL_TESTS})

    # Mark the new variables as advanced
    mark_as_advanced(FORCE P4EST_BUILD_EXAMPLES P4EST_BUILD_TESTING)

endif()

add_subdirectory( ${CMAKE_CURRENT_LIST_DIR}/src )

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
