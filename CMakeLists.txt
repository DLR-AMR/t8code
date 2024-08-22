cmake_minimum_required( VERSION 3.16 )

include(cmake/GitProjectVersion.cmake)

project( T8CODE
  DESCRIPTION "Parallel algorithms and data structures for tree-based AMR with arbitrary element shapes."
  LANGUAGES C CXX
  VERSION "${T8CODE_VERSION_MAJOR}.${T8CODE_VERSION_MINOR}.${T8CODE_VERSION_PATCH}" )
include( CTest )

option( T8CODE_BUILD_AS_SHARED_LIBRARY "Whether t8code should be built as a shared or a static library" ON )
option( T8CODE_BUILD_PEDANTIC "Compile t8code with `-Wall -pedantic -Werror` as done in the Github CI." OFF )
option( T8CODE_BUILD_TESTS "Build t8code's automated tests" ON )
option( T8CODE_BUILD_TUTORIALS "Build t8code's tutorials" ON )
option( T8CODE_BUILD_EXAMPLES "Build t8code's examples" ON )
option( T8CODE_BUILD_BENCHMARKS "Build t8code's benchmarks" ON )
option( T8CODE_ENABLE_LESS_TESTS "Tests not as thoroughly to speed up the test suite. Tests the same functionality. (WARNING: Use with care.)" OFF )

option( T8CODE_ENABLE_MPI "Enable t8code's features which rely on MPI" ON )
option( T8CODE_ENABLE_VTK "Enable t8code's features which rely on VTK" OFF )
option( T8CODE_ENABLE_OCC "Enable t8code's features which rely on OpenCASCADE" OFF )
option( T8CODE_ENABLE_NETCDF "Enable t8code's features which rely on netCDF" OFF )

option( T8CODE_USE_SYSTEM_SC "Use system-installed sc library" OFF )
option( T8CODE_USE_SYSTEM_P4EST "Use system-installed p4est library" OFF )

option( T8CODE_BUILD_DOCUMENTATION "Build t8code's documentation" OFF )

include(CMakeDependentOption)
cmake_dependent_option( T8CODE_BUILD_DOCUMENTATION_SPHINX "Build t8code's documentation using sphinx" OFF "T8CODE_BUILD_DOCUMENTATION" OFF )

set(T8CODE_CUSTOM_TEST_COMMAND "" CACHE STRING "Define custom test command, e.g.: mpirun -np 8 (overwrites standard mpirun -np 4 if build with mpi)")

if( NOT CMAKE_BUILD_TYPE )
    set( CMAKE_BUILD_TYPE "Release" )
endif()

set_property(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS Release RelWithDebInfo Debug)

set( CMAKE_C_STANDARD 11 )
set( CMAKE_C_STANDARD_REQUIRED ON )
set( CMAKE_C_EXTENSIONS OFF )

set( CMAKE_CXX_STANDARD 17 )
set( CMAKE_CXX_STANDARD_REQUIRED ON )
set( CMAKE_CXX_EXTENSIONS OFF )

set(CMAKE_MODULE_PATH "${PROJECT_SOURCE_DIR}/cmake" ${CMAKE_MODULE_PATH})

if( T8CODE_ENABLE_MPI )
    find_package(MPI COMPONENTS C REQUIRED)
    if( NOT MPIEXEC_EXECUTABLE )
        message( FATAL_ERROR "MPIEXEC was not found" )
    endif()
    set( mpi ON ) # This is very dirty and will be fixed in the libsc repo (https://github.com/cburstedde/libsc/pull/178)
    # set( SC_ENABLE_MPI ON ) # When the fix gets merged, replace the previous line with this one
endif()

if ( T8CODE_USE_SYSTEM_SC )
    find_package( SC REQUIRED )
else( T8CODE_USE_SYSTEM_SC )
    set (SC_LIBRARIES "SC::SC")
    add_subdirectory( ${CMAKE_CURRENT_LIST_DIR}/sc )
endif( T8CODE_USE_SYSTEM_SC )

if ( T8CODE_USE_SYSTEM_P4EST )
    find_package( P4EST REQUIRED )
else( T8CODE_USE_SYSTEM_P4EST )
    set (P4EST_LIBRARIES "P4EST::P4EST")
    add_subdirectory( ${CMAKE_CURRENT_LIST_DIR}/p4est )
endif( T8CODE_USE_SYSTEM_P4EST )

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
set(CMAKE_INSTALL_NAME_DIR ${CMAKE_INSTALL_PREFIX}/lib)
set(CMAKE_INSTALL_RPATH ${CMAKE_INSTALL_PREFIX}/lib)
set(CMAKE_INSTALL_RPATH_USE_LINK_PATH true)

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
