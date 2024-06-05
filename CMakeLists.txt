cmake_minimum_required( VERSION 3.16 )
project( T8CODE DESCRIPTION "Parallel algorithms and data structures for tree-based AMR with arbitrary element shapes." LANGUAGES C CXX )
include( CTest )

option( T8CODE_BUILD_AS_SHARED_LIBRARY "Whether t8code should be built as a shared or a static library" ON )
option( T8CODE_BUILD_TESTS "Build t8code's automated tests" ON )
option( T8CODE_BUILD_TUTORIALS "Build t8code's tutorials" ON )
option( T8CODE_BUILD_EXAMPLES "Build t8code's examples" ON )

option( T8CODE_ENABLE_MPI "Enable t8code's features which rely on MPI" ON )
option( T8CODE_ENABLE_VTK "Enable t8code's features which rely on VTK" OFF )

option( T8CODE_USE_SYSTEM_SC "Use system-installed sc library" OFF )
option( T8CODE_USE_SYSTEM_P4EST "Use system-installed p4est library" OFF )

option( BUILD_DOCUMENTATION "Build t8code's documentation" ON)

if( NOT CMAKE_BUILD_TYPE )
    set( CMAKE_BUILD_TYPE "Release" )
endif()

set_property(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS Release Debug)

set( CMAKE_C_STANDARD 11 )
set( CMAKE_C_STANDARD_REQUIRED ON )
set( CMAKE_C_EXTENSIONS OFF )
list( APPEND CMAKE_C_FLAGS "-Wall" )

set( CMAKE_CXX_STANDARD 17 )
set( CMAKE_CXX_STANDARD_REQUIRED ON )
set( CMAKE_CXX_EXTENSIONS OFF )
list( APPEND CMAKE_CXX_FLAGS "-Wall" )

if( T8CODE_ENABLE_MPI )
    find_package( MPI COMPONENTS C REQUIRED )
    if( NOT MPIEXEC_EXECUTABLE )
        message( FATAL_ERROR "MPIEXEC was not found" )
    endif()
    set( mpi ON ) # This is very dirty and will be fixed in the libsc repo (https://github.com/cburstedde/libsc/pull/178)
    # set( SC_ENABLE_MPI ON ) # When the fix gets merged, replace the previous line with this one
endif()

if( T8CODE_ENABLE_VTK )
    find_package( VTK REQUIRED )
endif()

# Override default for this libsc option
set( BUILD_SHARED_LIBS ON CACHE BOOL "Build libsc as a shared library" )

# Prevent `libsc` and `p4est` from overwriting the default install prefix.
set(CMAKE_INSTALL_PREFIX_INITIALIZED_TO_DEFAULT FALSE)

# Rpath options necessary for shared library install to work correctly in user projects.
set(CMAKE_INSTALL_NAME_DIR ${CMAKE_INSTALL_PREFIX}/lib)
set(CMAKE_INSTALL_RPATH ${CMAKE_INSTALL_PREFIX}/lib)
set(CMAKE_INSTALL_RPATH_USE_LINK_PATH true)

if ( T8CODE_USE_SYSTEM_SC )
    find_package( SC REQUIRED PATHS /path/to/system/sc )
else()
    add_subdirectory( ${CMAKE_CURRENT_LIST_DIR}/sc )
endif()

if ( T8CODE_USE_SYSTEM_P4EST )
    find_package( P4EST REQUIRED PATHS /path/to/system/p4est )
else()
    add_subdirectory( ${CMAKE_CURRENT_LIST_DIR}/p4est )
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

if ( BUILD_DOCUMENTATION ) 
    add_subdirectory( ${CMAKE_CURRENT_LIST_DIR}/doc )
endif()
