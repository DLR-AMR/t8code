cmake_minimum_required( VERSION 3.16 )
project( T8CODE DESCRIPTION "Parallel algorithms and data structures for tree-based AMR with arbitrary element shapes." LANGUAGES C CXX )
include( CTest )


option( T8CODE_BUILD_AS_SHARED_LIBRARY "Whether t8code should be built as a shared or a static library" ON )
option( T8CODE_BUILD_TESTS "Build t8code's automated tests" ON )
option( T8CODE_BUILD_TUTORIALS "Build t8code's tutorials" ON )
option( T8CODE_BUILD_EXAMPLES "Build t8code's examples" ON )

option( T8CODE_ENABLE_MPI "Enable t8code's features which rely on MPI" ON )
option( T8CODE_ENABLE_VTK "Enable t8code's features which rely on VTK" OFF )


if( NOT CMAKE_BUILD_TYPE )
    set( CMAKE_BUILD_TYPE "Release" )
endif()

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
    set( SC_ENABLE_MPI ON )
    set( mpi ON ) # This is very dirty, we should consider fixing this in the p4est repo
endif()

if( T8CODE_ENABLE_VTK )
    find_package( VTK REQUIRED )
endif()

add_subdirectory( ${CMAKE_CURRENT_LIST_DIR}/sc )
add_subdirectory( ${CMAKE_CURRENT_LIST_DIR}/p4est )
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