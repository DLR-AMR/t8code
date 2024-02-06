cmake_minimum_required( VERSION 3.16 )
project( T8CODE DESCRIPTION "Parallel algorithms and data structures for tree-based AMR with arbitrary element shapes." LANGUAGES C CXX )
include( CTest )


option( T8CODE_ENABLE_MPI "" ON )
option( T8CODE_ENABLE_DYNAMIC_LINKAGE "" ON )

option( T8CODE_BUILD_TESTS "" ON )
option( T8CODE_BUILD_TUTORIALS "" ON )
option( T8CODE_BUILD_EXAMPLES "" ON )


if( NOT CMAKE_BUILD_TYPE )
    set( CMAKE_BUILD_TYPE "Release" )
endif()

if( NOT CMAKE_INSTALL_PREFIX )
    set( CMAKE_INSTALL_PREFIX ${CMAKE_CURRENT_LIST_DIR}/build )
endif()

if( NOT CMAKE_CXX_STANDARD )
    set( CMAKE_CXX_STANDARD 17 )
endif()


if( T8CODE_ENABLE_MPI )
    find_package( MPI COMPONENTS C REQUIRED )
    if( NOT MPIEXEC_EXECUTABLE )
        message( FATAL_ERROR "MPIEXEC is missing, MPI will not work properly" )
    endif()
    set( SC_ENABLE_MPI ON )
    set( mpi ON ) # This is very dirty, we should consider fixing this in the p4est repo
endif()

add_subdirectory( ${CMAKE_CURRENT_LIST_DIR}/sc EXCLUDE_FROM_ALL )
add_subdirectory( ${CMAKE_CURRENT_LIST_DIR}/p4est EXCLUDE_FROM_ALL )
add_subdirectory( ${CMAKE_CURRENT_LIST_DIR}/src )

install( TARGETS T8 DESTINATION ${CMAKE_INSTALL_PREFIX}/lib )

if ( T8CODE_BUILD_TESTS )
    add_subdirectory( ${CMAKE_CURRENT_LIST_DIR}/test )
endif()

if ( T8CODE_BUILD_TUTORIALS )
    add_subdirectory( ${CMAKE_CURRENT_LIST_DIR}/tutorials )
endif()

if ( T8CODE_BUILD_EXAMPLES )
    add_subdirectory( ${CMAKE_CURRENT_LIST_DIR}/example )
endif()