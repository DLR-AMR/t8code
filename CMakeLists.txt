cmake_minimum_required( VERSION 3.16 )
project( T8CODE DESCRIPTION "Parallel algorithms and data structures for tree-based AMR with arbitrary element shapes." LANGUAGES C CXX )
include( CTest )

option( T8CODE_ENABLE_MPI "" ON )
option( T8CODE_BUILD_TESTS "" ON )
option( T8CODE_BUILD_TUTORIALS "" ON )

if( NOT CMAKE_CXX_STANDARD )
    set( CMAKE_CXX_STANDARD 17 )
endif()

if( NOT CMAKE_INSTALL_PREFIX )
    set( CMAKE_INSTALL_PREFIX ${CMAKE_CURRENT_LIST_DIR}/build )
endif()

if( T8CODE_ENABLE_MPI )
    find_package( MPI COMPONENTS C REQUIRED )
    if( NOT MPIEXEC_EXECUTABLE )
        message( FATAL_ERROR "MPIEXEC is missing, MPI will not work properly" )
    endif()
    set( SC_ENABLE_MPI ON )
    set( mpi ON )
endif()

add_subdirectory( sc )
add_subdirectory( p4est )
add_subdirectory( src )

if ( T8CODE_BUILD_TESTS )
    add_subdirectory( test )
endif()

if ( T8CODE_BUILD_TUTORIALS )
    add_subdirectory( tutorials )
endif()

install( TARGETS T8 DESTINATION ${CMAKE_INSTALL_PREFIX}/lib )