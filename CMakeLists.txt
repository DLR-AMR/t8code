cmake_minimum_required( VERSION 3.16 )
project( T8CODE DESCRIPTION "Parallel algorithms and data structures for tree-based AMR with arbitrary element shapes." LANGUAGES C CXX )
include( CTest )

option( T8CODE_BUILD_TESTS "" ON )
option( T8CODE_BUILD_TUTORIALS "" ON )
# option( T8CODE_ENABLE_MPI "" OFF )

if( NOT CMAKE_CXX_STANDARD )
    set( CMAKE_CXX_STANDARD 17 )
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