add_library( gtest ${CMAKE_CURRENT_LIST_DIR}/../thirdparty/googletest-mpi/gtest/gtest-all.cc )
target_include_directories( gtest PUBLIC ${CMAKE_CURRENT_LIST_DIR}/../thirdparty/googletest-mpi ${CMAKE_CURRENT_LIST_DIR}/.. )

add_executable( t8_gtest_cmesh_bcast t8_gtest_main.cxx t8_cmesh/t8_gtest_bcast.cxx )
target_link_libraries( t8_gtest_cmesh_bcast PRIVATE T8 gtest )
add_test( NAME t8_gtest_cmesh_bcast COMMAND ./t8_gtest_cmesh_bcast )