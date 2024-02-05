add_executable( t8_step0 ${CMAKE_CURRENT_LIST_DIR}/general/t8_step0_helloworld.cxx )
target_link_libraries( t8_step0 PRIVATE T8 )

add_executable( t8_step1 ${CMAKE_CURRENT_LIST_DIR}/general/t8_step1_coarsemesh.cxx )
target_link_libraries( t8_step1 PRIVATE T8 )

add_executable( t8_step2 ${CMAKE_CURRENT_LIST_DIR}/general/t8_step2_uniform_forest.cxx )
target_link_libraries( t8_step2 PRIVATE T8 )

add_executable( t8_step3 ${CMAKE_CURRENT_LIST_DIR}/general/t8_step3_adapt_forest.cxx
                         ${CMAKE_CURRENT_LIST_DIR}/general/t8_step3_main.cxx )
target_include_directories( t8_step3 PRIVATE ${CMAKE_CURRENT_LIST_DIR}/.. )
target_link_libraries( t8_step3 PRIVATE T8 )

add_executable( t8_step4 ${CMAKE_CURRENT_LIST_DIR}/general/t8_step3_adapt_forest.cxx
                         ${CMAKE_CURRENT_LIST_DIR}/general/t8_step4_partition_balance_ghost.cxx
                         ${CMAKE_CURRENT_LIST_DIR}/general/t8_step4_main.cxx )
target_include_directories( t8_step4 PRIVATE ${CMAKE_CURRENT_LIST_DIR}/.. )
target_link_libraries( t8_step4 PRIVATE T8 )

add_executable( t8_step5 ${CMAKE_CURRENT_LIST_DIR}/general/t8_step3_adapt_forest.cxx
                        #  ${CMAKE_CURRENT_LIST_DIR}/general/t8_step5_element_data_c_interface.c
                         ${CMAKE_CURRENT_LIST_DIR}/general/t8_step5_element_data.cxx
                         ${CMAKE_CURRENT_LIST_DIR}/general/t8_step5_main.cxx )
target_include_directories( t8_step5 PRIVATE ${CMAKE_CURRENT_LIST_DIR}/.. )
target_link_libraries( t8_step5 PRIVATE T8 )

add_executable( t8_step6 ${CMAKE_CURRENT_LIST_DIR}/general/t8_step3_adapt_forest.cxx
                         ${CMAKE_CURRENT_LIST_DIR}/general/t8_step6_stencil.cxx
                         ${CMAKE_CURRENT_LIST_DIR}/general/t8_step6_main.cxx )
target_include_directories( t8_step6 PRIVATE ${CMAKE_CURRENT_LIST_DIR}/.. )
target_link_libraries( t8_step6 PRIVATE T8 )

add_executable( t8_step7 ${CMAKE_CURRENT_LIST_DIR}/general/t8_step7_interpolation.cxx
                         ${CMAKE_CURRENT_LIST_DIR}/general/t8_step7_main.cxx )
target_include_directories( t8_step7 PRIVATE ${CMAKE_CURRENT_LIST_DIR}/.. )
target_link_libraries( t8_step7 PRIVATE T8 )