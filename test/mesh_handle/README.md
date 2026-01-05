# mesh_handle testing concept #

In this folder, we test the mesh handle functionality. 
The tests are structured as follow: 

- [t8_gtest_mesh_handle.cxx](t8_gtest_mesh_handle.cxx) tests some basic functionality, e.g., if the mesh handle generally works with more than one competence or the iterator. The handle's functionality is not tested in detail, but instead using some exemplary functions.
- [t8_gtest_custom_competence.cxx](t8_gtest_custom_competence.cxx) checks that user defined competences can be used for the mesh handle.
- [t8_gtest_compare_handle_to_forest.cxx](t8_gtest_compare_handle_to_forest.cxx) tests that the functionality of the handle gives the same results as if worked with the forest directly. This test checks each function. 
- [t8_gtest_cache_competence.cxx](t8_gtest_cache_competence.cxx) tests that all predefined caching competences work as intended.
- [t8_gtest_ghost.cxx](t8_gtest_ghost.cxx) checks ghosts and the neighbor algorithm. Furthermore tests if all functions work also for ghost cells if applicable.
- [t8_gtest_handle_data.cxx](t8_gtest_handle_data.cxx) defines tests for user and element data.
- [t8_gtest_adapt.cxx](t8_gtest_adapt.cxx) checks that the adapt routine works as intended and compares it with the result of an adapted forest.

