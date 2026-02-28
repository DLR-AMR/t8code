# t8_forest Tests

This folder contains tests for the `t8_forest` module, which implements the core functionality of `t8code`: the adaptive forest of octrees (or other element types). This includes mesh adaptation (refinement and coarsening), load balancing, ghost element exchange, and search operations.

## Files

- `t8_gtest_balance.cxx`: Tests the 2:1 balancing of the forest, which is a crucial step after adaptation to ensure a valid mesh.
- `t8_gtest_element_is_leaf.cxx`: Tests the function that determines if an element is a leaf (i.e., not refined further).
- `t8_gtest_element_volume.cxx`: Verifies the calculation of element volumes.
- `t8_gtest_find_owner.cxx`: Tests the logic for finding the MPI process that owns a given element.
- `t8_gtest_forest_commit.cxx`: Tests the `t8_forest_commit` function, which finalizes changes to the forest after adaptation or partitioning.
- `t8_gtest_forest_face_normal.cxx`: Tests the computation of face normals for elements in the forest.
- `t8_gtest_ghost_and_owner.cxx`: Tests the identification of ghost elements and their owners.
- `t8_gtest_ghost_delete.cxx`: Verifies the correct deletion of ghost elements.
- `t8_gtest_ghost_exchange.cxx`: Tests the exchange of ghost element data between processes, which is essential for parallel computations.
- `t8_gtest_half_neighbors.cxx`: Tests the retrieval of "half" neighbors, which are neighbors of the same size.
- `t8_gtest_partition_data.cxx`: Tests the partitioning of the forest data for load balancing.
- `t8_gtest_search.cxx`: Tests the functions for searching for elements within the forest (e.g., finding which element contains a given point).
- `t8_gtest_transform.cxx`: Tests coordinate transformations between the reference element and the physical element.
- `t8_gtest_user_data.cxx`: Tests the attachment and management of user data on forest elements.

## Main Features

These tests ensure the correctness and robustness of the dynamic mesh adaptation and parallel data management capabilities of `t8code`, which are fundamental to its use in AMR simulations.

## How to Run

These tests are part of the main test suite and are executed when you run the `t8_gtest_main` executable.
