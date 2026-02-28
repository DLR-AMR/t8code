# t8_data Tests

This folder contains tests for data handling and management within t8code. This includes tests for the data handler, which is responsible for managing user-defined data attached to mesh elements, as well as tests for shared memory operations.

## Files

- `t8_gtest_data_handler.cxx`: The main test file for the `t8_data_handler`. It verifies that user data can be correctly allocated, stored, and accessed on a per-element basis.
- `t8_gtest_shmem.cxx`: Contains tests for shared memory (`shmem`) functionalities, ensuring that data can be efficiently shared between processes on the same node.
- `t8_data_handler_specs.cxx` and `t8_data_handler_specs.hxx`: These files define specifications and configurations for the data handler tests, allowing for a variety of test scenarios.
- `t8_pseudo_trees.hxx`: A header file defining "pseudo-trees," which are simplified tree structures used for testing data management without needing a full forest.
- `t8_enlarged_stdtypes.hxx`: A header file providing enlarged standard data types, likely used to test data handling with non-standard data sizes.

## Main Features

The tests in this directory are crucial for verifying the correctness and efficiency of how `t8code` manages user data. This is a core feature for any simulation code, as it allows users to attach their application-specific data (like physical quantities) to the mesh.

## How to Run

These tests are part of the main test suite and are executed when you run the `t8_gtest_main` executable.
