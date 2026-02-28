# t8_cmesh_generator Tests

This folder contains tests for the `t8_cmesh_generator` module, which is responsible for programmatically generating various types of coarse meshes.

## Subfolders

- [t8_cmesh_parameterized_examples/](t8_cmesh_parameterized_examples/README.md): Contains parameterized examples of cmeshes used for testing.

## Files

- `t8_gtest_cmesh_generator_test.cxx`: This is the main test file for the cmesh generator. It verifies that the generator can create valid and correct coarse meshes based on different input parameters and configurations.
- `t8_cmesh_example_sets.hxx`: A header file defining sets of example cmeshes that are used as a basis for the generator tests.
- `t8_gtest_cmesh_cartestian_product.hxx`: A header file providing utilities to test the generation of cmeshes from Cartesian products of simpler meshes.
- `t8_gtest_cmesh_sum_of_sets.hxx`: A header file providing utilities to test the generation of cmeshes by combining different mesh sets.

## Main Features

The tests in this directory ensure that `t8code` can reliably generate a variety of coarse meshes, which is a fundamental capability for setting up different simulation scenarios without relying on externally created mesh files.

## How to Run

These tests are part of the main test suite and are executed when you run the `t8_gtest_main` executable.
