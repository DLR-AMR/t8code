# Test Folder

This folder contains all tests for the t8code library. The tests are written using the Google Test framework.

## Subfolders

The tests are organized into subfolders, each corresponding to a specific module or functionality of t8code.

- [api/](api/README.md): Tests for the public C-API of t8code.
- [t8_IO/](t8_IO/README.md): Tests for input/output operations.
- [t8_cmesh/](t8_cmesh/README.md): Tests for the cmesh module, which handles the coarse mesh.
- [t8_cmesh_generator/](t8_cmesh_generator/README.md): Tests for the generation of cmeshes.
- [t8_data/](t8_data/README.md): Tests for data handling and management.
- [t8_forest/](t8_forest/README.md): Tests for the adaptive forest of octrees.
- [t8_forest_incomplete/](t8_forest_incomplete/README.md): Tests for incomplete forest operations.
- [t8_geometry/](t8_geometry/README.md): Tests for geometry handling.
- [t8_schemes/](t8_schemes/README.md): Tests for different numerical schemes.
- [t8_types/](t8_types/README.md): Tests for custom data types used in t8code.
- [t8_vector_helper/](t8_vector_helper/README.md): Tests for vector helper functions.
- [testfiles/](testfiles/README.md): Contains files used by the tests, such as mesh files.

## Main Entry Points

The main entry point for running the tests is `t8_gtest_main.cxx`. The tests are built and run using CMake. To build the tests, enable the `T8_ENABLE_TESTING` option in your CMake configuration.

Individual tests are defined in `_test.cxx` files within the subfolders. The `t8_gtest_...` files in this directory provide the basic setup and custom assertions for the tests.

## Related Wiki Articles

For more information on how to build and run the tests, please refer to the [t8code wiki](https://github.com/DLR-AMR/t8code/wiki).
