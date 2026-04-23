# t8_geometry Tests

This folder contains tests for the geometry handling capabilities of `t8code`. These tests ensure that the mesh can correctly interact with and conform to various geometric representations, which is essential for simulating problems on complex domains.

## Subfolders

- [t8_geometry_implementations/](t8_geometry_implementations/README.md): Contains tests for specific geometry implementations, such as analytic geometries or those based on CAD kernels.

## Files

- `t8_gtest_geometry_handling.cxx`: Provides general tests for geometry handling, verifying the core functionalities of the geometry module.
- `t8_gtest_geometry_triangular_interpolation.cxx`: Tests the triangular interpolation functions, which are often used to represent and evaluate curved surfaces.
- `t8_gtest_point_inside.cxx`: Tests the crucial "point inside element" check, which determines if a given point lies within a mesh element, taking into account the underlying geometry.

## Main Features

The tests in this directory verify the accuracy of geometric queries and transformations. This includes ensuring that mesh elements correctly represent the underlying geometry and that operations like point location and interpolation are performed correctly.

## How to Run

These tests are part of the main test suite and are executed when you run the `t8_gtest_main` executable. Some tests may require specific geometry kernels (like OpenCASCADE) to be enabled in the CMake configuration.
