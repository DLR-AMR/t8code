# Specific Geometry Implementation Tests

This folder contains tests for different, specific implementations of the `t8_geometry` interface. Each file tests a particular way of representing and interacting with the domain geometry.

## Files

- `t8_gtest_geometry_cad.cxx`: Contains tests for geometries defined by a CAD kernel, such as OpenCASCADE. These tests verify that `t8code` can correctly interact with complex, industry-standard CAD models. Running these tests requires `t8code` to be compiled with support for a CAD library.
- `t8_gtest_geometry_lagrange.cxx`: Tests the implementation of high-order geometries using Lagrange polynomials. This is used for simulations that require a very accurate, curved representation of the domain boundary.
- `t8_gtest_geometry_linear.cxx`: Tests the simplest geometry implementation, where all surfaces are represented as flat, linear patches. This serves as a baseline and is used for many standard test cases.

## Main Purpose

These tests ensure that each supported geometry backend is working correctly. This is crucial for allowing users to choose the most appropriate geometry representation for their specific application, from simple analytical shapes to complex CAD models.

## How to Run

These tests are part of the main test suite. Note that some tests, particularly `t8_gtest_geometry_cad.cxx`, will only be built and executed if the corresponding optional dependencies (e.g., OpenCASCADE) are enabled in the CMake configuration.
