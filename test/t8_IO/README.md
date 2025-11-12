# t8_IO Tests

This folder contains tests for the input/output functionalities of t8code, specifically focusing on VTK file operations.

## Files

- `t8_gtest_vtk_reader.cxx`: Contains tests for reading VTK files. It verifies that t8code can correctly parse and load data from VTK file formats.
- `t8_gtest_vtk_writer.cxx`: Contains tests for writing VTK files. It ensures that t8code can correctly export data into various VTK file formats.

## Main Features

The tests in this directory ensure the reliability of t8code's I/O operations, which are crucial for data import/export and interoperability with other tools that use the VTK format, such as ParaView or VisIt.

## How to Run

These tests are part of the main test suite and are executed when you run the `t8_gtest_main` executable. Ensure that VTK support is enabled in the CMake configuration to build and run these tests.
