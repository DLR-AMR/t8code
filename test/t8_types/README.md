# t8_types Tests

This folder contains tests for the custom data types and data structures defined in `t8code`. These tests ensure that the fundamental building blocks for data storage and access are working correctly.

## Files

- `t8_gtest_type.cxx`: This file tests the basic custom data types used throughout `t8code`. This might include tests for their size, alignment, and basic operations to ensure platform-independent and consistent behavior.
- `t8_gtest_random_accessible.cxx`: This file tests data structures that are designed to provide random access to their elements. It likely verifies that elements can be accessed correctly and efficiently by their index.

## Main Purpose

The tests in this directory are fundamental to the stability of the entire `t8code` library. They verify the correctness of the low-level data types and structures upon which all the higher-level modules (like `t8_cmesh` and `t8_forest`) are built.

## How to Run

These tests are part of the main test suite and are executed when you run the `t8_gtest_main` executable.
