# t8_forest_incomplete Tests

This folder contains tests for handling "incomplete" forests. In `t8code`, an incomplete forest is one that may have unusual or non-standard properties, such as containing empty trees or having "holes" in its structure. These tests are designed to ensure the robustness of forest algorithms when faced with these edge cases.

## Files

- `t8_gtest_empty_global_tree.cxx`: Tests how the forest handles a cmesh tree that is globally empty (i.e., has no elements on any process).
- `t8_gtest_empty_local_tree.cxx`: Tests how the forest handles a cmesh tree that is locally empty (i.e., has no elements on the current process, but may have elements on other processes).
- `t8_gtest_iterate_replace.cxx`: Tests the iteration over a forest while simultaneously replacing elements, a complex scenario that can arise during adaptation.
- `t8_gtest_permute_hole.cxx`: Tests the handling of a "permuted hole," which likely refers to a scenario where a gap or missing element is created in the forest due to permutation or reordering operations.
- `t8_gtest_recursive.cxx`: Contains recursive tests, likely to probe deep adaptation scenarios or complex tree structures.

## Main Features

The tests in this directory are crucial for ensuring the stability and correctness of `t8code`'s core forest algorithms. By specifically targeting edge cases and unusual configurations, these tests help prevent bugs that might not be caught by standard adaptation and partitioning tests.

## How to Run

These tests are part of the main test suite and are executed when you run the `t8_gtest_main` executable.
