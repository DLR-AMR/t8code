# t8_schemes Tests

This folder contains tests for the `t8_schemes` module. In `t8code`, a "scheme" defines the connectivity and refinement rules for a specific element shape (e.g., tetrahedrons, prisms, pyramids). These tests verify the correctness of these complex, low-level algorithms that form the combinatorial backbone of `t8code`.

## Files

The tests in this directory cover various aspects of the refinement schemes:

### Tree Traversal and Relationships

- `t8_gtest_ancestor.cxx`: Tests finding the ancestor of an element.
- `t8_gtest_ancestor_id.cxx`: Tests finding the ancestor of an element by its ID.
- `t8_gtest_descendant.cxx`: Tests finding descendants of an element.
- `t8_gtest_face_descendant.cxx`: Tests finding descendants of a face.
- `t8_gtest_find_parent.cxx`: Tests finding the direct parent of an element.
- `t8_gtest_nca.cxx`: Tests finding the nearest common ancestor of two elements.
- `t8_gtest_root.cxx`: Tests identifying the root element of a tree.
- `t8_gtest_elements_are_family.cxx`: Verifies the logic to determine if elements share a direct lineage.
- `t8_gtest_successor.cxx`: Tests finding the successor element in a traversal.
- `t8_gtest_dfs_base.hxx` & `t8_gtest_bfs_base.hxx`: Base headers for Depth-First Search and Breadth-First Search traversal tests.

### Connectivity and Neighborhood

- `t8_gtest_child_parent_face.cxx`: Tests the relationship between a child element's face and its parent's face.
- `t8_gtest_face_corner.cxx`: Tests the mapping between faces and corners of an element.
- `t8_gtest_face_neigh.cxx`: Tests finding face-neighbors of an element.
- `t8_gtest_pyra_connectivity.cxx`: Specific connectivity tests for pyramid elements.

### Scheme Properties and Consistency

- `t8_gtest_scheme_consistency.cxx`: A crucial test that checks the internal consistency of the refinement schemes.
- `t8_gtest_default.cxx`: Tests the default scheme.
- `t8_gtest_equal.cxx`: Tests the equality comparison for scheme elements.
- `t8_gtest_pack_unpack.cxx`: Verifies the serialization (packing) and deserialization (unpacking) of scheme data.
- `t8_gtest_input_equal_output.cxx`: Checks that certain operations are idempotent or reversible.
- `t8_gtest_boundary_extrude.cxx`: Tests extrusion of boundary faces.

### Element Identification and Properties

- `t8_gtest_get_linear_id.cxx`: Tests retrieving the linear ID of an element.
- `t8_gtest_set_linear_id.cxx`: Tests setting the linear ID of an element.
- `t8_gtest_element_count_leaves.cxx`: Tests counting the number of leaf elements.
- `t8_gtest_element_ref_coords.cxx`: Tests getting the reference coordinates of an element.

## Main Purpose

These tests are fundamental to the reliability of `t8code`. They ensure that the rules for how elements are created, related, and subdivided are mathematically correct and consistently applied, which is essential for the integrity of the adaptive mesh.

## How to Run

These tests are part of the main test suite and are executed when you run the `t8_gtest_main` executable.
