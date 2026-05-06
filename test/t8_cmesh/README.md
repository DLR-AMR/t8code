# t8_cmesh Tests

This folder contains tests for the `t8_cmesh` module, which is responsible for handling the coarse mesh in `t8code`. The coarse mesh is the initial, non-adapted representation of the simulation domain.

## Files

The test files in this directory cover a wide range of `t8_cmesh` functionalities:

- `t8_gtest_attribute_gloidx_array.cxx`: Tests the handling of global index attributes on the coarse mesh.
- `t8_gtest_bcast.cxx`: Tests the broadcasting of cmesh data across MPI processes.
- `t8_gtest_cmesh_add_attributes_when_derive.cxx`: Verifies that attributes are correctly handled when a new cmesh is derived.
- `t8_gtest_cmesh_bounding_box.cxx`: Tests the computation of the cmesh's bounding box.
- `t8_gtest_cmesh_copy.cxx`: Tests the functionality for copying a cmesh.
- `t8_gtest_cmesh_face_is_boundary.cxx`: Tests the identification of boundary faces on the cmesh.
- `t8_gtest_cmesh_partition.cxx`: Verifies the partitioning of the cmesh across multiple processes.
- `t8_gtest_cmesh_readmshfile.cxx`: Tests reading a cmesh from a `.msh` file.
- `t8_gtest_cmesh_set_join_by_vertices.cxx`: Tests the logic for joining mesh entities by vertices.
- `t8_gtest_cmesh_set_partition_offsets.cxx`: Tests the setting of partition offsets.
- `t8_gtest_cmesh_tree_vertices_negative_volume.cxx`: A specific test to check for negative volume elements, which can indicate errors in mesh generation or adaptation.
- `t8_gtest_cmesh_vertex_conn.cxx`: Tests the vertex connectivity data structures.
- `t8_gtest_cmesh_vertex_conn_tree_to_vertex.cxx`: Tests the connectivity from trees to vertices.
- `t8_gtest_cmesh_vertex_conn_vertex_to_tree.cxx`: Tests the connectivity from vertices to trees.
- `t8_gtest_compute_first_element.cxx`: Tests the computation of the first element on a process.
- `t8_gtest_hypercube.cxx`: Contains tests specifically for hypercube-shaped cmeshes.
- `t8_gtest_multiple_attributes.cxx`: Tests the handling of multiple attributes on a cmesh.

## Main Features

These tests are crucial for ensuring the correctness of the coarse mesh representation and its manipulation, which forms the foundation for the adaptive forest of octrees in `t8code`.

## How to Run

These tests are part of the main test suite and are executed when you run the `t8_gtest_main` executable.
