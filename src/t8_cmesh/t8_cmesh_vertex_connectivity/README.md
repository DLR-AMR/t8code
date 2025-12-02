# cmesh Vertex Connectivity

This sub-module of `t8_cmesh` is responsible for building and querying explicit connectivity information based on the vertices of the coarse mesh. This information is not always needed by default, but it is crucial for certain algorithms that require knowledge about which elements share a common vertex.

## Key Data Structures and Functionality

The core of this module is the `t8_cmesh_vertex_connectivity_t` data structure, which stores the relationships between vertices and the coarse mesh trees (elements).

- `t8_cmesh_vertex_connectivity_types.hxx`: Defines the `t8_cmesh_vertex_connectivity_t` struct and other related types.
- `t8_cmesh_vertex_connectivity.cxx`/`.h`/`.hxx`: Contains the main functions to build, query, and destroy the vertex connectivity data structure. The `t8_cmesh_vertex_connectivity_new` function is the main entry point, which takes a `t8_cmesh` and constructs the full connectivity information.

The connectivity information is stored in two primary forms:

1.  **Tree-to-Vertex Mapping:** This allows for quickly finding the global vertex indices that make up a given coarse mesh tree.
    - `t8_cmesh_vertex_conn_tree_to_vertex.cxx`/`.hxx`: Implements the logic to build and access this part of the connectivity data.

2.  **Vertex-to-Tree Mapping:** This is the inverse mapping, which allows for finding all coarse mesh trees that share a specific vertex. This is often the more complex and useful part of the data structure.
    - `t8_cmesh_vertex_conn_vertex_to_tree.cxx`/`.hxx`: Implements the logic to build and access this inverse mapping.

## Main Purpose

The main purpose of this module is to provide answers to questions like:
- "Which elements are neighbors of element X around vertex Y?"
- "Which elements share this specific vertex?"

This information is vital for some numerical schemes, mesh quality checks, and advanced partitioning algorithms that operate on the coarse mesh. Since building this data structure has a computational cost, it is handled in a separate module and only computed when explicitly requested by the user or another algorithm.

## Developer Entry Points

A developer would typically call `t8_cmesh_vertex_connectivity_new` to build the connectivity data for a given `t8_cmesh`. They could then use the query functions provided in `t8_cmesh_vertex_connectivity.h` to access the connectivity information. Finally, they would call `t8_cmesh_vertex_connectivity_destroy` to free the associated memory.
