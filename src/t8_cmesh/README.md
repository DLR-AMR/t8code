# t8_cmesh Module

The `t8_cmesh` module is responsible for managing the coarse mesh, which is the initial, unadapted, and partitioned representation of the simulation domain. The `t8_cmesh_t` data structure, defined here, stores the geometry and connectivity of the coarse elements (e.g., tetrahedra, hexahedra) and how they are distributed among MPI processes.

The `t8_cmesh` serves as the foundation upon which the adaptive `t8_forest` is built.

## Subfolders
- [**t8_cmesh_vertex_connectivity/**](t8_cmesh_vertex_connectivity/README.md): Contains the implementation for building and querying vertex-based connectivity information for the coarse mesh.

## Key Files and Functionalities

The functionality of the `t8_cmesh` module is split across several files:

- **Core Cmesh Object:**
    - `t8_cmesh_types.h`: Defines the `t8_cmesh_t` struct and other related data types.
    - `t8_cmesh.cxx`: Implements the core functions for creating, querying, and destroying a `t8_cmesh`.

- **I/O and Mesh Generation:**
    - `t8_cmesh_readmshfile.cxx`: Handles reading coarse meshes from Gmsh `.msh` files.
    - `t8_cmesh_triangle.cxx`: Contains functionality for reading 2D meshes from the Triangle mesh generator.
    - `t8_cmesh_netcdf.c`: Implements reading/writing of cmesh data from/to NetCDF files.
    - `t8_cmesh_examples.cxx`/`.h`: Provides functions to generate various example cmeshes programmatically.
    - `t8_cmesh_save.cxx`/`.h`: Implements functionality for saving a `t8_cmesh` to a file.

- **Partitioning and Distribution:**
    - `t8_cmesh_partition.cxx`/`.h`: Implements the logic for partitioning a serial cmesh and distributing it across MPI processes.
    - `t8_cmesh_offset.c`/`.h`: Manages the offsets used to describe data distribution in parallel.

- **Geometry and CAD Interface:**
    - `t8_cmesh_geometry.cxx`/`.hxx`: Contains the logic for associating a `t8_geometry` object with the cmesh.
    - `t8_cmesh_cad.cxx`/`.hxx`: Implements the specific interface for handling geometries from CAD kernels.

- **Internal Operations:**
    - `t8_cmesh_commit.cxx`: Implements the `t8_cmesh_commit` function, which finalizes the cmesh after it has been created and partitioned.
    - `t8_cmesh_copy.c`/`.h`: Provides functions to create a copy of a `t8_cmesh`.
    - `t8_cmesh_trees.cxx`/`.h`: Manages the trees of the cmesh. Each tree corresponds to a single element in the original, unpartitioned mesh.
    - `t8_cmesh_helpers.cxx`/`.h`: A collection of utility and helper functions for internal use.
    - `t8_cmesh_stash.c`/`.h`: Implements a "stash" data structure, likely used for temporarily storing and communicating mesh data during construction.

## Main Entry Points for Developers

A developer working with the `t8_cmesh` would typically start by:
1. Creating a `t8_cmesh` instance, either from a file (`t8_cmesh_from_msh_file`) or programmatically (`t8_cmesh_new_...`).
2. Finalizing the mesh with `t8_cmesh_commit`.
3. Using the committed `t8_cmesh` to build a `t8_forest`.

The public API for these operations is primarily found in `src/t8_cmesh.h`. Internal implementation details are in the files listed above.
