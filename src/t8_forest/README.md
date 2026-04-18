# t8_forest Module

The `t8_forest` module is the heart of `t8code`. It implements the core data structure, `t8_forest_t`, which represents a collection (a "forest") of adaptive refinement trees. Each tree is built upon a coarse element from the `t8_cmesh`. This module contains the logic for all the fundamental operations of a dynamic, parallel adaptive mesh refinement (AMR) framework.

## Subfolders
- [**t8_forest_search/**](t8_forest_search/README.md): Contains the implementation of algorithms for searching within the forest, such as finding the process that owns a specific point in space.

## Key Files and Functionalities

The `t8_forest` is a complex module whose functionality is broken down into several logical components:

- **Core Forest Management:**
    - `t8_forest_types.h`: Defines the `t8_forest_t` struct, the central data structure that manages the entire adaptive mesh.
    - `t8_forest.h`: The main public C-API header for the forest. It contains the primary functions for creating, querying, and managing a `t8_forest`.
    - `t8_forest.cxx`: Implementation of the core forest management functions.
    - `t8_forest_private.h`/`.cxx`: Contains internal functions and data structures that are not part of the public API but are used by other forest components.

- **Dynamic Adaptation:**
    - `t8_forest_adapt.h`/`.cxx`: Implements the functions for mesh adaptation. This includes refining a set of elements (dividing them into children) and coarsening them (replacing children with their parent). This is where the `t8_forest_adapt` and `t8_forest_coarsen` functions are implemented.
    - `t8_forest_balance.h`/`.cxx`: Implements the 2:1 balancing algorithm. After adaptation, the forest must be "balanced" to ensure that no leaf element is adjacent to a non-parent element that is more than twice its size. This is essential for the validity of many numerical schemes.

- **Parallelism and Communication:**
    - `t8_forest_partition.h`/`.cxx`: Implements functions to re-partition the forest. When the mesh is adapted, the workload may become unbalanced across MPI processes. These functions redistribute the elements of the forest to re-balance the load.
    - `t8_forest_ghost.h`/`.cxx`: Implements the ghost layer functionality. It provides functions to build a layer of ghost elements along the boundaries of each process's subdomain and to exchange data stored on those elements. This is the primary mechanism for inter-process communication.

- **Iteration and Traversal:**
    - `t8_forest_iterate.h`/`.cxx`: Provides powerful and flexible iterators for traversing the elements and faces of the forest. These iterators are the main tool for developers implementing numerical solvers, as they provide access to the mesh entities and their neighborhood information.

- **I/O:**
    - `t8_forest_io.h`: Header file for I/O operations related to the forest.
    - `t8_forest_netcdf.cxx`: Implementation for saving a `t8_forest` to a NetCDF file.

- **Helper Headers:**
    - `t8_forest_general.h`: A collection of general-purpose forest functions.
    - `t8_forest_geometrical.h`: A collection of functions related to the geometry of the forest elements (e.g., computing volumes, face normals).
    - `t8_forest_profiling.h`: Contains tools for profiling the performance of forest operations.

## Developer Workflow

A typical developer workflow using the `t8_forest` module looks like this:

1.  **Initialization:** Create a `t8_forest` from a `t8_cmesh` using `t8_forest_new_from_cmesh`.
2.  **Adaptation Loop:**
    a.  **Solve:** Use the iterators (`t8_forest_iterate_...`) to perform calculations on the current mesh.
    b.  **Mark:** Decide which elements need to be refined or coarsened based on the solution.
    c.  **Adapt:** Call `t8_forest_adapt` and/or `t8_forest_coarsen` to modify the mesh.
    d.  **Balance:** Call `t8_forest_balance` to ensure the mesh is valid.
    e.  **Partition (Optional):** Call `t8_forest_partition` to re-balance the workload.
    f.  **Update Ghosts:** Re-build the ghost layer and exchange data for the new mesh configuration.
3.  **Finalization:** Destroy the forest using `t8_forest_unref`.

The main public API functions are declared in `t8_forest.h` and the other public header files in this directory.
