# VTK Library Interaction

This directory contains the source code that directly interacts with the third-party [VTK library](https://vtk.org/). The code in this folder is only compiled if `t8code` is configured with VTK support enabled (`-DT8_WITH_VTK=ON`).

The separation of this code allows the rest of the `t8_vtk` module to have a generic structure, while the files here handle the specific details of calling VTK API functions to create and populate VTK data structures.

## Key Files and Functionality

-   **Unstructured Grid Handling:**
    -   `t8_vtk_unstructured.cxx`/`.hxx`: This is a crucial component for visualization. It contains the code that takes the `t8_forest` data (points and element connectivity) and populates a `vtkUnstructuredGrid` object. This VTK data structure is capable of representing the mix of element types (tets, hexes, etc.) found in a `t8_forest`.

-   **PolyData Handling:**
    -   `t8_vtk_polydata.cxx`/`.hxx`: This contains code to create `vtkPolyData` objects. `vtkPolyData` is typically used for representing surface meshes (vertices and polygons). This might be used in `t8code` to visualize the boundary of a 3D domain or to represent a 2D mesh.

-   **Parallel I/O:**
    -   `t8_vtk_parallel.cxx`/`.hxx`: This component handles the specifics of writing parallel VTK files. It likely uses the VTK `vtkXMLPUnstructuredGridWriter` (for `.pvtu` files) or similar classes to coordinate the writing process across multiple MPI ranks, where each rank writes its own piece of the data.

-   **VTK File Reading:**
    -   `t8_vtk_reader.cxx`/`.hxx`: This contains the implementation for reading mesh data *from* a VTK file. It uses the VTK library's reader classes (e.g., `vtkXMLUnstructuredGridReader`) to parse a `.vtu` file and extract the mesh information, which can then be used to initialize a `t8_cmesh`.

## Main Purpose

The main purpose of this module is to encapsulate all direct dependencies on the VTK library API. By doing so, the rest of `t8code` remains independent of VTK, and `t8code` can be successfully compiled even if the VTK library is not installed on a user's system (albeit without visualization capabilities). When VTK support is enabled, these files provide the "bridge" between `t8code`'s internal data structures and the corresponding VTK data structures.
