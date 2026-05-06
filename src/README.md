# t8code Source

This directory contains the core source code of the `t8code` library. The code is organized into modules, each corresponding to a specific functionality.

## Subfolders (Modules)

- [**t8_cmesh/**](t8_cmesh/README.md): The coarse mesh module. It is responsible for storing and managing the initial, unadapted mesh.
- [**t8_data/**](t8_data/README.md): Data management module, including the `t8_data_handler` for user-defined data on elements.
- [**t8_forest/**](t8_forest/README.md): The core module of `t8code`. It implements the adaptive forest of trees, handling mesh refinement, coarsening, partitioning, and ghost element exchange.
- [**t8_geometry/**](t8_geometry/README.md): The geometry module, which provides an interface for handling various geometric representations of the domain.
- [**t8_schemes/**](t8_schemes/README.md): This module defines the refinement and connectivity rules for different element types (e.g., tetrahedra, hexahedra).
- [**t8_types/**](t8_types/README.md): Contains definitions for custom data types and structures used throughout `t8code`.
- [**t8_vector_helper/**](t8_vector_helper/README.md): Provides helper functions for vector manipulations.
- [**t8_vtk/**](t8_vtk/README.md): Contains code for writing data in the VTK format for visualization.

## Main Files

In addition to the modules in the subfolders, this directory contains some of the main header files and source files for the library:

- `t8.h`: The main public header file for the `t8code` C API. This is the primary entry point for users of the library.
- `t8.c`: Implementation of some of the general functions declared in `t8.h`.
- `t8_cmesh.h`, `t8_forest.h`, etc.: Public header files for the specific modules.
- `t8_eclass.h`, `t8_element.h`: Header files defining the element classes and element data structures.

## Developer Entry Points

For developers looking to understand the library, the main entry points are:
- `t8.h`: To understand the public API.
- `t8_forest.h`: To understand the main `t8_forest` data structure and its functions.
- The `README.md` file in each subdirectory for a detailed description of that module.

## Related Wiki Articles

For a higher-level overview of the architecture and concepts, please refer to the [t8code wiki](https://github.com/DLR-AMR/t8code/wiki).
