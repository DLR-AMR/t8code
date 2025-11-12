# t8_vtk Module

The `t8_vtk` module is responsible for writing `t8code` data, particularly the `t8_forest` and associated user data, into the Visualization Toolkit (VTK) file format. This enables users to visualize their adapted meshes and simulation results using popular visualization tools like [ParaView](https://www.paraview.org/) and [VisIt](https://visit-dav.github.io/visit-website/).

The module can write both serial (`.vtu`) and parallel (`.pvtu` + multiple `.vtu`) files.

## Subfolders

- [**t8_with_vtk/**](t8_with_vtk/README.md): Contains code that is only compiled when `t8code` is configured with VTK support. This code likely handles the direct interaction with the VTK library API.

## Key Files and Functionality

- **Main Writer Interface:**
    - `t8_vtk_writer.h`: The main public C-API header for the VTK writer. It declares the functions that a user calls to initiate the writing process, such as `t8_vtk_write_forest`.
    - `t8_vtk_writer.hxx`: The internal C++ header for the writer class.
    - `t8_vtk_writer.cxx`: The main implementation of the VTK writer logic. It orchestrates the process of gathering the mesh geometry and data and writing it to files.

- **ASCII Format Implementation:**
    - `t8_vtk_write_ASCII.hxx`/`.cxx`: These files contain the specific implementation for writing VTK files in the human-readable ASCII format. The VTK standard also supports a binary format, which would be implemented separately.

- **Helper Functions:**
    - `t8_vtk_writer_helper.hxx`/`.cxx`: These files contain helper functions that assist the main writer. This could include functions for triangulating non-tetrahedral elements (as VTK's unstructured grid format is based on simple cell types), collecting data from ghost elements, or managing parallel file information.

- **Type Definitions:**
    - `t8_vtk_types.h`: Defines data types and structures specific to the VTK writing process.

## Developer Workflow

A developer wanting to visualize their `t8_forest` and data would typically:
1.  Include `t8_vtk_writer.h`.
2.  After their simulation step, call `t8_vtk_write_forest`, passing in the `t8_forest` object.
3.  They can also provide pointers to their user data (managed by a `t8_data_handler`) to be included in the output file. The function allows specifying names for these data arrays.
4.  The module will then generate the appropriate `.vtu` or `.pvtu` files on disk.

This provides a straightforward way to get simulation data out of `t8code` and into a standard format for visualization and analysis.
