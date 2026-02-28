# cmesh I/O Examples

This folder provides examples on how to perform Input/Output operations with the `t8_cmesh` object. The main example demonstrates a general load-and-save capability.

## Files

- `t8_cmesh_load_save.cxx`: This program shows how to load a `t8_cmesh` from a file, and then save it back to another file. This is a useful pattern for converting between mesh formats or for inspecting a mesh.

## Subfolders

The subfolders contain examples that are specific to different mesh file formats or mesh generators that `t8code` can interface with.

- [**gmsh/**](gmsh/README.md): Examples for reading `.msh` files from Gmsh.
- [**netcdf/**](netcdf/README.md): Examples for reading meshes stored in the NetCDF format.
- [**tetgen/**](tetgen/README.md): Examples for reading meshes from TetGen.
- [**triangle/**](triangle/README.md): Examples for reading 2D meshes from Triangle.
- [**vtk/**](vtk/README.md): Examples for writing `t8_cmesh` data to VTK files for visualization.

## Main Purpose

These examples are designed to help developers understand how to get coarse mesh data into and out of `t8code`, which is a fundamental first step for many simulation workflows.
