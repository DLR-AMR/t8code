# Gmsh Forest I/O Example

This folder contains an example of a workflow involving a `t8_forest` and the Gmsh `.msh` file format.

## Files

- `t8_gmsh_to_vtk.cxx`: This program demonstrates a common use case: loading a mesh from a Gmsh file, creating a `t8_forest` from it, and then writing the resulting forest data to a VTK file for visualization. While the primary action is writing a VTK file, the initial step involves reading a `.msh` file, showcasing `t8code`'s ability to ingest Gmsh files as a starting point for forest creation.

## Main Purpose

This example is useful for users who have meshes in the Gmsh format and want to use them in `t8code` and then visualize the results. It bridges the gap between mesh generation with Gmsh and analysis/visualization with tools that support the VTK format (like ParaView or VisIt).

## How to Run

After building the examples, you can run this program by providing a path to a Gmsh `.msh` file.

```bash
./t8_gmsh_to_vtk <path_to_your_mesh.msh>
```

This will produce a VTK file (e.g., `output.vtk`) in the execution directory.
