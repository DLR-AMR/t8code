# Advection Solver Example

This folder contains a complete, stand-alone finite volume advection solver. It is one of the most comprehensive examples in `t8code` and serves as an excellent case study for how to build a real application on top of the library.

The solver advects a scalar quantity through a domain, demonstrating many of `t8code`'s key features in a practical context.

## Files

- `t8_advection.cxx`: The main source code for the advection solver. This file shows how to:
    - Initialize `t8code` and a `t8_cmesh`.
    - Create a `t8_forest` and adapt it based on a given feature (in this case, the advected scalar field).
    - Attach user data (the scalar quantity) to the mesh elements.
    - Use the `t8_ghost` mechanism to exchange data between MPI processes.
    - Implement a numerical scheme (finite volume) using the `t8_forest` iterators.
    - Write out visualization files (VTK) at regular intervals.
- `t8_advection_generate_channel.geo`: A Gmsh geometry file for generating a 3D channel mesh. This can be used as input to the advection solver.
- `t8_advection_generate_channel_2d.geo`: A Gmsh geometry file for generating a 2D channel mesh.

## Main Features Demonstrated

- **Adaptive Mesh Refinement (AMR):** The solver adaptively refines the mesh to better resolve the features of the advected scalar field.
- **Parallel Execution:** The example is fully parallelized with MPI and uses `t8code`'s ghost layer functionality for communication.
- **Hybrid Meshes:** It can run on meshes with different element types (tetrahedra, hexahedra, prisms, etc.).
- **User Data Management:** It shows how to associate application-specific data with each element in the forest.
- **File I/O:** It can take a mesh file as input and writes out VTK files for visualization.

## How to Run

After building the examples, the `t8_advection` executable will be available. You can run it with the `-h` flag to see all available command-line options.

```bash
./t8_advection -h
```

A typical use case is to run it on a mesh file. You can generate a mesh using the provided `.geo` files with Gmsh.

```bash
# First, generate a mesh with Gmsh
gmsh t8_advection_generate_channel.geo -3 -o channel.msh

# Then, run the advection solver on that mesh
mpirun -n 4 ./t8_advection --mesh-file channel.msh
```

The solver will produce a series of `.pvtu` and `.vtu` files that can be opened in ParaView or VisIt to visualize the advection process.
