# Forest I/O Examples

This folder contains examples related to Input/Output operations for the `t8_forest` object. These examples demonstrate how to save and load the state of an adaptive forest, which is a key feature for checkpointing and restarting simulations.

## Subfolders

- [**gmsh/**](gmsh/README.md): An example of how to save a `t8_forest` into the Gmsh `.msh` format.
- [**netcdf/**](netcdf/README.md): An example of how to save a `t8_forest` into the NetCDF format.

## Main Purpose

These examples are crucial for users who need to save the results of a simulation or checkpoint its state. Being able to save the adapted forest allows a simulation to be restarted exactly where it left off, or for the final adapted mesh to be analyzed or used in subsequent computations.
