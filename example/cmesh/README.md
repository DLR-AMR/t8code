# Coarse Mesh (cmesh) Examples

This folder contains examples that demonstrate how to create, manipulate, and query the `t8_cmesh` object. The `t8_cmesh` is the coarse, unadapted representation of the domain and serves as the basis for the adaptive `t8_forest`.

## Files

- `t8_cmesh_create_partitioned.cxx`: This example shows how to create a `t8_cmesh` that is already partitioned across multiple MPI processes. This is an advanced use case for situations where the coarse mesh itself is too large to fit on a single node.
- `t8_cmesh_geometry_examples.cxx`: Demonstrates how to associate different types of geometries with a `t8_cmesh`. This is key to running simulations on non-trivial, curved domains.
- `t8_cmesh_hypercube_pad.cxx`: An example of how to create a hypercube-shaped cmesh with padding, which can be useful for certain boundary condition treatments.
- `t8_cmesh_partition.cxx`: This program shows the standard workflow for partitioning a `t8_cmesh` that is initially loaded on a single process and then distributed across all available MPI processes.
- `t8_cmesh_set_join_by_vertices.cxx`: An example that demonstrates how to control the connectivity of the coarse mesh by defining how faces are joined based on their shared vertices.

## Main Purpose

These examples provide a guide to the various ways a `t8_cmesh` can be constructed and configured. Understanding these concepts is fundamental to setting up a simulation in `t8code`, as the `t8_cmesh` is the first major object that needs to be created.

## How to Run

After building the examples, you can run the individual executables. Some may be intended to be run with MPI.

```bash
# Run a simple example on a single process
./t8_cmesh_hypercube_pad

# Run a partitioning example with multiple processes
mpirun -n 4 ./t8_cmesh_partition
```
