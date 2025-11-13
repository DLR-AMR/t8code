# Forest Examples

This folder contains examples that demonstrate advanced usage of the `t8_forest` object. The `t8_forest` is the main data structure in `t8code` that manages the adaptively refined mesh.

## Files

- `t8_test_face_iterate.cxx`: This example demonstrates how to iterate over the faces of elements in a `t8_forest`. This is a common requirement for numerical methods like Discontinuous Galerkin (DG) or Finite Volume (FV) methods, where fluxes across faces need to be computed.
- `t8_test_ghost.cxx`: A key example that shows how to use the ghost layer functionality. It demonstrates the process of creating a layer of "ghost" elements around the boundary of each MPI process's local domain and then exchanging data with the real owner of those elements. This is fundamental for parallel computations.
- `t8_test_ghost_large_level_diff.cxx`: This is a more specific and advanced version of the ghost test. It sets up a scenario where neighboring elements have a large difference in their refinement levels and then tests the ghost exchange. This ensures that the ghost communication is robust even in cases of extreme mesh adaptation.

## Main Purpose

These examples are crucial for developers who are implementing parallel numerical solvers with `t8code`. They provide clear demonstrations of how to perform two of the most important operations on a `t8_forest`: iterating over faces and communicating data between processes using the ghost layer.

## How to Run

After building the examples, the executables can be run. The ghost examples are intended to be run with MPI.

```bash
# Run the face iteration example
./t8_test_face_iterate

# Run the ghost test with multiple processes
mpirun -n 4 ./t8_test_ghost
```
