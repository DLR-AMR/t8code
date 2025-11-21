# Gmsh I/O Examples

This folder contains examples that demonstrate how to read `.msh` files, the native format of the open-source mesh generator [Gmsh](https://gmsh.info/).

## Files

- `t8_read_msh_file.cxx`: A straightforward example that shows the basic process of reading a `.msh` file and initializing a `t8_cmesh` from it.
- `t8_load_and_refine_square_w_hole.cxx`: A more complex example that loads a hybrid 2D mesh (containing both quads and triangles) of a square with a circular hole. After loading, it uniformly refines the resulting forest to demonstrate a complete workflow.
- `circlesquare_hybrid_hole.msh`: The Gmsh mesh file used as input for the `t8_load_and_refine_square_w_hole` example.

## Main Purpose

These examples provide a practical guide for users who want to use `t8code` to process meshes created with Gmsh. This is a very common use case, as Gmsh is a powerful and widely used tool for generating unstructured meshes.

## How to Run

After building the examples, you can run the executables from the build directory. For example:

```bash
# Run the basic reader example
./t8_read_msh_file <path_to_your_mesh.msh>

# Run the load-and-refine example
./t8_load_and_refine_square_w_hole <path_to_circlesquare_hybrid_hole.msh>
```
