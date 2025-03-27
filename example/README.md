The programs here are a collection of examples on the usage of t8code.
A brief description of the folders and programs:

- `advect` -- A finite volume advection solver. This solver either computes on a unit cube, or on a user provided gmsh mesh. It works in 2D and 3D with all supported element types and hybrid meshes. Call `./t8_advection -h` for usage details.

- `IO` -- Options to read and write forests and cmeshes.

- `forest, cmesh, geometry` -- Advanced usage examples for the main t8code objects.

- `common` -- Common functions used in the examples.

- `version` -- A program to illustrate version functionalities of t8code.