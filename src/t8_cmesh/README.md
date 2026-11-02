# src/t8_cmesh

This folder contains the coarse mesh interface.

The coarse mesh (t8_cmesh) contains information about the domain topology and geometry.
It is usually created from an input mesh, either from a mesh generator or by hand.

See also our Wiki tutorials [Creating a coarse mesh](https://github.com/DLR-AMR/t8code/wiki/Step-1---Creating-a-coarse-mesh)
and [Building a Cmesh by hand](https://github.com/DLR-AMR/t8code/wiki/Building-a-Cmesh-by-hand).

Your entry point to familiarize yourself with the code should be [t8_cmesh.h](t8_cmesh.h) which is the main cmesh interface.

This folder contains the following sub-directories:

#### [src/t8_cmesh/t8_cmesh_internal](t8_cmesh_internal)

Contains cmesh functionality intended for internal use within t8code only.
See [cmesh internal README](t8_cmesh_internal/README.md).

#### [src/t8_cmesh/t8_cmesh_io](t8_cmesh_io)

Contains functionality for input and output of cmeshes, in particular for reading `.msh` files (usually from the Gmsh mesh generator).
See [cmesh IO README](t8_cmesh_io/README.md).

#### [src/t8_cmesh/t8_cmesh_vertex_connectivity](t8_cmesh_vertex_connectivity)

Contains functionality for handling the vertex connectivity of a cmesh.
See [cmesh vertex connectivity README](t8_cmesh_vertex_connectivity/README.md).


An incomplete list of the source files in this folder:

#### [t8_cmesh_examples.h](t8_cmesh_examples.h)

Contains many example cases for coarse mesh generation.

TODO: Complete this list