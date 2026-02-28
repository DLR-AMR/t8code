This folder contains the coarse mesh interface.

The coarse mesh (t8_cmesh) contains information about the domain topology and geometry.
It is usually created from an input mesh, either from a mesh generator or by hand.

See also our Wiki tutorials [Creating a coarse mesh](https://github.com/DLR-AMR/t8code/wiki/Step-1---Creating-a-coarse-mesh)
and [Building a Cmesh by hand](https://github.com/DLR-AMR/t8code/wiki/Building-a-Cmesh-by-hand).

Your entry point should be [t8_cmesh.h](https://github.com/DLR-AMR/t8code/blob/18e8af64286bef9370042a8c4d7885d279d2b933/src/t8_cmesh.h) which is the main cmesh interface.

An incomplete list of the source files in this folder:



#### [t8_cmesh_examples.h](https://github.com/DLR-AMR/t8code/blob/feature-folder_README/src/t8_cmesh/t8_cmesh_examples.h)

Contains many example cases for coarse mesh generation.

### Input/Output

#### [t8_cmesh_readmshfile.h](https://github.com/DLR-AMR/t8code/blob/feature-folder_README/src/t8_cmesh_readmshfile.h)

Input routines for reading `.msh` files (usually from the Gmsh mesh generator).

#### [t8_cmesh_tetgen.h](https://github.com/DLR-AMR/t8code/blob/18e8af64286bef9370042a8c4d7885d279d2b933/src/t8_cmesh_tetgen.h) and [t8_cmesh_triangle.h](https://github.com/DLR-AMR/t8code/blob/18e8af64286bef9370042a8c4d7885d279d2b933/src/t8_cmesh_triangle.h)

Input routines for the TRIANGLE and TETGEN mesh generators.
These have not been maintained for a long time and might be incorrect.

#### [t8_cmesh_netcdf.h](https://github.com/DLR-AMR/t8code/blob/18e8af64286bef9370042a8c4d7885d279d2b933/src/t8_cmesh_netcdf.h)

netcdf output

### Internal

#### [t8_cmesh_types.h](https://github.com/DLR-AMR/t8code/blob/18e8af64286bef9370042a8c4d7885d279d2b933/src/t8_cmesh/t8_cmesh_types.h)

Contains the class definition of `t8_cmesh` and various subclasses.
This is internal and should only be touched if you need to add member variables to the cmesh.

#### [t8_cmesh_stash.h](https://github.com/DLR-AMR/t8code/blob/feature-folder_README/src/t8_cmesh/t8_cmesh_stash.h)

The stash is an intermediate data structure that we use while a cmesh is being constructed.
It is very unhandy and contains lots of void pointers.
Do not touch unless really necessary.
