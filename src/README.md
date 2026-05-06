# t8code/src

This folder contains the main source code of t8code.

Here is an overview over the different modules:

### [src](.)

The main source folder. It holds main header files like [t8.h](t8.h), [t8_cmesh.hxx](t8_cmesh.hxx) and more.

#### [src/t8_cmesh](t8_cmesh)

Contains implementation details of algorithms and data structures related to the coarse mesh.
See [cmesh README](t8_cmesh/README.md) and the tutorials about coarse meshes: [Creating a coarse mesh](https://github.com/DLR-AMR/t8code/wiki/Step-1---Creating-a-coarse-mesh)
and [Building a Cmesh by hand](https://github.com/DLR-AMR/t8code/wiki/Building-a-Cmesh-by-hand).

#### [src/t8_data](t8_data)

Contains data handling classes and algorithms related to data containers, including arrays to handle element data.
See [data README](t8_data/README).

#### [src/t8_forest](t8_forest)

Contains implementation details of algorithms and data structures related to the forest, i.e. the actual computational mesh.
See [forest README](t8_forest/README) and the tutorials, [Step 2](https://github.com/DLR-AMR/t8code/wiki/Step-2---Creating-a-uniform-forest), [Step 3](https://github.com/DLR-AMR/t8code/wiki/Step-3---Adapting-a-forest), [Step 4](https://github.com/DLR-AMR/t8code/wiki/Step-4---Partition,-Balance,-Ghost).

#### [src/t8_geometry](t8_geometry)

Contains data handling classes and algorithms for geometry handling, i.e. mapping the elements into a physical domain space.
See [geometry README](t8_geometry/README) and the tutorial [Curved Meshes](https://github.com/DLR-AMR/t8code/wiki/Feature---Curved-meshes) (note that the geometry also handles linear, non-curved meshes).

#### [src/t8_schemes](t8_schemes)

Contains data handling classes and algorithms for schemes. Schemes are the implementation of the space-filling curve for each element shape.
A scheme describes how elements refine, construct their neighbors, etc.
See [scheme README](t8_schemes/README).

#### [src/t8_types](t8_types)

Custom data types and strong types. Also contains `t8_vec` vector implementation.
See [types README](t8_types/README).

#### [src/t8_vector_helper](t8_vector_helper)

Contains specialized algorithms on std::vector in particular `vector_split`.
See [vector helper README](t8_vector_helper/README).

#### [src/t8_vtk](t8_vtk)

Contains algorithms and data structures for the VTK interface.
See [vtk README](t8_vtk/README).



