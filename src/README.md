# t8code/src

This folder contains the main source code of t8code.

Here is an overview over the different modules:

### [src](.)

The main source folder. It holds the main header file [t8.h](t8.h) and the following sub-folders:

#### [src/t8_cad](t8_cad)

Contains the functionality of handling Open CASCADE shapes.
See [cad README](t8_cad/README.md).

#### [src/t8_cmesh](t8_cmesh)

Contains implementation details of algorithms and data structures related to the coarse mesh.
See [cmesh README](t8_cmesh/README.md) and the tutorials about coarse meshes: [Creating a coarse mesh](https://github.com/DLR-AMR/t8code/wiki/Step-1---Creating-a-coarse-mesh)
and [Building a Cmesh by hand](https://github.com/DLR-AMR/t8code/wiki/Building-a-Cmesh-by-hand).

#### [src/t8_data](t8_data)

Contains data handling classes and algorithms related to data containers, including arrays to handle element data.
See [data README](t8_data/README.md).

#### [src/t8_eclass](t8_eclass)

Contains the definitions of all possible element classes supported by t8code.
See [eclass README](t8_eclass/README.md).

#### [src/t8_element](t8_element)

Contains the definitions of the opaque element structure that acts like an abstract base class for all elements, along with some associated constants.
See [element README](t8_element/README.md).

#### [src/t8_forest](t8_forest)

Contains implementation details of algorithms and data structures related to the forest, i.e. the actual computational mesh.
See [forest README](t8_forest/README.md) and the tutorials, [Step 2](https://github.com/DLR-AMR/t8code/wiki/Step-2---Creating-a-uniform-forest), [Step 3](https://github.com/DLR-AMR/t8code/wiki/Step-3---Adapting-a-forest), [Step 4](https://github.com/DLR-AMR/t8code/wiki/Step-4---Partition,-Balance,-Ghost).

#### [src/t8_geometry](t8_geometry)

Contains data handling classes and algorithms for geometry handling, i.e. mapping the elements into a physical domain space.
See [geometry README](t8_geometry/README.md) and the tutorial [Curved Meshes](https://github.com/DLR-AMR/t8code/wiki/Feature---Curved-meshes) (note that the geometry also handles linear, non-curved meshes).

#### [src/t8_helper_functions](t8_helper_functions)

Contains some helper functions used by t8code internally, e.g., specialized algorithms on std::vector, such as `vector_split`.
See [helper functions README](t8_helper_functions/README.md).

#### [src/t8_misc](t8_misc)

Contains some files that did not fit into any other category, such as the version handling or some Windows-specific setups. 
See [miscellaneous README](t8_misc/README.md).

#### [src/t8_schemes](t8_schemes)

Contains data handling classes and algorithms for schemes. Schemes are the implementation of the space-filling curve for each element shape.
A scheme describes how elements refine, construct their neighbors, etc.
See [scheme README](t8_schemes/README.md).

#### [src/t8_types](t8_types)

Custom data types and strong types. Also contains `t8_vec` vector implementation.
See [types README](t8_types/README.md).


#### [src/t8_vtk](t8_vtk)

Contains algorithms and data structures for the VTK interface.
See [vtk README](t8_vtk/README.md).



