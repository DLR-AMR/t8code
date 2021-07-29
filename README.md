### Introduction

t8code (spoken as "tetcode") is a C/C++ library to manage parallel adaptive meshes with various element types.
t8code uses a collection (a forest) of multiple connected adaptive space-trees in parallel and scales to at least one million MPI ranks and over 1 Trillion mesh elements.
It is licensed under the GNU General Public License 2.0 or later. Copyright (c) 2015 the developers.

The t8code, or T8 for short, supports the following element types (also different types in the same mesh):

0D: vertices
1D: lines
2D: quadrilaterals and triangles
3D: hexahedra, tetrahedra, prisms (pyramids currently in development).

t8code offers the following functionalities:

- Create distributed adaptive meshes
- Adapt meshes according to user given refinement/coarsening criteria
- Establish a 2:1 balance
- (Re-)partition a mesh (and associated data) among MPI ranks
- Manage ghost (halo) elements and data
- Hierarchical search in the mesh


t8code uses space-filling curves (SFCs) to manage the adaptive refinement and efficiently.
A modular approach makes it possible to exchange the underlying SFC without changing the high-level algorithms.
Thus, we can use and compare different refinement schemes and users can implement their own refinement rules if so desired.

Currently, 
  - lines use a 1D Morton curve with 1:2 refinement
  - quadrilateral/hexahedral elements are inherited from the p4est submodule, using the Morton curve 1:4, 1:8 refinement; 
  - triangular/tetrahedral are implemented using the Tetrahedral Morton curve, 1:4, 1:8 refinement;
  - prisms are implemented using the triangular TM curve and a line curve, 1:8 refinement.
  - The code supports hybrid meshes including any of the above element types (of the same dimension).

You find more information on t8code in the [t8code Wiki](https://github.com/holke/t8code/wiki).

### Setup

We provide a short guide to install t8code. 

For a more detailed description, please see the [Installation guide](https://github.com/holke/t8code/wiki/Installation) in our Wiki.

  
  To setup the project perform the following steps
  
    1.) If you cloned from github, initialize and download the git submodules
       p4est and sc.
      - git submodule init
      - git submodule update      
    2.) Call the bootstrap script in the source directory
      - ./bootstrap        
    3.) Goto your installation folder and call configure and make
      - cd /path/to/install
      - /path/to/source/configure [OPTIONS]
      - make 
      - make check
      - make install

To see a list of possible configure options, call
 
 ./configure -h

or visit [the Wiki](https://github.com/holke/t8code/wiki/Configure-Options).

Most commonly used for t8code are

  --enable-mpi    (enables MPI parallelization)
  
  --enable-debug  (enables debugging mode - massively reduces performance)
  
  --with-LIB/--without-LIB (enable/disable linking with LIB)
  
  ### Getting started
  
  To get familiar with t8code and its algorithms and data structures we recommend excuting the tutorial examples in `example/tutorials`
  and read the corresponding Wiki pages starting with [Step 0 - Helloworld](https://github.com/holke/t8code/wiki/Step-0---Hello-World).
  
