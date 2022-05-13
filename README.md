[![t8code tests](https://github.com/holke/t8code/actions/workflows/tests.yml/badge.svg)](https://github.com/holke/t8code/actions/workflows/tests.yml)

### Introduction

t8code (spoken as "tetcode") is a C/C++ library to manage parallel adaptive meshes with various element types.
t8code uses a collection (a forest) of multiple connected adaptive space-trees in parallel and scales to at least one million MPI ranks and over 1 Trillion mesh elements.
It is licensed under the GNU General Public License 2.0 or later. Copyright (c) 2015 the developers.

t8code is intended to be used as a thirdparty library for numerical simulation codes or any other applications that require meshes.

<table>
    <tr>
        <td><img src="https://github.com/holke/t8code/blob/main/doc/pictures/cmesh_tet_holes.png?raw=true" height="200" /></td> 
        <td><img src="https://github.com/holke/t8code/blob/main/doc/pictures/flow_around_circle_sim2.jpg?raw=true" height="181" /></td>
    </tr>
      <tr>
        <td><img src="https://github.com/holke/t8code/blob/main/doc/pictures/mesh_3d_hybrid_cutout.jpg?raw=true" height="200" /></td>
        <td><img src="https://github.com/holke/t8code/blob/main/doc/pictures/AirplaneWithTetMesh.png?raw=true" height="200" /></td>
    </tr>
</table>

t8code, or T8 for short, supports the following element types (also different types in the same mesh):

- 0D: vertices
- 1D: lines
- 2D: quadrilaterals and triangles
- 3D: hexahedra, tetrahedra, prisms (pyramids currently in development).

Among others, t8code offers the following functionalities:

- Create distributed adaptive meshes over complex domain geometries
- Adapt meshes according to user given refinement/coarsening criteria
- Establish a 2:1 balance
- (Re-)partition a mesh (and associated data) among MPI ranks
- Manage ghost (halo) elements and data
- Hierarchical search in the mesh


t8code uses space-filling curves (SFCs) to manage the adaptive refinement and efficiently store the mesh elements and associated data.
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

#### Requirements

- [libsc](https://github.com/cburstedde/libsc) (Included in t8code's git repository)
- [p4est](https://github.com/cburstedde/p4est) (Included in t8code's git repository)
- automake
- libtool
- make

Optional
- The VTK library for advanced VTK output (basic VTK output is provided without linking against VTK)
- The netcdf library for netcdf file output

  
#### Steps
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
  
#### Example configurations

For a parallel release mode with local installation path `$HOME/t8code_install`:

`configure --enable-mpi CFLAGS=-O3 CXXFLAGS=-O3 --prefix=$HOME/t8code_install`

For a debugging mode with static linkage (makes using gdb and valgrind more comfortable):

`configure --enable-mpi --enable-debug --enable-static --disable-shared CFLAGS="-Wall -O0 -g" CXXFLAGS="-Wall -O0 -g"`
  
### Getting started
  
  To get familiar with t8code and its algorithms and data structures we recommend executing the tutorial examples in `example/tutorials`
  and read the corresponding Wiki pages starting with [Step 0 - Helloworld](https://github.com/holke/t8code/wiki/Step-0---Hello-World).
  
  A sophisticated example of a complete numerical simulation is our finite volume solver of the advection equation in `example/advection`.
  
  ### Publications
  
  An (incomplete) list of publications related to t8code:
    
  [1] Johannes Holke, Scalable algorithms for parallel tree-based adaptive mesh refinement with general element types, PhD thesis at University of Bonn, 2018,
      [Full text available](https://bonndoc.ulb.uni-bonn.de/xmlui/handle/20.500.11811/7661)
      
  [2] Carsten Burstedde and Johannes Holke, A Tetrahedral Space-Filling Curve for Nonconforming Adaptive Meshes, SIAM Journal on Scientific Computing, 2016, 10.1137/15M1040049
  
  [3] Carsten Burstedde and Johannes Holke, Coarse mesh partitioning for tree-based AMR, SIAM Journal on Scientific Computing, 2017, 10.1137/16M1103518
  
  [4] Johannes Holke and David Knapp and Carsten Burstedde, An Optimized, Parallel Computation of the Ghost Layer for Adaptive Hybrid Forest Meshes, Submitted to SIAM Journal on Scientific Computing, [Preprint available](https://arxiv.org/abs/1910.10641) 2019
  
  ### Citing t8code
  
  If you use t8code in any of your publications, please cite the [github repository](https://github.com/holke/t8code) and [1]. For publications specifically related to 
- the TM index, please cite [2].
- coarse mesh partitioning, please cite [3].
- construction and handling of the ghost layer, please cite [4].
