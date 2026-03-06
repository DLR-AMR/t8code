[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.7034838.svg)](https://doi.org/10.5281/zenodo.7034838)
[![t8code CI](https://github.com/DLR-AMR/t8code/actions/workflows/testsuite.yml/badge.svg)](https://github.com/DLR-AMR/t8code/actions/workflows/testsuite.yml)
[![codecov](https://codecov.io/gh/dlr-amr/t8code/branch/main/graph/badge.svg)](https://codecov.io/gh/dlr-amr/t8code)
[![docs](https://app.readthedocs.org/projects/t8code/badge/?version=latest)](https://t8code.readthedocs.io/en/latest/)

<p align="center">
  <img width="300px" src=t8code_logo.png>
</p>

### Introduction

t8code (spoken as "tetcode") is a C/C++ library to manage parallel adaptive meshes with various element types.
t8code uses a collection (a forest) of multiple connected adaptive space-trees in parallel and scales to at least one million MPI ranks and over 1 Trillion mesh elements.
It is licensed under the GNU General Public License 2.0 or later. Copyright (c) 2015 the developers.

t8code is intended to be used as a thirdparty library for numerical simulation codes or any other applications that require meshes.

<table>
    <tr>
        <td><img src="https://github.com/DLR-AMR/t8code/blob/main/doc/pictures/cmesh_tet_holes.png?raw=true" height="200" /></td> 
        <td><img src="https://github.com/DLR-AMR/t8code/blob/main/doc/pictures/flow_around_circle_sim2.jpg?raw=true" height="181" /></td>
    </tr>
      <tr>
        <td><img src="https://github.com/DLR-AMR/t8code/blob/main/doc/pictures/mesh_3d_hybrid_cutout.jpg?raw=true" height="200" /></td>
        <td><img src="https://github.com/DLR-AMR/t8code/blob/main/doc/pictures/AirplaneWithTetMesh.png?raw=true" height="200" /></td>
    </tr>
</table>

t8code, or T8 for short, supports the following element types (also different types in the same mesh):

- 0D: vertices
- 1D: lines
- 2D: quadrilaterals and triangles
- 3D: hexahedra, tetrahedra, prisms and pyramids

Among others, t8code offers the following functionalities:

- Create distributed adaptive meshes over complex domain geometries
- Adapt meshes according to user given refinement/coarsening criteria
- Establish a 2:1 balance
- (Re-)partition a mesh (and associated data) among MPI ranks
- Manage ghost (halo) elements and data
- Hierarchical search in the mesh
- Curved mesh elements

t8code uses space-filling curves (SFCs) to manage the adaptive refinement and efficiently store the mesh elements and associated data.
A modular approach makes it possible to exchange the underlying SFC without changing the high-level algorithms.
Thus, we can use and compare different refinement schemes and users can implement their own refinement rules if so desired.

Currently t8code offers the following implementations by default:
  - lines use a 1D Morton curve with 1:2 refinement
  - quadrilateral/hexahedral elements are inherited from the p4est submodule, using the Morton curve 1:4, 1:8 refinement; 
  - triangular/tetrahedral are implemented using the Tetrahedral Morton curve, 1:4, 1:8 refinement;
  - prisms are implemented using the triangular TM curve and a line curve, 1:8 refinement.
  - pyramids are implemented using the Pyramidal Morton curve and the TM curve for its tetrahedral children, 1:10 (for pyramids) / 1:8 (for tetrahedra) refinement.
  - The code supports hybrid meshes including any of the above element types (of the same dimension).

You find more information on t8code in the [t8code Wiki](https://github.com/DLR-AMR/t8code/wiki).

For a brief introduction in AMR and the algorithms used by t8code we recommend to read our [overview paper](https://elib.dlr.de/194377/1/t8code_overview_IMR2023.pdf).

### Setup

We provide a short guide to install t8code in our Wiki [Installation guide](https://github.com/DLR-AMR/t8code/wiki/Installation).

  
### Getting started
  
  To get familiar with t8code and its algorithms and data structures we recommend executing the tutorial examples in `tutorials`
  and read the corresponding Wiki pages starting with [Step 0 - Helloworld](https://github.com/DLR-AMR/t8code/wiki/Step-0---Hello-World).
  
  A sophisticated example of a complete numerical simulation is our finite volume solver of the advection equation in `example/advection`.


### Documentation

t8code uses [Doxygen](https://doxygen.nl/) to generate the code documentation. 
You can find the documentation on [readthedocs](https://t8code.readthedocs.io/en/latest/).
Follow the steps described in our Wiki [Documentation](https://github.com/DLR-AMR/t8code/wiki/Documentation) to create the documentation locally.


### License and contributing
t8code is licensed under GPLv2 (see [COPYING](COPYING)). We appreciate
contributions from the community and refer to [CONTRIBUTING.md](CONTRIBUTING.md)
for more details.

Note that we strive to be a friendly, inclusive open-source
community and ask all members of our community to adhere to our
[`CODE_OF_CONDUCT.md`](CODE_OF_CONDUCT.md). 

To get in touch, [open an issue](https://github.com/DLR-AMR/t8code/issues/new)
or write an email to one of the principal developers.

### Julia wrapper

We offer [T8code.jl](https://github.com/DLR-AMR/T8code.jl) - an official
[Julia](https://julialang.org/) package allowing to call t8code routines from
the [Julia](https://julialang.org/) programming language. From within a Julia
session do
```julia
julia> import Pkg; Pkg.add(["T8code", "MPI"])
```
to install the package on your system.

### Publications
  
  An (incomplete) list of publications related to t8code:
    
  [1] **Overview Paper**: 
  Holke, Johannes and Burstedde, Carsten and Knapp, David and Dreyer, Lukas and Elsweijer, Sandro and Ünlü, Veli and Markert, Johannes and Lilikakis, Ioannis and Böing, Niklas and Ponnusamy, Prasanna and Basermann, Achim  (2023) *t8code v. 1.0 - Modular Adaptive Mesh Refinement in the Exascale Era*. SIAM International Meshing Round Table 2023, 06.03.2023 - 09.03.2023, Amsterdam, Niederlande. 
  [Full text available](https://elib.dlr.de/194377/1/t8code_overview_IMR2023.pdf)
  
    
  [2] **Original PhD thesis**: 
  Holke, Johannes *Scalable algorithms for parallel tree-based adaptive mesh refinement with general element types*, PhD thesis at University of Bonn, 2018,
      [Full text available](https://bonndoc.ulb.uni-bonn.de/xmlui/handle/20.500.11811/7661)
   
      
  [3] **Tetrahedral and triangular Space-filling curve**:
  Burstedde, Carsten and Holke, Johannes *A Tetrahedral Space-Filling Curve for Nonconforming Adaptive Meshes*, SIAM Journal on Scientific Computing, 2016, [10.1137/15M1040049](https://epubs.siam.org/doi/10.1137/15M1040049)
  
  
  [4] **Coarse mesh partitioning**:
  Burstedde, Carsten and Holke, Johannes *Coarse mesh partitioning for tree-based AMR*, SIAM Journal on Scientific Computing, 2017, [10.1137/16M1103518](https://epubs.siam.org/doi/10.1137/16M1103518)
  
  
  [5] **Ghost computation**:
  Holke, Johannes and Knapp, David and Burstedde, Carsten *An Optimized, Parallel Computation of the Ghost Layer for Adaptive Hybrid Forest Meshes*, SIAM Journal on Scientific Computing, 2021, [10.1137/20M1383033](https://epubs.siam.org/doi/abs/10.1137/20M1383033)
 
  
  [6] **Geometry controlled refinement for hexahedra**:
  Elsweijer, Sandro and Holke, Johannes and Kleinert, Jan and Reith, Dirk  (2022) *Constructing a Volume Geometry Map for Hexahedra with Curved Boundary Geometries*.   In: SIAM International Meshing Roundtable Workshop 2022.  SIAM International Meshing Roundtable Workshop 2022, 22. - 25. Feb. 2022, [Full text available](https://elib.dlr.de/186570/1/ConstructingAVolumeGeometryMapForHexahedraWithCurvedBoundaryGeometries.pdf) 
  
  [7] **JOSS entry**:
  Holke, Johannes and Markert, Johannes, et. al. (2025) *t8code - modular adaptive mesh refinement in the exascale era*. In: Journal of Open Source Software, [Full text available](https://www.theoj.org/joss-papers/joss.06887/10.21105.joss.06887.pdf)
 
### Theses with t8code relations

  An (incomplete) list of theses written with or about t8code:
  

  [A] **Prism space-filling curve**: 
  Knapp, David (2017) *Adaptive Verfeinerung von Prismen*. Bachelor's thesis, Rheinische Friedrich-Wilhems-Universität Bonn.
  
  
  [B] **Pyramidal space-filling curve**: 
  Knapp, David (2020) *A space-filling curve for pyramidal adaptive mesh refinement*. Master's thesis, Rheinische Friedrich-Wilhems-Universität Bonn. [Full text available](https://www.researchgate.net/publication/346789160_A_space-filling_curve_for_pyramidal_adaptive_mesh_refinement)
  
   
  [C] **DG solver based on t8code**: 
  Dreyer, Lukas (2021) *The local discontinuous galerkin method for the advection-diffusion equation on adaptive meshes*.  Master's thesis, Rheinische Friedrich-Wilhems-Universität Bonn.
  [Full text available](https://elib.dlr.de/143969/1/masterthesis_dreyer.pdf) 
  
  
  [D] **Geometry controlled refinement for hexahedra (Part 1)**: 
  Elsweijer, Sandro (2021) *Curved Domain Adaptive Mesh Refinement with Hexahedra*.  Tech report, Hochschule Bonn-Rhein-Sieg.
  [Full text available](https://elib.dlr.de/186571/1/masterprojekt-2_elsweijer_ABGABEVERSION_TITEL.pdf)
  

  [E] **Subelement and resolving hanging faces in 2D**: 
  Becker, Florian (2021) *Removing hanging faces from tree-based adaptive meshes for numerical simulation*, Master's thesis, Universität zu Köln.
  [Full text available](https://elib.dlr.de/187499/1/RemovingHangingFacesFromTreeBasedAMR.pdf)
  
  
  [F] **Coarsening as post-processing to reduce simulation file size**: 
  Spataro, Luca  (2021) *Lossy data compression for atmospheric chemistry using adaptive mesh coarsening*.  Master's thesis, Technische Universität München.
  [Full text available](https://elib.dlr.de/144997/1/master-thesis-final-spataro.pdf)
 
  
  [G] **Geometry controlled refinement for hexahedra (Part 2)**: 
  Elsweijer, Sandro (2022) *Evaluation and generic application scenarios for curved hexahedral adaptive mesh refinement*.  Master's thesis, Hochschule Bonn-Rhein-Sieg.  [10.13140/RG.2.2.34714.11203](<https://doi.org/10.13140/RG.2.2.34714.11203>) [Full text available](https://elib.dlr.de/186561/1/sandro_elsweijer-evaluation_and_generic_application_scenarios_for_curved_hexahedral_adaptive_mesh_refinement.pdf)
  
  
  [H] **Multigrid and other preconditioners for DG**: 
  Böing, Niklas  (2022) *Evaluation of preconditioners for implicit solvers of local DG for the advection-diffusion equation* (*Untersuchung von Präkonditionierern für implizite Löser für das Local DG-Verfahren zur Lösung der Advektions-Diffusionsgleichung*).  Master's thesis, Universität zu Köln.
[Full text available](https://elib.dlr.de/186347/1/Untersuchung%20von%20Pr%C3%A4konditionierern%20f%C3%BCr%20implizite%20L%C3%B6ser%20f%C3%BCr%20das%20Local%20DG-Verfahren%20zur%20L%C3%B6sung%20der%20Advektions-Diffusionsgleichung.pdf)
  

  [I] **Removing elements from the mesh (cutting holes)**: 
  Lilikakis, Ioannis  (2022) *Algorithms for tree-based adaptive meshes with incomplete trees*.  Master's thesis, Universität zu Köln.    
 [Full text may be available in future](https://elib.dlr.de/191968/)
 
  [J] **Curved tetrahedra**:
  Fussbroich Jakob (2023) *Towards high-order, hybrid adaptive mesh refinement: Implementation and evaluation of curved unstructured mesh elements*. Master's thesis. Technische Hochschule Köln.
  [Full text available](https://elib.dlr.de/200442/)

  [K] **Hanging node resolution 3D**:
  Tabea Leistikow (2024) *Derivation and implementation of a hanging nodes resolution scheme for hexahedral non-conforming meshes in t8code*. Master's thesis, Universität zu Köln.
  Full text currently not available.

  ### Citing t8code
  
  If you use t8code in any of your publications, please cite the [github repository](https://doi.org/10.5281/zenodo.7034838), [1] and [2]. For publications specifically related to 
- **the tetrahedral index**, please cite [3].
- **coarse mesh partitioning**, please cite [4].
- **construction and handling of the ghost layer**, please cite [5].
- **geometry controlled refinement**, please cite [6] (general) and [J] (tetrahedral).
- **hanging node resolution and/or subelements**, please cite [E] and [K].

If you use any functionality described in the theses, we encourage you to cite them as well.
