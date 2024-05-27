.. t8code documentation master file, created by
   sphinx-quickstart on Wed May 22 13:19:02 2024.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Welcome to t8code's documentation!
==================================
This is the documentation for t8code version |release|

Author
======
    Johannes Holke, David Knapp, Lukas Dreyer, Niklas Böing, Sandro Elsweijer, Ioannis Lilikakis, Veli Ünlü, Carsten Burstedde, Johannes Markert et. al. 

Copyright
=========
    GNU General Public License version 2 (or, at your option, any later version)

t8code (spoken as "tetcode") is a C/C++ library to manage parallel adaptive meshes with various element types. t8code uses a collection (a forest) of multiple connected adaptive space-trees in parallel and scales to at least one million MPI ranks and over 1 Trillion mesh elements.

t8code is intended to be used as a thirdparty library for numerical simulation codes or any other applications that require meshes.

The best place to start learning about t8code is the Wiki, especially the Installation Guide and the Tutorials.

You find the source code to t8code on github.

t8code supports the following element types (also different types in the same mesh):

   - 0D: vertices
   - 1D: lines
   - 2D: quadrilaterals and triangles
   - 3D: hexahedra, tetrahedra, prisms, and pyramids (currently in development).

Among others, t8code offers the following functionalities:

   - Create distributed adaptive meshes over complex domain geometries
   - Adapt meshes according to user given refinement/coarsening criteria
   - Establish a 2:1 balance
   - (Re-)partition a mesh (and associated data) among MPI ranks
   - Manage ghost (halo) elements and data
   - Hierarchical search in the mesh

t8code uses space-filling curves (SFCs) to manage the adaptive refinement and efficiently store the mesh elements and associated data. A modular approach makes it possible to exchange the underlying SFC without changing the high-level algorithms. Thus, we can use and compare different refinement schemes and users can implement their own refinement rules if so desired.

Currently,

   - lines use a 1D Morton curve with 1:2 refinement
   - quadrilaterals/hexahedra are inherited from the p4est submodule, using the Morton curve 1:4, 1:8 refinement;
   - triangles/tetrahedra are implemented using the Tetrahedral Morton curve, 1:4, 1:8 refinement;
   - prisms are implemented using the triangular TM curve and a line curve, 1:8 refinement.
   - pyramids are implemented via a pyramidal Morton code.
   - The code supports hybrid meshes including any of the above element types (of the same dimension).

You find more information on t8code in the t8code Wiki.


.. toctree::
   :maxdepth: 2
   :caption: Contents:



Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`

Docs
====

.. doxygenindex::
   :allow-dot-graphs:
.. doxygenpage:: dotgraphs