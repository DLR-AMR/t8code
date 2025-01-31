---
title: 't8code - modular adaptive mesh refinement in the exascale era'
tags:
  - C
  - C++
  - adaptive mesh refinement
  - exascale
  - hybrid meshes
  - modularity
  - high performance computing
authors:
  - name: Johannes Holke
    orcid: 0000-0002-2783-3286
    affiliation: 1
    equal-contrib: true
  - name: Johannes Markert
    orcid: 0000-0001-6297-9494
    corresponding: true
    equal-contrib: true
    affiliation: 1
  - name: David Knapp
    orcid: 0000-0002-6305-1572
    equal-contrib: true
    affiliation: 1
  - name: Lukas Dreyer
    orcid: 0000-0001-7484-3674
    equal-contrib: true
    affiliation: 1
  - name: Sandro Elsweijer
    orcid: 0000-0002-2753-3873
    equal-contrib: true
    affiliation: 1
  - name: Niklas Böing
    equal-contrib: true
    affiliation: 1
  - name: Ioannis Lilikakis
    equal-contrib: true
    affiliation: 1
  - name: Jakob Fussbroich
    orcid: 0000-0003-0784-2182
    equal-contrib: true
    affiliation: 1
  - name: Tabea Leistikow
    equal-contrib: true
    affiliation: 1
  - name: Florian Becker
    orcid: 0000-0002-8384-9282
    equal-contrib: true
    affiliation: 1
    orcid: 0000-0002-8384-9282
  - name: Uenlue, Veli
    equal-contrib: true
    affiliation: 1
  - name: Albers, Ole
    orcid: 0000-0002-4950-7297
    equal-contrib: true
    affiliation: 1
  - name: Carsten Burstedde
    orcid: 0000-0001-9843-1041
    affiliation: 2
  - name: Achim Basermann
    orcid: 0000-0003-3637-3231
    affiliation: 1
  - name: Hergl, Chiara
    orcid: 0000-0002-4016-9113
    affiliation: 1
  - name: Julia, Weber
    affiliation: 1
  - name: Schoenlein, Kathrin
    affiliation: 1
  - name: Ackerschott, Jonas
    affiliation: 10
  - name: Evgenii, Andreev
    affiliation: 10
  - name: Csati, Zoltan
    affiliation: 4
  - name: Dutka, Alexandra
    affiliation: 4
  - name: Geihe, Benedict
    affiliation: 5
  - name: Kestener, Pierre
    affiliation: 6
  - name: Kirby, Andrew
    affiliation: 7
  - name: Ranocha, Hendrik
    affiliation: 8
  - name: Schlottke-Lakemper, Michael
    affiliation: 9

affiliations:
 - name: German Aerospace Center (DLR), Institute of Software Technology, Department High-Performance Computing, Cologne, Germany
   index: 1
 - name: Rheinische Friedrich-Wilhelms-Universität Bonn, Institute for Numerical Simulations and Hausdorff Center for Mathematics,  Germany
   index: 2
 - name: Forschungszentrum Juelich (JSC)
   index: 3
 - name: CERFACS
   index: 4
 - name: University of Cologne
   index: 5
 - name: CEA
   index: 6
 - name: University of Wyoming
   index: 7
 - name: University of Hamburg
   index: 8
 - name: RWTH Aachen
   index: 9
 - name: unaffiliated
   index: 10

date: 29 October 2024
bibliography: paper.bib

# Optional fields if submitting to a AAS journal too, see this blog post:
# https://blog.joss.theoj.org/2018/12/a-new-collaboration-with-aas-publishing
# aas-doi: 10.3847/xxxxx <- update this with the DOI from AAS once you know it.
# aas-journal: Astrophysical Journal <- The name of the AAS journal.
---

# Summary

In this paper, we present our scalable dynamic adaptive mesh refinement (AMR)
library `t8code`, which was officially released in 2022 [@Holke_t8code_2022].
`t8code` is written in C/C++, open source, and readily available at
[dlr-amr.github.io/t8code](https://dlr-amr.github.io/t8code/). It is developed
and maintained at the [Institute of Software
Technology](https://www.dlr.de/sc/en/) of the German Aerospace Center (DLR).
AMR is a widely used method of locally adapting the mesh resolution according
to an adequate error indicator in grid-based applications - especially in the
context of computational fluid dynamics. Our software library provides fast and
memory efficient parallel algorithms for dynamic AMR to handle tasks such as
mesh adaptation, load-balancing, ghost computation, feature search and more.
`t8code` can manage meshes with over one trillion mesh elements
[@holke_optimized_2021] and scales up to one million parallel processes
[@holke_scalable_2018]. It is intended to be used as mesh management backend in
scientific and engineering simulation codes paving the way towards
high-performance applications of the upcoming exascale era.

# Statement of Need

Adaptive Mesh Refinement has been established as a successful approach for
scientific and engineering simulations over the past decades
[@TEUNISSEN2019106866; @DEAL2; @DOERFLER1996;
@BABUVSKA1978]. By modifying the mesh resolution locally according to
problem specific indicators, the computational power is efficiently
concentrated where needed and the overall memory usage is reduced by orders of
magnitude. However, managing adaptive meshes and associated data is a very
challenging task, especially for parallel codes. Implementing fast and scalable
AMR routines generally leads to a large development overhead motivating the
need for external mesh management libraries like `t8code`. Our target
audiences are scientists and application developers working on grid-based
simulation and visualization frameworks who are looking for a comprehensive and
versatile mesh management solution. Besides offering AMR we also aim to lower
the threshold to parallelize their codes by solely interacting with t8code's
API. Alternative AMR libraries with a similar range of features are p4est
[@BursteddeWilcoxGhattas11], libMesh [@libMeshPaper], PARAMESH [@macneice2000paramesh], and
SAMRAI [@gunney2013scalable]. 

In contrast to the other AMR solutions, only `t8code` natively supports
recursive refinement on a wide range of element types: vertices, lines,
quadrilaterals, triangles, hexahedra, tetrahedra, prisms, and pyramids.
Additionally, extensions to other refinement patterns and element shapes are
straightforwardly supported due to `t8code`'s modular code structure and clear
distinction between low- and high-level mesh operations.  This gives our AMR
solution an unique position in the market catering for a wide range of use
cases. Currently, `t8code` is optimized for grid-based applications using
face-to-face connectivity between elements, such as Finite-Volume and
Discontinuous Galerkin methods. In the future, we plan to support node-to-node
connectivity and hanging nodes resolution to further increase the range of
applications, such as Finite Element methods.

# Exemplary application

\autoref{fig:visploremesh} depicts an examplary adapted mesh managed by
`t8code` using two different element types: quads and triangles. Shown is the
temperature profile of a convection simulation of a model planet's mantle
(source: Institute of Planetary Research, DLR). The original, uniform mesh
consists of over 158 million quad cells allocating 6.818 GB of memory.  By
applying AMR to the data the memory usage could be reduced down to 20\% with
a compression error of less than 1\%. The error measure was chosen to be the
norm of the variance between refinement resp. coarsening steps. That is,
starting from the uniform mesh at highest refinement level ($l = 8$), the mesh
was successively coarsened until the disagreement from the original data reached
1\%. It should be noted that `t8code`'s primary objective is to provide
flexible adaptive mesh management. The layout of the data inside an element and
its interpretation regarding, for example, when and how to refine/coarsen is up
to the application.

![Visualization of a planetary mantle convection simulation (source: Institute
of Planetary Research, DLR). Shown is the 2D slice of the temperature profile.
Left: original uniform data. The highlighting of the grid lines was omitted for
visual clarity. Middle: adapted mesh with quad elements. Right: adapted mesh
with triangle elements. The original data living on a uniform quad mesh was
first transferred to a triangle mesh and adapted afterwards. This shows the
versatility of t8code regarding to the choice of mesh elements. \label{fig:visploremesh}](pics/Gaia_original_vs_AMR.png)

# Fundamental Concepts

`t8code` is based on the forest-of-trees approach. The starting point
for usage of `t8code` is an unstructured conformal input mesh, which
we denote a coarse mesh. This coarse mesh describes the geometry of the
computational domain and is usually provided by a mesh generator such as 
Gmsh [@geuzaine2009gmsh]. Each of the coarse mesh cells is then viewed as the
root of a refinement tree. These trees are refined recursively in a structured
pattern, resulting in a collection of trees, which we call a forest. `t8code`
stores only a minimal amount of information about the finest elements of the mesh -
the leaves of the trees - in order to reconstruct the whole forest.

By enumerating the leaves in a recursive refinement pattern we obtain a
space-filling curve (SFC) logic. Via these SFCs, all elements in a refinement
tree are assigned an integer-based index and are stored in linear order.
Element coordinates or element neighbors do not need to be stored explicitly
but can be reconstructed from the SFC index. Fast bitwise SFC operations ensure
optimal runtimes and diminish the need for memory lookups. Moreover, the SFC is
used to distribute the forest mesh across multiple processes, so that each
process only stores a unique portion of the SFC. See
\autoref{fig:SpaceFillingCurves} for an illustration of the concept.

While being successfully applied to quadrilateral
and hexahedral meshes [@burstedde_p4est_2011; @weinzierl_peano_2019],
these SFC techniques are extended by `t8code` in a modular fashion, such that arbitrary
element shapes are supported. We achieve this modularity through a novel
decoupling approach that separates high-level (mesh global) algorithms from
low-level (element local) implementations. All high-level algorithms can
be applied to different implementations of element shapes and refinement
patterns. A mix of different element shapes in the same mesh is also
supported.

Mesh adaptation as it is done in t8code leads to hanging nodes. Numerical
methods have to specifically handle these non-conforming interfaces.
Finite-Volume schemes or Discontinuous Galerkin methods naturally treat this
problem via so-called mortar methods. In the future, it is planned to also
support hanging nodes resolving routines by inserting transition elements
conformally connecting elements at different refinement levels.

![Left: Exemplary t8code forest mesh consisting of two trees (k0, k1)
distributed over three parallel processes p0 to p2. The SFC is represented by a
black curve tracing only the finest elements (leaves) of each tree. Right:
Sketch of the associated mixed shape (a triangle and a quad) mesh refined up to
level three. \label{fig:SpaceFillingCurves}](pics/t8code_sfc_hybrid_tree_vs_mesh.png)

# Performance

`t8code` supports distributed coarse meshes of arbitrary size and complexity,
which we tested for up to 370 million coarse mesh cells
[@burstedde_coarse_2017].  Moreover, we conducted various performance studies
on the JUQUEEN and the JUWELS supercomputers at the Jülich Supercomputing
Center. In \autoref{tab:t8code_runtimes}, [@holke_optimized_2021] we show that
`t8code`'s ghost routine is exceptionally fast with proper scaling of up to 1.1
trillion mesh elements. Computing ghost layers around parallel domains is
usually the most expensive of all mesh operations. To put these results into
perspective, we conducted scaling tests on the terrabyte cluster at Leibniz
Supercomputing Centre comparing the ghost layer creation runtimes of p4est and
t8code. In \autoref{fig:ghost_layer_runtimes} the measured runtimes of both
libraries are plotted over the number of processes. The p4est library has been
established as one of the most performant meshing libraries
[@BursteddeWilcoxGhattas11] specializing on adaptive quadrilateral and
hexahedral meshes. Clearly, t8code shows near perfect scaling for tetrahedral
meshes on par with p4est. The absolute runtime of t8code is around 1.5 times
the runtime of p4est measured on a per ghost element basis. This is expected
since the ghost layer algorithm is more complex and thus a bit less optimized,
while supporting a wider range of element types.

Furthermore, in a prototype code [@Dreyer2021] implementing a high-order
Discontinuous Galerkin (DG) method for advection-diffusion equations on
dynamically adaptive hexahedral meshes we can report of a 12 times speed-up
compared to non-AMR meshes with only an overall 15\% runtime contribution of
`t8code`. In \autoref{fig:t8code_runtimes} we compare the runtimes over number
of processes of the DG solver and the summed mesh operations done by t8code
which are ghost computation, ghost data exchange, partitioning (load
balancing), refinement, and coarsening as well as balancing ensuring only a
difference of one refinement level among element's face neighbors. From the
graphs in \autoref{fig:t8code_runtimes} we clearly see that `t8code` only takes
around 15\% to 20\% of overall runtime compared to the solver.

+----------------+-------------------+--------------------+--------+
| \# Process     | \# Elements       | \# Elem. / process |  Ghost |
+:==============:+:=================:+:==================:+:======:+
|     49,152     | 1,099,511,627,776 |     22,369,621     | 2.08 s |
+----------------+-------------------+--------------------+--------+
|     98,304     | 1,099,511,627,776 |     11,184,811     | 1.43 s |
+================+===================+====================+========+
| Table 1: Runtimes on JUQUEEN for the ghost layer                 |
| computation for a distributed mesh consisting of 1.1 trillion    |
| elements. \label{tab:t8code_runtimes}                            |
+================+===================+====================+========+

![Runtimes of ghost layer creation on the terrabyte cluster for p4est and
t8code. The meshes have been refined into a Menger sponge for hexahedral mesh
with p4est (max. level 12) and a Sierpinski sponge for the tetrahedral mesh in
t8code (max. level 13) to create a fractal pattern with billions of elements as
a stress test. To make the two runs comparable the runtimes have been divided
by the average local number of ghost elements on a MPI rank.
\label{fig:ghost_layer_runtimes}
](pics/plot-timings-per-num-ghosts.png){width="90%"}

![Runtimes on JUQUEEN of the solver and summed mesh operations of our DG
prototype code coupled with `t8code`. Mesh operations are ghost computation,
ghost data exchange, partitioning (load balancing), refinement and coarsening
as well as balancing (max. difference of one level of refinement of neighboring
elements). t8code only takes around 15\% to 20\% of the overall runtime.
\label{fig:t8code_runtimes}
](pics/t8code-runtimes-simple.png){width="90%"}

# Research Projects

Even though `t8code` is a newcomer to the market, it is already in use as the
mesh management backend in various research projects, most notably in the earth
system modeling (ESM) community. In the
[ADAPTEX](https://dlr-amr.github.io/adaptex/) project `t8code` is integrated
with the [Trixi framework](https://trixi-framework.github.io/)
[@schlottkelakemper2020trixi] - a modern computational fluid dynamics code
written in [Julia](https://julialang.org/). Over the next years several ESM
applications are planned to couple to this combination, including
[MESSy](https://messy-interface.org),
[MPTrac](https://helmholtz.software/software/mptrac), and
[SERGHEI](https://helmholtz.software/software/serghei).  Moreover, `t8code`
also plays an important role in several DLR funded research projects, e.g.,
[VisPlore](https://www.dlr.de/en/research-and-transfer/projects-and-missions/visplore)
(massive data visualization), [HYTAZER](https://elib.dlr.de/201347/) (hydrogen
tank certification), and
[Greenstars](https://www.dlr.de/en/ra/research-transfer/projects/dlr-projects/green-propulsion-free-flight-demonstrator-the-lander)
(additive rocket engine manufacturing).

# Further Information

For further information beyond this short note and also for code examples, we
refer to our
[Documentation](https://dlr-amr.github.io/t8code/pages/documentation.html) and
[Wiki](https://github.com/DLR-AMR/t8code/wiki) reachable via our homepage
[dlr-amr.github.io/t8code](https://dlr-amr.github.io/t8code/) and our technical
publications on `t8code` [@holke_scalable_2018; @burstedde_coarse_2017;
@holke_optimized_2021; @burstedde_tetrahedral_2016; @Knapp20;
@Becker_hanging_faces; @elsweijer_curved_2021; @elsweijer_evaluation_2022; @Dreyer2021;
@Lilikakis_removing; @Holke_t8code_2022; @Fussbroich_towards_2023].

# Acknowledgements

Johannes Holke thanks the Bonn International School Graduate School of
Mathematics (BIGS) for funding the initial development of `t8code`.  Further
development work was funded by the German Research Foundation as part of
project 467255783, the European Union via NextGenerationEU and the German
Federal Ministry of Research and Education (BMBF) as part of the ADAPTEX and
PADME-AM projects. Development work was performed as part of the Helmholtz School
for Data Science in Life, Earth and Energy (HDS-LEE) and received funding from
the Helmholtz Association of German Research Centres. The development team of
`t8code` thanks the Institute of Software Technology and the German Aerospace
Center (DLR).

The authors state that there are no conflicts of interest.

# References
