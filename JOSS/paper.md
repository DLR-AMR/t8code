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
  - name: Johannes Markert
    orcid: 0000-0001-6297-9494
    corresponding: true
    equal-contrib: true
    affiliation: 1
  - name: David Knapp
    orcid: 0000-0002-6305-1572
    equal-contrib: true
    affiliation: 1
  - given-names: Lukas Dreyer
    orcid: 0000-0001-7484-3674
    equal-contrib: true
    affiliation: 1
  - given-names: Sandro Elsweijer
    orcid: 0000-0002-2753-3873
    equal-contrib: true
    affiliation: 1
  - given-names: Niklas Böing
    equal-contrib: true
    affiliation: 1
  - given-names: Chiara Hergl
    orcid: 0000-0002-4016-9113
    equal-contrib: true
    affiliation: 1
  - given-names: Prasanna Ponnusamy
    equal-contrib: true
    affiliation: 1
  - given-names: Jakob Fussbroich
    orcid: 0000-0003-0784-2182
    equal-contrib: true
    affiliation: 1
  - given-names: Tabea Leistikow
    equal-contrib: true
    affiliation: 1
  - given-names: Florian Becker
    orcid: 0000-0002-8384-9282
    equal-contrib: true
    affiliation: 1
  - given-names: Ioannis Lilikakis
    equal-contrib: true
    affiliation: 1
  - name: Carsten Burstedde
    orcid: 0000-0001-9843-1041
    affiliation: 2
affiliations:
 - name: German Aerospace Center (DLR), Institute for Software Technology, Cologne, Germany
   index: 1
 - name: Rheinische Friedrich-Wilhelms-Universität Bonn, Institute for Numerical Simulations and Hausdorff Center for Mathematics,  Germany
   index: 2
date: 10 June 2024
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
[www.dlr-amr.github.io/t8code](www.dlr-amr.github.io/t8code). It is developed
and maintained at the [Institute for Software Technology](https://www.dlr.de/sc/en/)
of the German Aerospace Center (DLR). The software library provides fast and memory
efficient parallel algorithms for dynamic AMR to handle tasks such as mesh
adaptation, load-balancing, ghost computation, feature search and more.
`t8code` can manage meshes with over one trillion mesh elements
[@holke_optimized_2021] and scales up to one million parallel processes
[@holke_scalable_2018]. It is intended to be used as mesh management backend in
scientific and engineering simulation codes paving the way towards
high-performance applications of the upcoming exascale era.

# Statement of Need

Adaptive Mesh Refinement has been established as a successful approach
for scientific and engineering simulations over the past decades
[@TEUNISSEN2019106866; @10.1145/1268776.1268779; @doi:10.1137/0733054;
@doi:10.1137/0715049]. By modifying the mesh resolution locally according to
problem specific indicators, the computational power is efficiently
concentrated where needed and the overall memory usage is reduced by orders of
magnitude. However, managing adaptive meshes and associated data is a very
challenging task, especially for parallel codes. Implementing fast and scalable
AMR routines generally leads to a large development overhead motivating the
need for external mesh management libraries like `t8code`.

Currently, `t8code`'s AMR routines support a wide range of element types:
vertices, lines, quadrilaterals, triangles, hexahedra, tetrahedra, prisms, and
pyramids. Additionally, implementation of other refinement patterns and element
shapes is possible.
See \autoref{fig:visploremesh} for an examplary adapted mesh managed by `t8code` for visualizing
earth mantle convection data.

![2D slice of an adapted `t8code` mesh for a visualization of earth mantle convection data.
\label{fig:visploremesh}](pics/visplore_magma_tilted_grid.png){width="70%"}

# Fundamental Concepts

`t8code` is based on the forest-of-trees approach. Starting point
for the usage of `t8code` is an unstructured input mesh, which
we denote a coarse mesh. This coarse mesh describes the geometry of the
computational domain. Each of the coarse mesh cells are then viewed as the
root of a refinement tree. These trees are refined recursively in a structured
pattern, resulting in a collection of trees, which we call a forest. `t8code`
stores only a minimal amount of information about the finest elements of the mesh -
the leaves of the trees - in order to reconstruct the whole forest.

By enumerating the leaves in a recursive refinement pattern we obtain a
space-filling curve (SFC) logic. Via these SFCs, all elements in a refinement
tree are assigned an integer-based index and are stored in linear order.
Element coordinates or element neighbors do not need to be stored explicitly
but can be reconstructed from the SFC index. Fast bitwise SFC operations ensure
optimal runtimes and diminish the need for memory lookups. Moreover, the SFC
is used to distribute the forest mesh across multiple processes, so that each
process only stores a unique portion of the SFC.  See
\autoref{fig:SpaceFillingCurves}.

While being successfully applied to quadrilateral
and hexahedral meshes [@burstedde_p4est_2011; @weinzierl_peano_2019],
these SFC techniques are extended by `t8code` in a modular fashion, such that arbitrary
element shapes are supported. We achieve this modularity through a novel
decoupling approach that separates high-level (mesh global) algorithms from
low-level (element local) implementations. All high-level algorithms can
be applied to different implementations of element shapes and refinement
patterns. A mix of different element shapes in the same mesh is also
supported.

![Left: Quad-tree of an exemplary forest mesh consisting of two trees (k0, k1)
distributed over three parallel processes p0 to p2. The SFC is represented by a
black curve tracing only the finest elements (leaves) of each tree. Right:
Sketch of the associated mixed shape mesh refined up to level three. Bottom
left: The elements saved by p1 and the associated ghost elements (non process
local neighbors).
\label{fig:SpaceFillingCurves}](pics/t8code_sfc_hybrid.png)

# Performance

`t8code` supports distributed coarse meshes of arbitrary size and complexity,
which we tested for up to 370 million coarse mesh cells
[@burstedde_coarse_2017].  Moreover, we conducted various performance studies
on the JUQUEEN and the JUWELS supercomputers at the Jülich Supercomputing
Center. `t8code`'s ghost and partition routines are exceptionally fast with
proper scaling of up to 1.1 trillion mesh elements; see
\autoref{tab:t8code_runtimes}, [@holke_optimized_2021].  Furthermore, in a
prototype code [@Dreyer2021] implementing a high-order discontinuous Galerkin
method (DG) for advection-diffusion equations on dynamically adaptive
hexahedral meshes we obverve a 12 times speed-up compared to non-AMR meshes
with only an overall 15\% runtime contribution of `t8code`; see
\autoref{fig:t8code_runtimes}. 

+----------------+-------------------+--------------------+--------+-----------+
| \# Process     | \# Elements       | \# Elem. / process | Ghost  | Partition |
+:==============:+:=================:+:==================:+:======:+:=========:+
|     49,152     | 1,099,511,627,776 |     22,369,621     | 2.08 s |   0.73 s  |
+----------------+-------------------+--------------------+--------+-----------+
|     98,304     | 1,099,511,627,776 |     11,184,811     | 1.43 s |   0.33 s  |
+================+===================+====================+========+===========+
| Table 1: Runtimes on JUQUEEN for the ghost layer and partitioning operations |
| for a distributed mesh consisting of 1.1 trillion elements.                  |
| \label{tab:t8code_runtimes}                                                  |
+================+===================+====================+========+===========+

![Runtimes on JUQUEEN of the different components of our DG prototype code
coupled with `t8code`. Note that the operations associated with dynamical mesh
adaptation (adapt, balance, partition, and ghost) utilize only around 15\% of
the total runtime largely independent of the number of
processes.\label{fig:t8code_runtimes}
](pics/t8code_runtimes_2.png){width="90%"}

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
(massive data visualization), [HYTAZER](https://elib.dlr.de/201347/) (hydrogen tank certification), [Greenstars](https://www.dlr.de/en/ra/research-transfer/projects/dlr-projects/green-propulsion-free-flight-demonstrator-the-lander)
(additive rocket engine manufacturing) and [PADME-AM](https://m2p2023.cimne.com/event/contribution/24a9d174-9411-11ed-b949-000c29ddfc0c) (simulation assisted
additive manufacturing).

# Further Information

For further information beyond this short note and also for code examples, we
refer to our
[Documentation](https://dlr-amr.github.io/t8code/pages/documentation.html) and
[Wiki](https://github.com/DLR-AMR/t8code/wiki) reachable via our homepage
[dlr-amr.github.io/t8code](https://dlr-amr.github.io/t8code/) and our technical
publications on `t8code` [@holke_scalable_2018; @burstedde_coarse_2017;
@holke_optimized_2021; @burstedde_tetrahedral_2016; @Knapp20;
@Becker_hanging_faces; @elsweijer_curved_2021; @Dreyer2021;
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
`t8code` thanks the Institute for Software Technology and the German Aerospace
Center (DLR).

The authors state that there are no conflicts of interest.

# References
