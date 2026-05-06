# Triangle Default Scheme

This directory contains the implementation of the default refinement scheme for a **triangle**. This scheme is fundamental for 2D simulations and for representing the faces of 3D elements like tetrahedra and prisms.

## Functionality

The files in this directory define the precise combinatorial rules for the refinement of a triangular element. Specifically, they encode the following information:

-   A triangle is refined into 4 child triangles.
-   The orientation, position, and numbering of these 4 children relative to their parent.
-   The relationship between the parent's edges and the children's edges.
-   The connectivity between the child elements themselves.

This information is stored in a series of static lookup tables.

## Files

-   `t8_default_tri.cxx`/`.hxx`: These files contain the primary logic that uses the `t8_scheme_builder` to construct the final scheme object from the raw connectivity data.
-   `t8_dtri.h`: Public header for the default triangle scheme, providing access to its specific properties.
-   `t8_dtri_connectivity.c`/`.h`: These files contain the explicit, hard-coded lookup tables that define the connectivity of the 4-child refinement pattern.
-   `t8_dtri_bits.c`/`.h`: These files contain low-level, bitwise manipulation functions that implement the logic of the space-filling curve for the triangle's topology.

## Main Purpose

This is a low-level, internal component of `t8code`. The `t8_forest` module consults these lookup tables whenever it needs to perform a topological operation on a triangular element or face. The use of pre-computed tables makes these operations extremely fast. Developers typically do not need to interact with these files directly.
