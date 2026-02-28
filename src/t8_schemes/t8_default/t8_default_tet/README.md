# Tetrahedron Default Scheme

This directory contains the implementation of the default refinement scheme for a **tetrahedron**. This is one of the most important schemes for 3D simulations.

## Functionality

The files in this directory define the precise combinatorial rules for the refinement of a tetrahedral element. Specifically, they encode the following information:

-   A tetrahedron is refined into 8 child tetrahedra.
-   The orientation, position, and numbering of these 8 children relative to their parent.
-   The relationship between the parent's triangular faces and the children's faces.
-   The connectivity between the child elements themselves.

This information is stored in a series of static lookup tables.

## Files

-   `t8_default_tet.cxx`/`.hxx`: These files contain the primary logic that uses the `t8_scheme_builder` to construct the final scheme object from the raw connectivity data.
-   `t8_dtet.h`: Public header for the default tetrahedron scheme, providing access to its specific properties.
-   `t8_dtet_connectivity.c`/`.h`: These files contain the large, explicit, hard-coded lookup tables that define the connectivity of the 8-child refinement pattern.
-   `t8_dtet_bits.c`/`.h`: These files contain low-level, bitwise manipulation functions that implement the logic of the space-filling curve for the tetrahedron's topology.
-   `t8_dtri_to_dtet.h`: A helper header that defines the mapping from the 2D triangle scheme (used for faces) to the 3D tetrahedron scheme. This is crucial for correctly handling face-neighbor relationships.

## Main Purpose

This is a low-level, internal component of `t8code`. The `t8_forest` module consults these lookup tables whenever it needs to perform a topological operation (like finding a neighbor or a child) on a tetrahedral element. The use of pre-computed tables makes these operations extremely fast. Developers typically do not need to interact with these files directly.
