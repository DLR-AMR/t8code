# Quadrangle Default Scheme

This directory contains the implementation of the default refinement scheme for a **quadrangle (quad)** element.

## Functionality

The files in this directory define the precise combinatorial rules for the refinement of a quadrilateral element. Specifically, they encode the following information:

-   A quad is refined into 4 child quads.
-   The orientation, position, and numbering of these 4 children relative to their parent.
-   The relationship between the parent's edges and the children's edges.
-   The connectivity between the child elements themselves.

This information is stored in a series of static lookup tables.

## Files

-   `t8_default_quad.cxx`/`.hxx`: These files contain the primary static arrays that define the quad scheme and use the `t8_scheme_builder` to construct the final scheme object.
-   `t8_dquad.h`: Public header for the default quad scheme, providing access to its specific properties.
-   `t8_default_quad_bits.cxx`/`.hxx`: These files contain low-level, bitwise manipulation functions. The refinement pattern for a quad is isomorphic to a 2D Morton (or Z-order) curve, which is implemented here with efficient bitwise operations.

## Main Purpose

This is a low-level, internal component of `t8code`. The `t8_forest` module consults these lookup tables whenever it needs to perform a topological operation (like finding a neighbor or a child) on a quad element. This is essential for 2D simulations or for handling quad faces in 3D. Developers typically do not need to interact with these files directly.
