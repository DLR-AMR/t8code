# Prism Default Scheme

This directory contains the implementation of the default refinement scheme for a **prism** element.

## Functionality

The files in this directory define the precise combinatorial rules for the refinement of a prism-shaped element. Specifically, they encode the following information:

-   A prism is refined into 8 child prisms.
-   The orientation, position, and numbering of these 8 children relative to their parent.
-   The relationship between the parent's faces and the children's faces.
-   The connectivity between the child elements themselves.

This information is stored in a series of static lookup tables.

## Files

-   `t8_default_prism.cxx`/`.hxx`: These files contain the primary static arrays that define the prism scheme. They use the `t8_scheme_builder` to construct the final scheme object from these arrays.
-   `t8_dprism.h`: Public header for the default prism scheme, providing access to its specific properties.
-   `t8_dprism_bits.c`/`.h`: These files contain low-level, bitwise manipulation functions that implement the logic of the space-filling curve for the prism's unique topology (a product of a triangle and a line).

## Main Purpose

This is a low-level, internal component of `t8code`, essential for supporting hybrid meshes containing prisms. The `t8_forest` module consults these lookup tables whenever it needs to perform a topological operation on a prism element. Developers typically do not need to interact with these files directly.
