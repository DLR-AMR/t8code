# Pyramid Default Scheme

This directory contains the implementation of the default refinement scheme for a **pyramid** element. The refinement of a pyramid is particularly complex compared to other element types.

## Functionality

The files in this directory define the precise combinatorial rules for the refinement of a pyramid-shaped element. Unlike other element types which refine into children of the same shape, the `t8code` pyramid scheme refines a parent pyramid into a collection of **10 child elements: 6 tetrahedra and 4 pyramids**.

This directory encodes all the information about this complex refinement, including:
-   The type, orientation, position, and numbering of the 10 children relative to their parent.
-   The relationship between the parent's faces and the children's faces.
-   The connectivity between the child elements themselves.

This information is stored in a series of static lookup tables.

## Files

-   `t8_default_pyramid.cxx`/`.hxx`: These files contain the primary logic that uses the `t8_scheme_builder` to construct the final scheme object from the raw connectivity data.
-   `t8_dpyramid.h`: Public header for the default pyramid scheme, providing access to its specific properties.
-   `t8_dpyramid_connectivity.c`/`.h`: These files contain the large, explicit, hard-coded lookup tables that define the complex connectivity of the 10-child refinement pattern.
-   `t8_dpyramid_bits.c`/`.h`: These files contain low-level, bitwise manipulation functions that implement the logic of the space-filling curve adapted for the pyramid's topology.

## Main Purpose

This is a low-level, internal component of `t8code`, essential for supporting hybrid meshes containing pyramids. The complexity of the files reflects the inherent difficulty in conformally refining a pyramid. The `t8_forest` module consults these lookup tables whenever it needs to perform a topological operation on a pyramid element. Developers typically do not need to interact with these files directly.
