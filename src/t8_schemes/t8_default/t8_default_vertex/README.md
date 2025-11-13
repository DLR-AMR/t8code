# Vertex Default Scheme

This directory contains the implementation of the default refinement scheme for a **0D Vertex** element.

## Functionality

This is the simplest possible scheme. A vertex represents a point in space and has no dimension. The files in this directory define the combinatorial rules for a vertex, which are trivial:

-   A vertex is refined into a single child vertex (itself).
-   It has no faces or edges.

This scheme is necessary for the completeness and consistency of the `t8_schemes` framework, providing a base case for recursive algorithms.

## Files

-   `t8_default_vertex.cxx`/`.hxx`: These files contain the minimal data and logic required to build the vertex scheme object.
-   `t8_dvertex.h`: Public header for the default vertex scheme.

## Main Purpose

This is a low-level, internal component of `t8code` that handles the topological properties of a 0D element. It is primarily used to represent the corners of 1D, 2D, and 3D elements. Developers typically do not need to interact with these files directly.
