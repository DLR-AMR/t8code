# Standalone Schemes

This directory contains a "standalone" implementation of the `t8code` refinement schemes. This version is designed to be compiled and used independently of the main `t8code` library, with minimal dependencies.

## Main Purpose

The primary purpose of this standalone version is to allow external tools, testing frameworks, or other applications to access the combinatorial information of `t8code`'s refinement patterns without needing to link against the full library (which includes MPI, `p4est`, etc.).

Potential use cases include:
-   **Unit Testing:** The schemes' logic can be tested in a lightweight, isolated environment.
-   **External Tools:** A post-processor or a mesh converter could use this to understand the structure of a `t8code` mesh.
-   **Visualization and Debugging:** A simple program can be built using this module to visualize the refinement of a single element, which is useful for debugging and educational purposes.

## Files

-   `t8_standalone.hxx`: The main public header for the standalone schemes. It provides the interface to access the scheme information.
-   `t8_standalone.cxx`: The implementation file for the standalone scheme interface.
-   `t8_standalone_elements.hxx`: Defines the element types (`eclass`) in a context that is independent of the main `t8code` library.
-   `t8_standalone_implementation.hxx`: This header contains the actual implementation. It likely includes the raw connectivity data from the `t8_default` schemes and uses the `t8_scheme_builder` to construct the scheme objects, but without any of the other `t8code` components.

## How to Use

A developer would include `t8_standalone.hxx` in their application and link against the object file compiled from `t8_standalone.cxx`. They could then get a scheme object (e.g., for a tetrahedron) and query it for its properties, just as one would with the full `t8_schemes` module.
