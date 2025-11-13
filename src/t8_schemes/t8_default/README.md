# Default Schemes

This directory provides the concrete implementations for `t8code`'s standard, built-in refinement schemes. Each supported element shape (tetrahedron, hexahedron, etc.) has its own subdirectory containing the raw, hard-coded connectivity data that defines its refinement pattern.

These schemes are initialized once when `t8code` is set up and are then made available to the rest of the library.

## Subfolders

Each subfolder contains the raw connectivity data (in the form of static C arrays) for a specific element shape, known in `t8code` as an `eclass`.

- [**t8_default_common/**](t8_default_common/README.md): Common data and definitions shared by all default schemes.
- [**t8_default_line/**](t8_default_line/README.md): Connectivity data for a 1D Line element.
- [**t8_default_quad/**](t8_default_quad/README.md): Connectivity data for a 2D Quadrangle element.
- [**t8_default_tri/**](t8_default_tri/README.md): Connectivity data for a 2D Triangle element.
- [**t8_default_hex/**](t8_default_hex/README.md): Connectivity data for a 3D Hexahedron element.
- [**t8_default_prism/**](t8_default_prism/README.md): Connectivity data for a 3D Prism element.
- [**t8_default_pyramid/**](t8_default_pyramid/README.md): Connectivity data for a 3D Pyramid element.
- [**t8_default_tet/**](t8_default_tet/README.md): Connectivity data for a 3D Tetrahedron element.
- [**t8_default_vertex/**](t8_default_vertex/README.md): Connectivity data for a 0D Vertex element.

## Key Files

- `t8_default_c_interface.h`: The public C-API header that provides access to the default schemes. It contains functions like `t8_scheme_new_default_tet()` which return a `t8_scheme_cxx_t` object for a given element type.
- `t8_default.hxx`: The internal C++ header that orchestrates the creation of the default scheme objects from the raw data in the subdirectories.
- `t8_default.cxx`: The C++ implementation file for creating the default schemes.

## Main Purpose

This directory is the ultimate source of truth for the combinatorial topology of `t8code`. The static arrays in the subdirectories are the lookup tables that the `t8_forest` module consults for all its operations. For example, when a hexahedron is refined, `t8_forest` queries the `t8_default_hex` scheme to find out the coordinates and face connectivity of the 8 child hexahedra.

A deep understanding of these files is only necessary for developers who are debugging the lowest levels of the refinement logic or who are interested in the mathematical details of the space-filling curve and refinement patterns.
