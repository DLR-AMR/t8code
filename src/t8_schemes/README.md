# t8_schemes Module

The `t8_schemes` module is a fundamental part of `t8code`'s combinatorial engine. A "scheme" defines the rules for how an element of a certain shape is refined and how its descendants (children, faces, etc.) are numbered and related to each other. For example, the scheme for a tetrahedron defines that it is refined into 8 child tetrahedra and describes the orientation and connectivity of those children.

`t8code` uses this information for all its tree-based navigation and neighborhood queries. By abstracting these rules into a scheme, the core algorithms in `t8_forest` can be written in a generic way that is independent of any specific element shape.

## Subfolders

- [**t8_default/**](t8_default/README.md): Contains the source code for the default, built-in refinement schemes used by `t8code` for standard element shapes like tetrahedra, hexahedra, prisms, and pyramids.
- [**t8_standalone/**](t8_standalone/README.md): Contains a standalone implementation of the schemes that can be compiled and used independently of the rest of `t8code`. This is useful for testing, debugging, and potentially for other tools that need to work with `t8code`'s refinement patterns.

## Key Files and Functionality

- **Scheme Interface:**
    - `t8_scheme.h`: The public C-API header. It defines the `t8_scheme_cxx_t` object and provides functions to query properties of a scheme, such as the number of children an element has, the number of faces, etc.
    - `t8_scheme.hxx`: The internal C++ header that defines the `t8_scheme_cxx` class, which is the C++ interface for a refinement scheme.
    - `t8_scheme.cxx`: The implementation of the C-API functions, which wrap the C++ class methods.

- **Scheme Construction:**
    - `t8_scheme_builder.hxx`: A header file containing a builder class used to construct a `t8_scheme_cxx` object from its raw connectivity data. This is an internal tool used when defining a new scheme.

## Main Purpose

The main purpose of this module is to provide a lookup-table-based mechanism for all combinatorial queries related to mesh refinement. Instead of having complex geometric logic scattered throughout the code, `t8code` can simply ask the scheme: "If I refine this face, what are the child faces?" or "What is the parent of child number 5?". This makes the core algorithms clean, fast, and extensible to new element types.

## Developer Entry Points

Most developers will not need to interact with this module directly. The `t8_forest` module uses the schemes internally. The only time a developer would work with this module is if they wanted to add support for a completely new element type with a custom refinement pattern. In that case, they would need to:
1.  Define the connectivity tables for their new element type.
2.  Use the `t8_scheme_builder` to construct a new `t8_scheme_cxx` instance.
3.  Register this new scheme with `t8code`.

The public query interface is in `t8_scheme.h`.
