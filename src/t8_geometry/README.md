# t8_geometry Module

The `t8_geometry` module provides the framework for representing and interacting with the geometric description of the simulation domain. `t8code` is designed to work with domains that are more complex than simple cubes, including those with curved boundaries. This module provides a flexible, abstract interface that allows the core mesh algorithms to query geometric information without being tied to a specific geometry representation.

The central concept is the `t8_geometry_t` object, which is an abstract representation of a geometry. Concrete implementations (e.g., for an analytical sphere, or a complex shape from a CAD file) inherit from a common base class and implement a standard set of geometric queries.

## Subfolders
- [**t8_geometry_implementations/**](t8_geometry_implementations/README.md): Contains the source code for various concrete implementations of the `t8_geometry` interface.

## Key Files and Functionalities

- **Abstract Geometry Interface:**
    - `t8_geometry.h`: The main public C-API header for the geometry module. It defines the `t8_geometry_t` type and the functions that users can call to interact with a geometry object.
    - `t8_geometry.cxx`: The implementation of the public API functions. These functions typically delegate their work to the virtual functions defined in the base class.
    - `t8_geometry_base.h`/`.cxx`/`.hxx`: Defines the abstract base class for all geometry implementations. It specifies the virtual function table (vtable) that each concrete geometry must provide. This includes functions for evaluating surface parameterizations, checking if a point is inside the domain, and other geometric queries.

- **Geometry Management:**
    - `t8_geometry_handler.cxx`/`.hxx`: Implements a handler class that is responsible for managing the lifecycle of geometry objects.
    - `t8_geometry_hash.hxx`: Provides hashing functions for geometries, which can be used to efficiently compare or store them.

- **Specific Geometry Types:**
    - `t8_geometry_with_vertices.h`/`.cxx`/`.hxx`: A specialized geometry implementation that is defined explicitly by a set of vertices. This is often used to represent piecewise linear geometries.

- **Utilities:**
    - `t8_geometry_helpers.c`/`.h`: A collection of helper and utility functions for geometric calculations.
    - `t8_geometry_extended.hxx`: Contains declarations for more advanced or extended geometry functionalities, often for internal C++ usage.

## Developer Workflow

A developer using this module would typically:
1.  Create an instance of a concrete geometry (e.g., `t8_geometry_new_analytic` or `t8_geometry_new_occ`). The functions to do this are found in the `t8_geometry_implementations` module.
2.  Associate this `t8_geometry_t` object with a `t8_cmesh`.
3.  The `t8_forest` will then use this geometry information during adaptation and other operations to ensure that the mesh conforms to the domain boundaries. For example, when a new child element is created on a curved boundary, its nodes will be projected onto the true surface representation provided by the `t8_geometry` object.

The public API is defined in `t8_geometry.h`. The interface that new geometry implementations must satisfy is defined in `t8_geometry_base.h`.
