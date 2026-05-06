# Concrete Geometry Implementations

This folder contains the source code for several concrete implementations of the `t8_geometry` abstract base class. Each implementation represents a different way to define a geometric domain. Users can choose the implementation that best suits their needs.

The public functions to create instances of these geometries (e.g., `t8_geometry_new_analytic`) are declared in the `.h` files in this directory.

## Geometry Implementations

- **`t8_geometry_analytic`**:
    - `t8_geometry_analytic.h`/`.cxx`/`.hxx`
    - This implementation is for geometries that can be described by a mathematical formula. The user provides function pointers that define the geometry, such as a level-set function (`phi(x,y,z) = 0` for the surface) and its gradient. This is useful for simple, well-defined shapes like spheres, cylinders, or tori.

- **`t8_geometry_cad`**:
    - `t8_geometry_cad.h`/`.cxx`/`.hxx`
    - This is a powerful implementation that interfaces with a CAD (Computer-Aided Design) kernel, such as OpenCASCADE. It allows `t8code` to use complex, industry-standard CAD models as the simulation domain. This is the preferred method for complex, real-world geometries. Note that this requires `t8code` to be compiled with CAD support.

- **`t8_geometry_lagrange`**:
    - `t8_geometry_lagrange.h`/`.cxx`/`.hxx`
    - This implementation represents the geometry using high-order Lagrange polynomials. This is used to represent curved surfaces with high fidelity and is often used in high-order numerical methods.

- **`t8_geometry_linear`**:
    - `t8_geometry_linear.h`/`.cxx`/`.hxx`
    - This represents the geometry as a collection of flat, linear facets. It is a piecewise linear representation of the domain.
    - `t8_geometry_linear_axis_aligned.h`/`.cxx`/`.hxx`: A specialized and optimized version for linear geometries that are aligned with the coordinate axes.

- **`t8_geometry_zero`**:
    - `t8_geometry_zero.h`/`.cxx`/`.hxx`
    - This is a "null" or "do-nothing" geometry implementation. It essentially provides no geometric information. This is used as a default when no other geometry is specified, and it effectively means the domain is defined purely by the straight-edged elements of the coarse mesh.

- **`t8_geometry_examples`**:
    - `t8_geometry_examples.h`/`.cxx`/`.hxx`
    - This file does not provide a new geometry implementation itself, but rather contains functions that create instances of various example geometries (like a brick, sphere, etc.) using the other implementations. This is used by the example programs and for testing.

## Developer Entry Points

To use one of these geometries, a developer would include the relevant header file (e.g., `t8_geometry_analytic.h`) and call its constructor function (e.g., `t8_geometry_new_analytic(...)`), passing the required parameters. The resulting `t8_geometry_t*` can then be passed to `t8_cmesh_set_geometry`.
