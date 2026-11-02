# src/t8_forest/t8_forest_boundary_conditions

This folder provides `t8_forest_get_boundary_conditions` and `t8_forest_get_boundary_condition` for
C++ and C. They return the boundary condition labels of the faces of a forest element.

Most of the boundary condition module is implemented in the cmesh. Therefore, a more detailed description can be found
over in the [README of the cmesh boundary condition module](../../t8_cmesh/t8_cmesh_boundary_conditions/README.md)

## Files

#### [t8_forest_boundary_conditions.hxx](t8_forest_boundary_conditions.hxx)

This header contains the public C++ interface. `t8_forest_get_boundary_conditions` returns the labels
of all faces of an element and `t8_forest_get_boundary_condition` returns the label of a single face.
If an element face lies inside its tree, it has no label and an empty `std::optional` is returned.

#### [t8_forest_boundary_conditions_c_interface.h](t8_forest_boundary_conditions_c_interface.h)

This header provides the same functions for C. They return `const char *` labels, which are `NULL`
for element faces that lie inside their tree.
