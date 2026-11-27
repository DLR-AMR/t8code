# src/t8_cmesh/t8_cmesh_io

*NOTE: This README file is still to be written / finished.*

#### [t8_cmesh_readmshfile.h](t8_cmesh_readmshfile.h)

Input routines for reading `.msh` files (usually from the Gmsh mesh generator).

#### [t8_cmesh_cgns_reader.hxx](t8_cmesh_cgns_reader.hxx)

Reader for CGNS files (only built with `T8CODE_ENABLE_CGNS`). `t8_cmesh_from_cgns_file` reads a file
into a committed cmesh. The class `t8_cmesh_cgns_reader` additionally exposes the zones and boundary
conditions that were found. The reader supports unstructured zones with linear elements, zone interfaces
and periodic interfaces (`ZoneGridConnectivity`, type `Abutting1to1`) and boundary conditions (`ZoneBC`).
The zone index of each tree is stored as a tree attribute (`t8_cmesh_get_cgns_zone`), the boundary
conditions are stored as boundary condition labels of the tree faces
(see [t8_cmesh_boundary_conditions](../t8_cmesh_boundary_conditions/README.md)).