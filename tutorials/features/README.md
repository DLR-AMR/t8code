# Feature Tutorials

This directory contains tutorials that focus on specific, advanced features of `t8code`. Unlike the general tutorials, these are not necessarily meant to be followed in a specific order. Instead, they serve as deep dives into particular topics.

---

### [Feature: Curved Meshes](https://github.com/DLR-AMR/t8code/wiki/Feature---Curved-meshes)

This tutorial demonstrates one of `t8code`'s most powerful features: the ability to create and manage adaptive meshes that conform to curved domain boundaries.

- **File:** `t8_features_curved_meshes.cxx`
- **Goal:** To show how to:
    1.  Create a `t8_geometry` object that describes a curved domain (in this case, a sphere).
    2.  Associate this geometry with a `t8_cmesh`.
    3.  Create a `t8_forest` from this geometry-aware cmesh.
    4.  Observe how `t8code` automatically projects the nodes of refined elements onto the true curved surface, resulting in a mesh that accurately represents the geometry, rather than being a coarse, blocky approximation.

#### Input Geometry Files

The tutorial can be run with different coarse meshes to show its flexibility. This directory includes several Gmsh geometry (`.geo`) files that can be used to generate suitable input meshes.

-   `t8_features_curved_meshes_generate_cmesh_hex.geo`: A 3D coarse mesh of hexahedra.
-   `t8_features_curved_meshes_generate_cmesh_tet.geo`: A 3D coarse mesh of tetrahedra.
-   `t8_features_curved_meshes_generate_cmesh_quad.geo`: A 2D coarse mesh of quadrilaterals.
-   `t8_features_curved_meshes_generate_cmesh_tri.geo`: A 2D coarse mesh of triangles.

#### How to Run

1.  **Generate a mesh:** Use Gmsh to create a `.msh` file from one of the `.geo` files.
    ```bash
    gmsh t8_features_curved_meshes_generate_cmesh_tet.geo -3 -o tet_mesh.msh
    ```
2.  **Run the tutorial:** Execute the tutorial program, passing the path to the generated mesh file.
    ```bash
    ./t8_features_curved_meshes --mesh-file tet_mesh.msh
    ```

The program will output a VTK file (e.g., `t8_features_curved_mesh.vtu`) which you can open in ParaView or VisIt to see the resulting curved, adaptive mesh.
