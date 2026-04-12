# Geometry Examples

This folder contains examples that demonstrate how to use different geometry representations within `t8code`. Associating a geometry with a mesh is essential for simulations on domains that are not simple cubes.

## Files

- `t8_example_geometries.cxx`: This program showcases how to create and use various types of geometries. It likely includes examples of:
    - **Analytical Geometries:** Geometries defined by mathematical formulas, such as a sphere or a cylinder.
    - **CAD Geometries:** Geometries loaded from a CAD file (e.g., via OpenCASCADE), which allows for arbitrarily complex shapes.
    - **Linear Geometries:** The default, where all boundaries are treated as flat polygons.

The example probably demonstrates how `t8code` can use this geometric information to correctly handle mesh-to-geometry conformity, which is crucial for accurate boundary conditions and simulations.

## Main Purpose

This example is a guide for users who need to run simulations on complex domains. It shows how to set up the geometry object and associate it with the `t8_cmesh`, which is a prerequisite for creating a geometry-conforming adaptive mesh.

## How to Run

After building the examples, you can run the executable. Note that to run the CAD-based examples, you must have compiled `t8code` with support for a CAD kernel like OpenCASCADE.

```bash
./t8_example_geometries
```

The program will likely print information to the console about the geometries it creates and the operations it performs.
