# Test Files

This folder contains various data files that are used as input for the tests in the `t8code` test suite. These files provide a consistent basis for testing functionalities like file I/O and mesh generation.

## Files

The files in this directory are examples of different mesh and data formats that `t8code` can interact with.

### MSH Files (`.msh`)

These are mesh files in the Gmsh `.msh` format. They are used to test the `t8_cmesh_readmshfile` functionality.
- `test_msh_file_vers2_ascii.msh`: An ASCII-formatted Gmsh version 2 file.
- `test_msh_file_vers2_bin.msh`: A binary-formatted Gmsh version 2 file.
- `test_msh_file_vers4_ascii.msh`: An ASCII-formatted Gmsh version 4 file.
- `test_msh_file_vers4_bin.msh`: A binary-formatted Gmsh version 4 file.

### VTK Files (`.vtu`, `.vtp`, `.pvtu`, `.pvtp`)

These are files in the Visualization Toolkit (VTK) format. They are used to test `t8code`'s VTK reader and writer.
- `test_vtk_cube.vtp`: A single VTK PolyData file representing a cube.
- `test_vtk_tri.vtu`: A single VTK UnstructuredGrid file.
- `test_parallel_file.pvtu`: A parallel VTK UnstructuredGrid file, which points to the individual partition files.
- `test_parallel_file_0.vtu`, `test_parallel_file_1.vtu`: The individual data files for the parallel unstructured grid.
- `test_polydata.pvtp`: A parallel VTK PolyData file.
- `test_polydata_0.vtp`, `test_polydata_1.vtp`: The individual data files for the parallel PolyData.

### Other Input Files (`.inp`)

- `test_cube_unstructured_1.inp`, `test_cube_unstructured_2.inp`: These appear to be input files for generating unstructured cube meshes, possibly in a different format.

## Main Purpose

Having a dedicated set of test files ensures that the tests are reproducible and can be run without external dependencies on specific mesh generation tools. They are essential for verifying the I/O capabilities of `t8code`.
