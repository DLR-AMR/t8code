# t8code Benchmarks

This directory contains a collection of benchmark programs designed to measure the performance and scalability of various components within the `t8code` library. These benchmarks are crucial for performance regression testing, identifying bottlenecks, and demonstrating the efficiency of `t8code`'s core algorithms.

## Building the Benchmarks

To build the benchmarks, you need to enable the `T8_ENABLE_BENCHMARKS` option in your CMake configuration.

## Benchmark Programs

Each file is a standalone program that times a specific operation or workflow.

- **`t8_time_forest_partition.cxx`**: This benchmark measures the time it takes to partition a `t8_forest`. It  creates a forest, possibly adapts it to create a load imbalance, and then times the `t8_forest_partition` call. This is a key benchmark for assessing parallel scalability.

- **`t8_time_fractal.cxx`**: This program benchmarks the creation of a "fractal" mesh. This typically involves recursively adapting the mesh based on a geometric criterion, which heavily exercises the `t8_forest_adapt` and `t8_forest_balance` functions. It's a good measure of the raw speed of mesh adaptation.

- **`t8_time_new_refine.c`**: A C-based benchmark that specifically times the refinement process (`t8_forest_adapt`).

- **`t8_time_partition.cxx`**: This is another benchmark for partitioning to time the partition. It might test a different scenario or a different aspect of partitioning compared to `t8_time_forest_partition.cxx`.

- **`t8_time_prism_adapt.cxx`**: This benchmark specifically measures the performance of adaptation on a mesh composed of prisms. This is important for ensuring that the performance of hybrid mesh elements is also tracked.

- **`t8_time_set_join_by_vertices.cxx`**: This benchmark times the `t8_cmesh_set_join_by_vertices` function. This function is part of setting up the coarse mesh connectivity and can have a performance impact, especially for very large coarse meshes.

## How to Run

After building the benchmarks, you can run the individual executables. Most of them are intended to be run with MPI to test parallel performance.

```bash
# Example of running a benchmark with multiple processes
mpirun -n 16 ./t8_time_forest_partition
```

The programs will typically print timing information and other performance metrics to the console. The exact parameters (e.g., initial mesh size, refinement levels) can often be controlled via command-line arguments. Run a benchmark with the `-h` flag to see its available options.
