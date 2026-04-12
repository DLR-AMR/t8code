# General Tutorials

This directory contains the general, step-by-step tutorials for learning `t8code`. They are designed to be followed in order, as each step builds upon the concepts introduced in the previous one.

For a detailed walkthrough of each tutorial, please see the corresponding article on the [t8code Wiki](https://github.com/DLR-AMR/t8code/wiki/Tutorial---Overview).

**Note on File Structure:** Many tutorials are split into multiple files. Typically, there is a `_main.cxx` file that contains the main function and entry point, and other `.cxx` and `.h` files that contain the core logic for the tutorial. This is done to keep the code organized and focused.

---

### [Step 0: Hello World](https://github.com/DLR-AMR/t8code/wiki/Step-0---Hello-World)
- **File:** `t8_step0_helloworld.cxx`
- **Goal:** The absolute basics: how to initialize and finalize the `t8code` library and its dependency, `p4est`.

### [Step 1: Creating a Coarse Mesh](https://github.com/DLR-AMR/t8code/wiki/Step-1---Creating-a-coarse-mesh)
- **File:** `t8_step1_coarsemesh.cxx`
- **Goal:** Learn how to create a `t8_cmesh_t` object, the coarse mesh that serves as the basis for the adaptive forest. This tutorial creates a cmesh from a built-in example and writes it to a VTK file.

### [Step 2: Creating a Uniform Forest](https://github.com/DLR-AMR/t8code/wiki/Step-2---Creating-a-uniform-forest)
- **File:** `t8_step2_uniform_forest.cxx`
- **Goal:** Learn how to create a `t8_forest_t` from a `t8_cmesh`. This tutorial creates a uniformly refined forest where all elements have the same refinement level.

### [Step 3: Adapting a Forest](https://github.com/DLR-AMR/t8code/wiki/Step-3---Adapting-a-forest)
- **Files:** `t8_step3_main.cxx`, `t8_step3_adapt_forest.cxx`, `t8_step3.h`
- **Goal:** Introduce adaptive mesh refinement (AMR). This tutorial shows how to provide a callback function to `t8_forest_adapt` to selectively refine elements based on a user-defined criterion.

### [Step 4: Partition, Balance, Ghost](https://github.com/DLR-AMR/t8code/wiki/Step-4---Partition,-Balance,-Ghost)
- **Files:** `t8_step4_main.cxx`, `t8_step4_partition_balance_ghost.cxx`, `t8_step4.h`
- **Goal:** Explain the key concepts for parallel execution. This tutorial covers partitioning the forest for load balancing, balancing the forest to maintain the 2:1 constraint, and creating the ghost layer for communication.

### [Step 5: Storing Element Data](https://github.com/DLR-AMR/t8code/wiki/Step-5---Store-element-data)
- **Files:** `t8_step5_main.cxx`, `t8_step5_element_data.cxx`, `t8_step5.h`, `t8_step5_element_data_c_interface.c`
- **Goal:** Show how to associate application-specific data with mesh elements using the `t8_data_handler`. It also demonstrates how to exchange this data between processes using the ghost layer.

### [Step 6: Computing Stencils](https://github.com/DLR-AMR/t8code/wiki/Step-6-Computing-stencils)
- **Files:** `t8_step6_main.cxx`, `t8_step6_stencil.cxx`, `t8_step6.h`
- **Goal:** Demonstrate how to use the `t8_forest` iterators to traverse the mesh and access data from face-neighbors. This is a fundamental operation for implementing numerical schemes like finite volume or finite difference methods.

### [Step 7: Interpolation](https://github.com/DLR-AMR/t8code/wiki/Tutorial-Step-7---Interpolation)
- **Files:** `t8_step7_main.cxx`, `t8_step7_interpolation.cxx`, `t8_step7.h`
- **Goal:** Show how to handle data after the mesh has been adapted. This tutorial demonstrates how to define an interpolation operator to transfer data from parent elements to their new children after refinement.

---

## Other General Tutorials

### [Building a Custom cmesh](https://github.com/DLR-AMR/t8code/wiki/Build-Cmesh)
- **Files:** `t8_tutorial_build_cmesh_main.cxx`, `t8_tutorial_build_cmesh.cxx`, `t8_tutorial_build_cmesh.h`
- **Goal:** Learn how to programmatically create your own `t8_cmesh` from scratch by defining the vertices and element connectivity manually.

### [Building a Custom Scheme](https://github.com/DLR-AMR/t8code/wiki/Tutorial:-Build-Scheme)
- **File:** `t8_tutorial_build_scheme.cxx`
- **Goal:** An advanced tutorial that shows how to define a new refinement scheme for a custom element type.

### [Search](https://github.com/DLR-AMR/t8code/wiki/Tutorial:-Search)
- **File:** `t8_tutorial_search.cxx`
- **Goal:** Demonstrate how to use the `t8_forest_search` functions to find which elements in the forest contain a given set of points.
