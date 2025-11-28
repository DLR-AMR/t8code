/**
 * @file t8_mra_unified_plotting_example.cxx
 * @brief Example showing MRA with VTK plotting
 *
 * Demonstrates:
 * 1. Triangle and quad MRA with unified interface
 * 2. VTK output generation at different stages
 * 3. Adaptive refinement/coarsening with visualization
 * 4. Comparison between uniform and adapted meshes
 */

#ifdef T8_ENABLE_MRA

#include "t8_mra/t8_mra_unified_triangle.hpp"
#include "t8_mra/t8_mra_unified_cartesian.hpp"
#include "t8_mra/t8_mra_vtk_unified.hpp"
#include "t8_cmesh.h"
#include "t8_cmesh/t8_cmesh_examples.h"

#include <iostream>
#include <cmath>
#include <string>

//=============================================================================
// Test Functions
//=============================================================================

/**
 * @brief Gaussian bump function
 */
template <int U>
auto
gaussian_bump ()
{
  return [] (double x, double y) -> std::array<double, U> {
    const double cx = 0.5, cy = 0.5;
    const double sigma = 0.15;
    const double r2 = (x - cx) * (x - cx) + (y - cy) * (y - cy);
    return { std::exp (-r2 / (2 * sigma * sigma)) };
  };
}

/**
 * @brief Sine wave function
 */
template <int U>
auto
sine_wave ()
{
  return
    [] (double x, double y) -> std::array<double, U> { return { std::sin (4 * M_PI * x) * std::sin (4 * M_PI * y) }; };
}

/**
 * @brief Step function (discontinuous)
 */
template <int U>
auto
step_function ()
{
  return [] (double x, double y) -> std::array<double, U> { return { (x > 0.5 && y > 0.5) ? 1.0 : 0.0 }; };
}

//=============================================================================
// VTK Output Helpers
//=============================================================================

/**
 * @brief Write VTK output for a given MRA state
 */
template <typename MRA>
void
write_vtk_output (MRA &mra, const std::string &prefix, int step)
{
  const std::string filename = prefix + "_step" + std::to_string (step);

  std::cout << "  Writing VTK: " << filename << ".vtu" << std::endl;

  // Write high-order Lagrange VTK (P-1 is polynomial degree)
  t8_mra::write_forest_lagrange_vtk_unified (mra, filename.c_str (), MRA::P_DIM - 1);
}

//=============================================================================
// Example 1: Triangle MRA with Adaptation and Plotting
//=============================================================================

void
example_triangle_adaptive_with_plotting ()
{
  std::cout << "\n";
  std::cout << "════════════════════════════════════════════════════════════\n";
  std::cout << "  Triangle MRA: Adaptive Refinement with VTK Output\n";
  std::cout << "════════════════════════════════════════════════════════════\n\n";

  // MRA parameters
  constexpr int U = 1;
  constexpr int P = 3;
  const int min_level = 0;
  const int max_level = 6;
  const double c_thresh = 1.0;
  const int gamma = 1;
  const int dunavant_rule = 10;
  const bool balanced = false;

  // Create multiscale object
  t8_mra::multiscale<T8_ECLASS_TRIANGLE, U, P> mra (max_level, c_thresh, gamma, dunavant_rule, balanced,
                                                    sc_MPI_COMM_WORLD);

  // Create mesh
  t8_cmesh_t cmesh = t8_cmesh_new_hypercube (T8_ECLASS_TRIANGLE, sc_MPI_COMM_WORLD, 0, 0, 0);
  auto *scheme = t8_scheme_new_default ();

  // Initialize with Gaussian bump
  std::cout << "1. Initializing uniform mesh at level " << max_level << "...\n";
  auto func = gaussian_bump<U> ();
  mra.initialize_data (cmesh, scheme, max_level, func);

  auto num_elements = t8_forest_get_global_num_leaf_elements (mra.get_forest ());
  std::cout << "   Elements: " << num_elements << "\n";
  std::cout << "   Total DOF: " << (num_elements * mra.DOF) << "\n\n";

  // Write initial uniform solution
  write_vtk_output (mra, "unified/triangle_uniform", 0);

  // Perform adaptive coarsening
  std::cout << "2. Performing adaptive coarsening...\n";
  mra.coarsening_new (min_level, max_level);

  num_elements = t8_forest_get_global_num_leaf_elements (mra.get_forest ());
  std::cout << "\n   After coarsening:\n";
  std::cout << "   Elements: " << num_elements << "\n";
  std::cout << "   Total DOF: " << (num_elements * mra.DOF) << "\n\n";

  // Write coarsened solution
  write_vtk_output (mra, "unified/triangle_coarsened", 1);

  // Perform adaptive refinement
  // std::cout << "3. Performing adaptive refinement...\n";
  // mra.refinement_new (min_level, max_level);
  //
  // num_elements = t8_forest_get_global_num_leaf_elements (mra.get_forest ());
  // std::cout << "\n   After refinement:\n";
  // std::cout << "   Elements: " << num_elements << "\n";
  // std::cout << "   Total DOF: " << (num_elements * mra.DOF) << "\n\n";
  //
  // // Write refined solution
  // write_vtk_output (mra, "unified/triangle_refined", 2);

  // Cleanup
  mra.cleanup ();
  t8_cmesh_destroy (&cmesh);
  t8_scheme_unref (const_cast<t8_scheme **> (&scheme));

  std::cout << "✓ Triangle example completed!\n\n";
}

//=============================================================================
// Example 2: Quad MRA with Adaptation and Plotting
//=============================================================================

void
example_quad_adaptive_with_plotting ()
{
  std::cout << "\n";
  std::cout << "════════════════════════════════════════════════════════════\n";
  std::cout << "  Quad MRA: Adaptive Refinement with VTK Output\n";
  std::cout << "════════════════════════════════════════════════════════════\n\n";

  // MRA parameters
  constexpr int U = 1;
  constexpr int P = 3;
  const int min_level = 0;
  const int max_level = 6;
  const double c_thresh = 1.0;
  const int gamma = 1;
  const int num_quad_points_1d = 4;
  const bool balanced = false;

  // Create multiscale object
  t8_mra::multiscale<T8_ECLASS_QUAD, U, P> mra (max_level, c_thresh, gamma, num_quad_points_1d, balanced,
                                                sc_MPI_COMM_WORLD);

  // Create mesh
  t8_cmesh_t cmesh = t8_cmesh_new_hypercube (T8_ECLASS_QUAD, sc_MPI_COMM_WORLD, 0, 0, 0);
  auto *scheme = t8_scheme_new_default ();

  // Initialize with Gaussian bump
  std::cout << "1. Initializing uniform mesh at level " << max_level << "...\n";
  auto func = gaussian_bump<U> ();
  mra.initialize_data (cmesh, scheme, max_level, func);

  auto num_elements = t8_forest_get_global_num_leaf_elements (mra.get_forest ());
  std::cout << "   Elements: " << num_elements << "\n";
  std::cout << "   Total DOF: " << (num_elements * mra.DOF) << "\n\n";

  // Write initial uniform solution
  write_vtk_output (mra, "unified/quad_uniform", 0);

  // Perform adaptive coarsening
  std::cout << "2. Performing adaptive coarsening...\n";
  mra.coarsening_new (min_level, max_level);

  num_elements = t8_forest_get_global_num_leaf_elements (mra.get_forest ());
  std::cout << "\n   After coarsening:\n";
  std::cout << "   Elements: " << num_elements << "\n";
  std::cout << "   Total DOF: " << (num_elements * mra.DOF) << "\n\n";

  // Write coarsened solution
  write_vtk_output (mra, "unified/quad_coarsened", 1);

  // Perform adaptive refinement
  // std::cout << "3. Performing adaptive refinement...\n";
  // mra.refinement_new (min_level, max_level);
  //
  // num_elements = t8_forest_get_global_num_leaf_elements (mra.get_forest ());
  // std::cout << "\n   After refinement:\n";
  // std::cout << "   Elements: " << num_elements << "\n";
  // std::cout << "   Total DOF: " << (num_elements * mra.DOF) << "\n\n";
  //
  // // Write refined solution
  // write_vtk_output (mra, "unified/quad_refined", 2);

  // Cleanup
  mra.cleanup ();
  t8_cmesh_destroy (&cmesh);
  t8_scheme_unref (const_cast<t8_scheme **> (&scheme));

  std::cout << "✓ Quad example completed!\n\n";
}

//=============================================================================
// Example 3: Full Adaptation Cycle with Multiple Iterations
//=============================================================================

void
example_full_adaptation_cycle ()
{
  std::cout << "\n";
  std::cout << "════════════════════════════════════════════════════════════\n";
  std::cout << "  Full Adaptation Cycle: Multiple Iterations\n";
  std::cout << "════════════════════════════════════════════════════════════\n\n";

  // MRA parameters
  constexpr int U = 1;
  constexpr int P = 3;
  const int min_level = 2;
  const int max_level = 6;
  const double c_thresh = 0.03;
  const int gamma = 1;
  const int dunavant_rule = 5;
  const bool balanced = false;

  // Create multiscale object
  t8_mra::multiscale<T8_ECLASS_TRIANGLE, U, P> mra (max_level, c_thresh, gamma, dunavant_rule, balanced,
                                                    sc_MPI_COMM_WORLD);

  // Create mesh
  t8_cmesh_t cmesh = t8_cmesh_new_hypercube (T8_ECLASS_TRIANGLE, sc_MPI_COMM_WORLD, 0, 0, 0);
  auto *scheme = t8_scheme_new_default ();

  // Initialize with sine wave
  std::cout << "Initializing with sine wave function...\n";
  auto func = sine_wave<U> ();
  mra.initialize_data (cmesh, scheme, max_level, func);

  auto num_elements = t8_forest_get_global_num_leaf_elements (mra.get_forest ());
  std::cout << "Initial elements: " << num_elements << "\n\n";

  write_vtk_output (mra, "unified/cycle_initial", 0);

  // Perform multiple adaptation cycles
  const int num_cycles = 3;
  for (int cycle = 1; cycle <= num_cycles; ++cycle) {
    std::cout << "─────────────────────────────────────\n";
    std::cout << "Adaptation Cycle " << cycle << ":\n";
    std::cout << "─────────────────────────────────────\n";

    // Full adaptation (coarsening + refinement)
    mra.adapt (min_level, max_level);

    num_elements = t8_forest_get_global_num_leaf_elements (mra.get_forest ());
    std::cout << "Elements after cycle " << cycle << ": " << num_elements << "\n\n";

    write_vtk_output (mra, "unified/cycle_adapted", cycle);
  }

  // Cleanup
  mra.cleanup ();
  t8_cmesh_destroy (&cmesh);
  t8_scheme_unref (const_cast<t8_scheme **> (&scheme));

  std::cout << "✓ Full adaptation cycle completed!\n\n";
}

//=============================================================================
// Example 4: Comparison - Different Test Functions
//=============================================================================

void
example_comparison_test_functions ()
{
  std::cout << "\n";
  std::cout << "════════════════════════════════════════════════════════════\n";
  std::cout << "  Comparison: Different Test Functions\n";
  std::cout << "════════════════════════════════════════════════════════════\n\n";

  // MRA parameters
  constexpr int U = 1;
  constexpr int P = 3;
  const int min_level = 2;
  const int max_level = 5;
  const double c_thresh = 0.05;
  const int gamma = 1;
  const int num_quad_points_1d = 4;
  const bool balanced = false;

  auto *scheme = t8_scheme_new_default ();

  // Test functions
  struct TestCase
  {
    std::string name;
    std::function<std::array<double, U> (double, double)> func;
  };

  std::vector<TestCase> test_cases = {
    { "gaussian", gaussian_bump<U> () },
    { "sine", sine_wave<U> () },
    { "step", step_function<U> () },
  };

  for (const auto &test : test_cases) {
    std::cout << "Testing function: " << test.name << "\n";

    // Create quad MRA
    t8_mra::multiscale<T8_ECLASS_QUAD, U, P> mra (max_level, c_thresh, gamma, num_quad_points_1d, balanced,
                                                  sc_MPI_COMM_WORLD);

    // Create mesh
    t8_cmesh_t cmesh = t8_cmesh_new_hypercube (T8_ECLASS_QUAD, sc_MPI_COMM_WORLD, 0, 0, 0);

    // Initialize
    mra.initialize_data (cmesh, scheme, max_level, test.func);

    const auto num_init = t8_forest_get_global_num_leaf_elements (mra.get_forest ());
    write_vtk_output (mra, "unified/compare_" + test.name + "_uniform", 0);

    // Adapt
    mra.adapt (min_level, max_level);

    const auto num_adapted = t8_forest_get_global_num_leaf_elements (mra.get_forest ());
    std::cout << "  Elements: " << num_init << " -> " << num_adapted;
    std::cout << " (compression: " << (100.0 * (1.0 - (double) num_adapted / num_init)) << "%)\n";

    write_vtk_output (mra, "unified/compare_" + test.name + "_adapted", 1);

    // Cleanup
    mra.cleanup ();
    t8_cmesh_destroy (&cmesh);
  }

  t8_scheme_unref (const_cast<t8_scheme **> (&scheme));

  std::cout << "\n✓ Comparison example completed!\n\n";
}

//=============================================================================
// Main
//=============================================================================

int
main (int argc, char **argv)
{
  int mpiret;
  sc_MPI_Comm comm;

  mpiret = sc_MPI_Init (&argc, &argv);
  SC_CHECK_MPI (mpiret);

  comm = sc_MPI_COMM_WORLD;
  sc_init (comm, 1, 1, nullptr, SC_LP_ESSENTIAL);
  t8_init (SC_LP_PRODUCTION);

  std::cout << "\n";
  std::cout << "╔════════════════════════════════════════════════════════════╗\n";
  std::cout << "║     Unified MRA: Adaptive Refinement with VTK Output      ║\n";
  std::cout << "║                                                            ║\n";
  std::cout << "║  This example demonstrates:                               ║\n";
  std::cout << "║  • Triangle and quad MRA with unified interface           ║\n";
  std::cout << "║  • Adaptive refinement and coarsening                     ║\n";
  std::cout << "║  • VTK output generation for visualization                ║\n";
  std::cout << "║  • Comparison of different test functions                 ║\n";
  std::cout << "╚════════════════════════════════════════════════════════════╝\n";

  // Run examples
  example_triangle_adaptive_with_plotting ();
  example_quad_adaptive_with_plotting ();
  // example_full_adaptation_cycle ();
  // example_comparison_test_functions ();

  std::cout << "════════════════════════════════════════════════════════════\n";
  std::cout << "✓ All plotting examples completed!\n";
  std::cout << "════════════════════════════════════════════════════════════\n\n";

  std::cout << "Generated VTK files:\n";
  std::cout << "  Triangle:\n";
  std::cout << "    - triangle_uniform_step0.vtu\n";
  std::cout << "    - triangle_coarsened_step1.vtu\n";
  std::cout << "    - triangle_refined_step2.vtu\n\n";
  std::cout << "  Quad:\n";
  std::cout << "    - quad_uniform_step0.vtu\n";
  std::cout << "    - quad_coarsened_step1.vtu\n";
  std::cout << "    - quad_refined_step2.vtu\n\n";
  std::cout << "  Cycles:\n";
  std::cout << "    - cycle_initial_step0.vtu\n";
  std::cout << "    - cycle_adapted_step1.vtu ... step3.vtu\n\n";
  std::cout << "  Comparison:\n";
  std::cout << "    - compare_gaussian_uniform_step0.vtu\n";
  std::cout << "    - compare_gaussian_adapted_step1.vtu\n";
  std::cout << "    - compare_sine_*.vtu\n";
  std::cout << "    - compare_step_*.vtu\n\n";

  std::cout << "Open these files in ParaView to visualize the results!\n\n";

  sc_finalize ();
  return 0;
}

#endif  // T8_ENABLE_MRA
