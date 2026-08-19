/**
 * @file t8_mra_example.cxx
 * @brief MRA tutorial: adaptive coarsening/refinement with VTK output
 *
 * Demonstrates:
 * 1. Full adaptation cycle (top-down): initialize -> coarsen -> refine -> coarsen
 * 2. Bottom-up initialization: adaptive grid without building the uniform grid
 * 3. Custom adaptation criteria
 * 4. Triangle vs quad comparison on the same data
 * 5. 3D (hex) adaptation
 * 6. Two state variables (U = 2) with different jump locations
 */

#include "t8.h"
#if T8_ENABLE_MRA

#include "t8_mra/t8_mra.hxx"
#include "t8_cmesh/t8_cmesh.h"
#include "t8_cmesh/t8_cmesh_examples.h"

#include <cmath>
#include <filesystem>
#include <iostream>
#include <string>

#include <sc_refcount.h>
#include <sc_flops.h>
#include <sc_statistics.h>
#include <sc_options.h>

#include <t8_vtk/t8_vtk_writer.h>

//=============================================================================
// Output Helpers
//=============================================================================


 const bool vtk_activated = true;

void timing_start(sc_flopinfo_t& fi, sc_flopinfo_t& snapshot)
{
  sc_flops_start (&fi);
  sc_flops_snap (&fi, &snapshot);
}

void timing_end(const std::string& name, sc_flopinfo_t& fi, sc_flopinfo_t& snapshot, sc_statinfo_t* stats)
{
  sc_flops_shot (&fi, &snapshot);
  sc_stats_set1 (&stats[0], snapshot.iwtime, name.c_str());
  sc_stats_compute (sc_MPI_COMM_WORLD, 1, stats);
}

void timing_print_all(sc_statinfo_t* stats, int nstats)
{
  for(int i=0; i<nstats; i++)
  {
    sc_stats_print (-1, SC_LP_INFO, 1, &(stats[i]), 1, 1);
  }
}

/**
 * @brief Rank-0-only stdout
 *
 * Every printed value in the examples is global (global element counts,
 * parameters), so the other ranks would only duplicate the lines.
 */
struct root_ostream
{
  bool root;

  template <typename T>
  const root_ostream &
  operator<< (const T &v) const
  {
    if (root)
      std::cout << v;
    return *this;
  }
};

static root_ostream
root_out ()
{
  int rank;
  sc_MPI_Comm_rank (sc_MPI_COMM_WORLD, &rank);
  return { rank == 0 };
}

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

/**
 * @brief Two polynomials separated by a jump along a quarter circle
 */
template <int U>
auto
quarter_circle ()
{
  return [] (double x, double y) -> std::array<double, U> {
    double r = x * x + y * y;
    return { (r < 0.25) ? (x * y + x + 3.) : (x * x * y - 2. * x * y * y + 3. * x) };
  };
}

/**
 * @brief Two components with quarter-circle jumps around the bottom-left
 * (u0) and top-right (u1) corner
 */
template <int U>
auto
two_quarter_circles ()
{
  return [] (double x, double y) -> std::array<double, U> {
    const double xm = 1. - x;
    const double ym = 1. - y;
    const double r0 = x * x + y * y;
    const double r1 = xm * xm + ym * ym;
    return { (r0 < 0.25) ? (x * y + x + 3.) : (x * x * y - 2. * x * y * y + 3. * x),
             (r1 < 0.25) ? (xm * ym + xm + 3.) : (xm * xm * ym - 2. * xm * ym * ym + 3. * xm) };
  };
}

/**
 * @brief 3D Gaussian bump function
 */
template <int U>
auto
gaussian_bump_3d ()
{
  return [] (double x, double y, double z) -> std::array<double, U> {
    const double cx = 0.5, cy = 0.5, cz = 0.5;
    const double sigma = 0.15;
    const double r2 = (x - cx) * (x - cx) + (y - cy) * (y - cy) + (z - cz) * (z - cz);
    return { std::exp (-r2 / (2 * sigma * sigma)) };
  };
}

//=============================================================================
// Helpers
//=============================================================================

/**
 * @brief Write high-order Lagrange VTK output (polynomial degree P-1)
 */
template <typename MRA>
void
write_vtk_output (MRA &mra, const std::string &filename)
{

  if(not vtk_activated) return;

  root_out () << "  Writing VTK: " << filename << ".vtu\n";

  t8_mra::write_forest_lagrange_vtk (mra, filename.c_str (), MRA::P_DIM - 1);
}

/**
 * @brief Print element and DOF count of the current grid
 */
template <typename MRA>
t8_gloidx_t
print_grid_stats (MRA &mra, const std::string &label)
{
  const auto num_elements = t8_forest_get_global_num_leaf_elements (mra.get_forest ());
  root_out () << "  " << label << ": " << num_elements << " elements, " << (num_elements * MRA::DOF * MRA::U_DIM)
              << " DOF\n";
  return num_elements;
}

//=============================================================================
// Example 1: Full Adaptation Cycle (top-down)
//=============================================================================

/**
 * This example functions runs through one full top-down adaptation cycle,
 * by going through the following steps:
 *
 *  Step 0: Initialize on the uniform max_level grid.
 *  Step 1: Coarsen away the non-significant details.
 *  Step 2: Refine via Harten's prediction (grading band at the jump).
 *  Step 3: Coarsen again: the zero-detail children created by Harten's
 *          refinement carry no information, so the resulting grid is the same
 *          as for step 1.
 *  Step 4: Establish a 2:1 balance of the grid.
 */
void
example_adaptation_cycle ()
{
  root_out () << "\n=== 1. Triangle: full adaptation cycle (top-down) ===\n";

  // Define parameters used in this example.
  constexpr int U = 1;
  constexpr int P = 3;
  const int min_level = 0;
  const int max_level = 4;
  const double c_thresh = 1.0;

  int mpisize, mpirank;
  sc_MPI_Comm_size (sc_MPI_COMM_WORLD, &mpisize);
  sc_MPI_Comm_rank (sc_MPI_COMM_WORLD, &mpirank);
  int sqrt_mpisize = static_cast<int>(std::sqrt(mpisize));

  // t8_debugf("sqrt_mpisize = %i", sqrt_mpisize);

  // Set up instance of MRA class with the given (template) parameters.
  t8_mra::multiscale<T8_ECLASS_TRIANGLE, U, P> mra (max_level, sc_MPI_COMM_WORLD);

  // Define cmesh and scheme.
  const double boundary_coords[24] = { 0, 0, 0, 1, 0, 0, 0, 1, 0, 1, 1, 0, 0, 0, 1, 1, 0, 1, 0, 1, 1, 1, 1, 1 };
  // const double boundary_coords[8] = { 0, 0, 1, 0, 0, 1, 1, 1};
  // t8_cmesh_t cmesh = t8_cmesh_new_hypercube (T8_ECLASS_TRIANGLE, sc_MPI_COMM_WORLD, 0, 0, 0);

  // const  t8_cmesh_get_num_trees (cmesh)

  const int n_x = 2;
  const int n_y = 2;
  const int num_global_trees = n_x * n_y * 2;


  // const t8_gloidx_t local_tree_offset = mpirank * num_global_trees / mpisize;

  const t8_gloidx_t local_tree_offset = 0;

  t8_productionf("Local tree offset:%i\n", local_tree_offset );


  root_out () << "\nCreate cmesh...\n";

  // t8_cmesh_t cmesh = t8_cmesh_new_hypercube_pad (T8_ECLASS_TRIANGLE, sc_MPI_COMM_WORLD, boundary_coords, 10, 10, 1, false);
  t8_cmesh_t cmesh = t8_cmesh_new_hypercube_pad_ext(T8_ECLASS_TRIANGLE, sc_MPI_COMM_WORLD, boundary_coords, n_x, n_y, 1, 0, 0, 0, false, false, local_tree_offset);

  root_out () << "\n... done!\n";

  if(vtk_activated) t8_cmesh_vtk_write_file(cmesh, "mra_output/mra_cmesh");

  auto *scheme = t8_scheme_new_default ();

  sc_flopinfo_t fi, snapshot;
  sc_statinfo_t stats[5];

  timing_start(fi, snapshot);
  // Step 0: Initialize on the uniform max_level grid.
  mra.initialize_data (cmesh, scheme, max_level, quarter_circle<U> ());

  // mra.initialize_data_adaptive (cmesh, scheme, max_level, quarter_circle<U> (),
  //                               t8_mra::hard_thresholding { .c_thresh = c_thresh });


  timing_end("MRA Init", fi, snapshot, &(stats[0]));

  print_grid_stats (mra, "Uniform level " + std::to_string (max_level));
  write_vtk_output (mra, "mra_output/01_cycle_step0_uniform");
  // write_vtk_output (mra, "mra_output/01_cycle_step0_bottom_up");

  mra.repartition();

  write_vtk_output (mra, "mra_output/01_cycle_step0_uniform_repartitioned");

  // Step 1: Coarsen away the non-significant details.
  timing_start(fi, snapshot);
  mra.coarsen (min_level, max_level, t8_mra::hard_thresholding { .c_thresh = c_thresh });
  timing_end("MRA Coarsen", fi, snapshot, &(stats[1]));
  const auto num_coarse = print_grid_stats (mra, "After coarsening");
  write_vtk_output (mra, "mra_output/01_cycle_step1_coarsened");

  // // Step 2: Refine via Harten's prediction (grading band at the jump).
  // timing_start(fi, snapshot);
  // mra.refine (min_level, max_level, t8_mra::harten_prediction { .c_thresh = c_thresh });
  // timing_end("MRA Refine", fi, snapshot, &(stats[2]));
  // print_grid_stats (mra, "After refinement");
  // write_vtk_output (mra, "mra_output/01_cycle_step2_refined");

  // // Step 3: Coarsen again.
  // timing_start(fi, snapshot);
  // mra.coarsen (min_level, max_level, t8_mra::hard_thresholding { .c_thresh = c_thresh });
  // timing_end("MRA Coarsen 2", fi, snapshot, &(stats[3]));
  // const auto num_recoarse = print_grid_stats (mra, "After second coarsening");
  // write_vtk_output (mra, "mra_output/01_cycle_step3_coarsened");

  // // Log output: Compare number of elements for the two coarse meshes (step 1 and step 3).
  // //             They should match since the Harten's prediction only added zero-detail children.
  // root_out () << "  Round-trip: " << num_coarse << " -> " << num_recoarse
  //             << (num_coarse == num_recoarse ? " (exact)\n" : "\n");

  // // Step 4: Establish a 2:1 balance of the grid.
  // timing_start(fi, snapshot);
  // mra.balance ();
  // timing_end("MRA Balance", fi, snapshot, &(stats[4]));
  // const auto num_balanced = print_grid_stats (mra, "After balancing");
  // write_vtk_output (mra, "mra_output/01_cycle_step4_balance");

  // timing_print_all(stats,5);

  // Clean up memory.
  mra.cleanup ();
}

//=============================================================================
// Example 2: Bottom-Up Initialization
//=============================================================================

/**
 * Build the adaptive grid directly from the initial data: project on level 1,
 * then refine level by level only where the details are significant. The
 * uniform max_level grid of example 1 is never built.
 */
void
example_bottom_up ()
{
  root_out () << "\n=== 2. Triangle: bottom-up initialization ===\n";

  // Define the parameters used in this example.
  constexpr int U = 1;
  constexpr int P = 3;
  const int max_level = 7;
  const double c_thresh = 1.0;

  // Set up instance of MRA class with the given (template) parameters.
  t8_mra::multiscale<T8_ECLASS_TRIANGLE, U, P> mra (max_level, sc_MPI_COMM_WORLD);

  // Define cmesh and scheme.
  t8_cmesh_t cmesh = t8_cmesh_new_hypercube (T8_ECLASS_TRIANGLE, sc_MPI_COMM_WORLD, 0, 0, 0);
  auto *scheme = t8_scheme_new_default ();

  // Initialize MRA grid through a level-by-level adaptive refimenet, i.e., with a bottom-up strategy
  // that does not use the uniform max-level grid.
  mra.initialize_data_adaptive (cmesh, scheme, max_level, quarter_circle<U> (),
                                t8_mra::hard_thresholding { .c_thresh = c_thresh });

  // Print out statistics of the resulting adaptive grid.
  const auto num_adaptive = print_grid_stats (mra, "Adaptive grid");
  const auto num_trees = t8_forest_get_num_global_trees (mra.get_forest ());
  const auto num_uniform = num_trees * static_cast<t8_gloidx_t> (std::pow (4, max_level));
  root_out () << "  Uniform level " << max_level << " grid (never built): " << num_uniform << " elements\n";
  root_out () << "  Compression: " << (100.0 * (1.0 - static_cast<double> (num_adaptive) / num_uniform)) << " %\n";

  // Write VTK output.
  write_vtk_output (mra, "mra_output/02_bottom_up");

  // Clean up memory.
  mra.cleanup ();
}

//=============================================================================
// Example 3: Custom Adaptation Criterion
//=============================================================================

/**
 * @brief Hard thresholding with an enforced minimum refinement level
 *
 * Any type satisfying the coarsening_criterion concept can be passed to
 * coarsen(). This one composes the default thresholding with a level floor:
 * families below floor_level are always kept refined.
 */
struct thresholding_with_floor
{
  t8_mra::hard_thresholding thresholding;
  unsigned int floor_level = 4;

  /**
   * \brief Prepare the threshold criterion.
   *
   * \tparam MRA the specific MRA type used
   * \param [in] mra an instance of MRA
   */
  template <typename MRA>
  void
  prepare (MRA &mra)
  {
    thresholding.prepare (mra);
  }

  /**
   * \brief This function is used to decide whether a given levelmultiindex is significant or not.
   *
   *  In this example, it combines a hard thresholding with a floor level.
   *
   * \tparam MRA          the specific MRA type used.
   * \param [in] mra  an instance of MRA
   * \param [in] lmi  a levelmultiindex
   * \return True if the \a lmi is considered significant, false if not.
   */
  template <typename MRA>
  bool
  significant (MRA &mra, const typename MRA::levelmultiindex &lmi)
  {
    return lmi.level () < floor_level || thresholding.significant (mra, lmi);
  }
};

/**
 * \brief This example demonstrates the usage of a custom adaptation criterion.
 *
 *  Specifically, it uses a custom criterion that combines extends the hard thresholding
 *  provided by t8_mra by a floor level, i.e., a minimum level below which we do not coarsen.
 */
void
example_custom_criterion ()
{
  root_out () << "\n=== 3. Quad: custom coarsening criterion (level floor) ===\n";

  // Define parameters used in this example.
  constexpr int U = 1;
  constexpr int P = 3;
  const int min_level = 0;
  const int max_level = 6;
  const double c_thresh = 1.0;
  const unsigned int floor_level = 4;

  // Set up two MRA instances: One for a plain hard thresholding and one for the thresholding with floor level.
  // NOTE: The difference will be introduced below by passing different criteria to the coarsen function.
  t8_mra::multiscale<T8_ECLASS_QUAD, U, P> mra_plain (max_level, sc_MPI_COMM_WORLD);
  t8_mra::multiscale<T8_ECLASS_QUAD, U, P> mra_floor (max_level, sc_MPI_COMM_WORLD);

  // Define a cmesh and a shceme.
  t8_cmesh_t cmesh = t8_cmesh_new_hypercube (T8_ECLASS_QUAD, sc_MPI_COMM_WORLD, 0, 0, 0);
  auto *scheme = t8_scheme_new_default ();
  // Increase reference counter for cmesh and t8scheme to allow both MRA instances to use them.
  t8_cmesh_ref (cmesh);
  t8_scheme_ref (const_cast<t8_scheme *> (scheme));

  // Use a Gaussian bump as example.
  auto func = gaussian_bump<U> ();

  // Part 1: Plain hard thresholding
  //         Initialize the grid at max_level, coarsen insignificant parts, and write (vtk) output.
  mra_plain.initialize_data (cmesh, scheme, max_level, func);
  mra_plain.coarsen (min_level, max_level, t8_mra::hard_thresholding { .c_thresh = c_thresh });
  print_grid_stats (mra_plain, "Plain thresholding");
  write_vtk_output (mra_plain, "mra_output/03_criterion_plain");

  // Part 2: Hard thresholding combined with floor level
  //         Initialize the grid at max_level, coarsen insignificant parts, and write (vtk) output.
  mra_floor.initialize_data (cmesh, scheme, max_level, func);
  mra_floor.coarsen (min_level, max_level, thresholding_with_floor { { .c_thresh = c_thresh }, floor_level });
  print_grid_stats (mra_floor, "With level floor " + std::to_string (floor_level));
  write_vtk_output (mra_floor, "mra_output/03_criterion_floor");

  // Memory clean up.
  mra_plain.cleanup ();
  mra_floor.cleanup ();
}

//=============================================================================
// Example 4: Triangle vs Quad
//=============================================================================

/**
 * \brief Perform one example run (initialize + coarsen) for a specific shape.
 *
 * \tparam Shape              The element shape, given as t8_eclass.
 *
 * \param [in] name       The name used for output.
 * \param [in] func       The function considered in the MRA.
 * \param [in] max_level  The maximum level.
 * \param [in] c_thresh   The threshold value.
 *
 * \return The number of elements in the resulting adaptive grid.
 */
template <t8_eclass Shape>
t8_gloidx_t
run_shape (const std::string &name, auto &&func, int max_level, double c_thresh)
{
  // Define template parameters.
  constexpr int U = 1;
  constexpr int P = 3;

  // Set up MRA instance.
  t8_mra::multiscale<Shape, U, P> mra (max_level, sc_MPI_COMM_WORLD);

  // Define a cmesh and a scheme.
  t8_cmesh_t cmesh = t8_cmesh_new_hypercube (Shape, sc_MPI_COMM_WORLD, 0, 0, 0);
  auto *scheme = t8_scheme_new_default ();

  // Initialize the MRA grid on max level.
  mra.initialize_data (cmesh, scheme, max_level, func);

  // Store number of elements in uniform max-level grid.
  const auto num_uniform = t8_forest_get_global_num_leaf_elements (mra.get_forest ());

  // Coarsen insignificant elements.
  mra.coarsen (0, max_level, t8_mra::hard_thresholding { .c_thresh = c_thresh });

  // Store number of elements in resulting grid.
  const auto num_adapted = print_grid_stats (mra, name + " adapted");

  // Log output.
  root_out () << "    " << num_uniform << " -> " << num_adapted
              << " elements (compression: " << (100.0 * (1.0 - static_cast<double> (num_adapted) / num_uniform))
              << " %)\n";
  write_vtk_output (mra, "mra_output/04_compare_" + name);

  // Clean up memory.
  mra.cleanup ();

  // Return the number of elements.
  return num_adapted;
}

/**
 * Same data, same threshold, same domain: compare how triangle and quad
 * grids adapt. The triangle hypercube consists of two base trees, the quad
 * of one; DOF per element differ (P(P+1)/2 vs P^2), so compare total DOF.
 */
void
example_triangle_vs_quad ()
{
  root_out () << "\n=== 4. Triangle vs quad on the same data ===\n";

  // Set example parameters.
  const int max_level = 6;
  const double c_thresh = 1.0;
  auto func = quarter_circle<1> ();

  // Run the same example twice: For triangular elements and for quadrilaterals.
  run_shape<T8_ECLASS_TRIANGLE> ("triangle", func, max_level, c_thresh);
  run_shape<T8_ECLASS_QUAD> ("quad", func, max_level, c_thresh);
}

//=============================================================================
// Example 5: Hex (3D)
//=============================================================================

/**
 * \brief This example demonstrates the MRA functionality for hexahedral elements.
 *
 *  Specifically, the MRA mesh is first created on the maximum level (step 0)
 *  and then coarsened in step 1, before step 2 applies a refinement based on
 *  Harten's prediction.
 */
void
example_hex_3d ()
{
  root_out () << "\n=== 5. Hex: 3D adaptation ===\n";

  // Define parameters used in this example.
  constexpr int U = 1;
  constexpr int P = 3;
  const int min_level = 0;
  const int max_level = 4;
  const double c_thresh = 1.0;

  // Set up MRA instance.
  t8_mra::multiscale<T8_ECLASS_HEX, U, P> mra (max_level, sc_MPI_COMM_WORLD);

  // Define mesh and scheme.
  t8_cmesh_t cmesh = t8_cmesh_new_hypercube (T8_ECLASS_HEX, sc_MPI_COMM_WORLD, 0, 0, 0);
  auto *scheme = t8_scheme_new_default ();

  // Step 0; Initialize the MRA grid on maximum level.
  mra.initialize_data (cmesh, scheme, max_level, gaussian_bump_3d<U> ());
  print_grid_stats (mra, "Uniform level " + std::to_string (max_level));
  write_vtk_output (mra, "mra_output/05_hex_step0_uniform");

  // Step 1: Coarsen the grid.
  mra.coarsen (min_level, max_level, t8_mra::hard_thresholding { .c_thresh = c_thresh });
  print_grid_stats (mra, "After coarsening");
  write_vtk_output (mra, "mra_output/05_hex_step1_coarsened");

  // Step 2: Refine by Harten's prediction (around features).
  mra.refine (min_level, max_level, t8_mra::harten_prediction { .c_thresh = c_thresh });
  print_grid_stats (mra, "After refinement");
  write_vtk_output (mra, "mra_output/05_hex_step2_refined");

  // Clean up memory.
  mra.cleanup ();

  // Print hint to user.
  root_out () << "  Use ParaView 'Clip' / 'Slice' filters to see the internal structure.\n";
}

//=============================================================================
// Example 6: Two State Variables
//=============================================================================

/**
 * U = 2: each component carries its own quarter-circle jump (bottom-left vs
 * top-right corner). Significance is the maximum over the components, so
 * the grid refines along both arcs.
 */
void
example_two_components ()
{
  root_out () << "\n=== 6. Triangle: two state variables ===\n";

  // Set the example's parameters
  constexpr int U = 2;
  constexpr int P = 3;
  const int min_level = 0;
  const int max_level = 7;
  const double c_thresh = 0.2;

  // Create the MRA instance.
  t8_mra::multiscale<T8_ECLASS_TRIANGLE, U, P> mra (max_level, sc_MPI_COMM_WORLD);

  // Define a cmesh and a scheme.
  t8_cmesh_t cmesh = t8_cmesh_new_hypercube (T8_ECLASS_TRIANGLE, sc_MPI_COMM_WORLD, 0, 0, 0);
  auto *scheme = t8_scheme_new_default ();

  // Step 0: Initialize the MRA grid adaptively with a bottom-up strategy.
  mra.initialize_data_adaptive (cmesh, scheme, max_level, two_quarter_circles<U> ());
  print_grid_stats (mra, "Uniform level " + std::to_string (max_level));
  write_vtk_output (mra, "mra_output/06_two_components_step0_initial");

  // Step 1: Coarsen the grid.
  mra.coarsen (min_level, max_level, t8_mra::hard_thresholding { .c_thresh = c_thresh });
  print_grid_stats (mra, "After coarsening");
  write_vtk_output (mra, "mra_output/06_two_components_step1_coarsened");

  // Step 2: Refine using Harten's prediction.
  mra.refine (min_level, max_level, t8_mra::harten_prediction { .c_thresh = c_thresh });
  print_grid_stats (mra, "After refinement");
  write_vtk_output (mra, "mra_output/06_two_components_step2_refined");

  // Clean up memory.
  mra.cleanup ();

  // Print hint to user.
  root_out () << "  Color by u0 / u1 in ParaView: the grid follows both jumps.\n";
}

//=============================================================================
// Main
//=============================================================================

/**
 * \brief Main routine running all six examples covered in this file.
 *
 * The examples are:
 *  1. Full adaptation cycle (top-down): initialize -> coarsen -> refine -> coarsen
 *  2. Bottom-up initialization: adaptive grid without building the uniform grid
 *  3. Custom adaptation criteria
 *  4. Triangle vs quad comparison on the same data
 *  5. 3D (hex) adaptation
 *  6. Two state variables (U = 2) with different jump locations
 *
 * \param [in] argc Standard argc (not used, only passed to sc_MPI_Init)
 * \param [in] argv Standard argv (not used, only passed to sc_MPI_Init)
 * \return 0
 */
int
main (int argc, char **argv)
{
  // Initialize SC library for MPI handling.
  int mpiret = sc_MPI_Init (&argc, &argv);
  SC_CHECK_MPI (mpiret);
  sc_init (sc_MPI_COMM_WORLD, 1, 1, nullptr, SC_LP_ESSENTIAL);

  // Initialize t8code.
  t8_init (SC_LP_DEBUG);

  // Rank 0 creates the output directory; everyone writes into it
  int mpirank;
  sc_MPI_Comm_rank (sc_MPI_COMM_WORLD, &mpirank);
  if (mpirank == 0)
    std::filesystem::create_directory ("mra_output");
  sc_MPI_Barrier (sc_MPI_COMM_WORLD);
  if (!std::filesystem::exists ("mra_output"))
    t8_errorf ("Could not create directory");

  // Example 1: Full adaptation cycle (top-down): initialize -> coarsen -> refine -> coarsen.
  example_adaptation_cycle ();

  // // Example 2: Bottom-up initialization: adaptive grid without building the uniform grid.
  // example_bottom_up ();

  // // Example 3: Custom adaptation criteria.
  // example_custom_criterion ();

  // // Example 4: Triangle vs quad comparison on the same data.
  // example_triangle_vs_quad ();

  // // Example 5: 3D (hex) adaptation.
  // example_hex_3d ();

  // // Example 6: Two state variables (U = 2) with different jump locations.
  // example_two_components ();

  // Log output.
  root_out () << "\nAll examples completed. Output in mra_output/ (open in ParaView).\n";

  // Finalize SC library.
  sc_finalize ();
  mpiret = sc_MPI_Finalize ();
  SC_CHECK_MPI (mpiret);

  return 0;
}

#endif  // T8_ENABLE_MRA
