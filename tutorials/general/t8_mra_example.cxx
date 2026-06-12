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
#ifdef T8_ENABLE_MRA

#include "t8_mra/t8_mra.hxx"
#include "t8_cmesh.h"
#include "t8_cmesh/t8_cmesh_examples.h"

#include <cmath>
#include <filesystem>
#include <iostream>
#include <string>

//=============================================================================
// Output Helpers
//=============================================================================

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
 * Initialize on the uniform max_level grid, then coarsen away the
 * non-significant details, refine via Harten's prediction (grading band at
 * the jump) and coarsen again: the zero-detail children created by the
 * refinement carry no information, so the grid returns to the coarsened one.
 */
void
example_adaptation_cycle ()
{
  root_out () << "\n=== 1. Triangle: full adaptation cycle (top-down) ===\n";

  constexpr int U = 1;
  constexpr int P = 3;
  const int min_level = 0;
  const int max_level = 7;
  const double c_thresh = 1.0;

  t8_mra::multiscale<T8_ECLASS_TRIANGLE, U, P> mra (max_level, sc_MPI_COMM_WORLD);

  t8_cmesh_t cmesh = t8_cmesh_new_hypercube (T8_ECLASS_TRIANGLE, sc_MPI_COMM_WORLD, 0, 0, 0);
  auto *scheme = t8_scheme_new_default ();
  /* The forest takes ownership of cmesh and scheme; keep our own references
   * since we destroy/unref them explicitly below. */
  t8_cmesh_ref (cmesh);
  t8_scheme_ref (const_cast<t8_scheme *> (scheme));

  mra.initialize_data (cmesh, scheme, max_level, quarter_circle<U> ());
  print_grid_stats (mra, "Uniform level " + std::to_string (max_level));
  write_vtk_output (mra, "mra_output/01_cycle_step0_uniform");

  mra.coarsen (min_level, max_level, t8_mra::hard_thresholding { .c_thresh = c_thresh });
  const auto num_coarse = print_grid_stats (mra, "After coarsening");
  write_vtk_output (mra, "mra_output/01_cycle_step1_coarsened");

  mra.refine (min_level, max_level, t8_mra::harten_prediction { .c_thresh = c_thresh });
  print_grid_stats (mra, "After refinement");
  write_vtk_output (mra, "mra_output/01_cycle_step2_refined");

  mra.coarsen (min_level, max_level, t8_mra::hard_thresholding { .c_thresh = c_thresh });
  const auto num_recoarse = print_grid_stats (mra, "After second coarsening");
  write_vtk_output (mra, "mra_output/01_cycle_step3_coarsened");

  root_out () << "  Round-trip: " << num_coarse << " -> " << num_recoarse
              << (num_coarse == num_recoarse ? " (exact)\n" : "\n");

  mra.balance ();
  const auto num_balanced = print_grid_stats (mra, "After balancing");
  write_vtk_output (mra, "mra_output/01_cycle_step4_balance");

  mra.cleanup ();
  t8_cmesh_destroy (&cmesh);
  t8_scheme_unref (const_cast<t8_scheme **> (&scheme));
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

  constexpr int U = 1;
  constexpr int P = 3;
  const int max_level = 13;
  const double c_thresh = 1.0;

  t8_mra::multiscale<T8_ECLASS_TRIANGLE, U, P> mra (max_level, sc_MPI_COMM_WORLD);

  t8_cmesh_t cmesh = t8_cmesh_new_hypercube (T8_ECLASS_TRIANGLE, sc_MPI_COMM_WORLD, 0, 0, 0);
  auto *scheme = t8_scheme_new_default ();
  t8_cmesh_ref (cmesh);
  t8_scheme_ref (const_cast<t8_scheme *> (scheme));

  mra.initialize_data_adaptive (cmesh, scheme, max_level, quarter_circle<U> (),
                                t8_mra::hard_thresholding { .c_thresh = c_thresh });

  const auto num_adaptive = print_grid_stats (mra, "Adaptive grid");
  const auto num_trees = t8_forest_get_num_global_trees (mra.get_forest ());
  const auto num_uniform = num_trees * static_cast<t8_gloidx_t> (std::pow (4, max_level));
  root_out () << "  Uniform level " << max_level << " grid (never built): " << num_uniform << " elements\n";
  root_out () << "  Compression: " << (100.0 * (1.0 - static_cast<double> (num_adaptive) / num_uniform)) << " %\n";

  write_vtk_output (mra, "mra_output/02_bottom_up");

  mra.cleanup ();
  t8_cmesh_destroy (&cmesh);
  t8_scheme_unref (const_cast<t8_scheme **> (&scheme));
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

  template <typename MRA>
  void
  prepare (MRA &mra)
  {
    thresholding.prepare (mra);
  }

  template <typename MRA>
  bool
  significant (MRA &mra, const typename MRA::levelmultiindex &lmi)
  {
    return lmi.level () < floor_level || thresholding.significant (mra, lmi);
  }
};

void
example_custom_criterion ()
{
  root_out () << "\n=== 3. Quad: custom coarsening criterion (level floor) ===\n";

  constexpr int U = 1;
  constexpr int P = 3;
  const int min_level = 0;
  const int max_level = 6;
  const double c_thresh = 1.0;
  const unsigned int floor_level = 4;

  t8_mra::multiscale<T8_ECLASS_QUAD, U, P> mra_plain (max_level, sc_MPI_COMM_WORLD);
  t8_mra::multiscale<T8_ECLASS_QUAD, U, P> mra_floor (max_level, sc_MPI_COMM_WORLD);

  t8_cmesh_t cmesh = t8_cmesh_new_hypercube (T8_ECLASS_QUAD, sc_MPI_COMM_WORLD, 0, 0, 0);
  auto *scheme = t8_scheme_new_default ();
  /* Two forests consume one reference each; keep our own on top */
  t8_cmesh_ref (cmesh);
  t8_cmesh_ref (cmesh);
  t8_scheme_ref (const_cast<t8_scheme *> (scheme));
  t8_scheme_ref (const_cast<t8_scheme *> (scheme));

  auto func = gaussian_bump<U> ();

  mra_plain.initialize_data (cmesh, scheme, max_level, func);
  mra_plain.coarsen (min_level, max_level, t8_mra::hard_thresholding { .c_thresh = c_thresh });
  print_grid_stats (mra_plain, "Plain thresholding");
  write_vtk_output (mra_plain, "mra_output/03_criterion_plain");

  mra_floor.initialize_data (cmesh, scheme, max_level, func);
  mra_floor.coarsen (min_level, max_level, thresholding_with_floor { { .c_thresh = c_thresh }, floor_level });
  print_grid_stats (mra_floor, "With level floor " + std::to_string (floor_level));
  write_vtk_output (mra_floor, "mra_output/03_criterion_floor");

  mra_plain.cleanup ();
  mra_floor.cleanup ();
  t8_cmesh_destroy (&cmesh);
  t8_scheme_unref (const_cast<t8_scheme **> (&scheme));
}

//=============================================================================
// Example 4: Triangle vs Quad
//=============================================================================

template <t8_eclass Shape>
t8_gloidx_t
run_shape (const std::string &name, auto &&func, int max_level, double c_thresh)
{
  constexpr int U = 1;
  constexpr int P = 3;

  t8_mra::multiscale<Shape, U, P> mra (max_level, sc_MPI_COMM_WORLD);

  t8_cmesh_t cmesh = t8_cmesh_new_hypercube (Shape, sc_MPI_COMM_WORLD, 0, 0, 0);
  auto *scheme = t8_scheme_new_default ();
  t8_cmesh_ref (cmesh);
  t8_scheme_ref (const_cast<t8_scheme *> (scheme));

  mra.initialize_data (cmesh, scheme, max_level, func);
  const auto num_uniform = t8_forest_get_global_num_leaf_elements (mra.get_forest ());

  mra.coarsen (0, max_level, t8_mra::hard_thresholding { .c_thresh = c_thresh });

  const auto num_adapted = print_grid_stats (mra, name + " adapted");
  root_out () << "    " << num_uniform << " -> " << num_adapted
              << " elements (compression: " << (100.0 * (1.0 - static_cast<double> (num_adapted) / num_uniform))
              << " %)\n";
  write_vtk_output (mra, "mra_output/04_compare_" + name);

  mra.cleanup ();
  t8_cmesh_destroy (&cmesh);
  t8_scheme_unref (const_cast<t8_scheme **> (&scheme));

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

  const int max_level = 6;
  const double c_thresh = 1.0;
  auto func = quarter_circle<1> ();

  run_shape<T8_ECLASS_TRIANGLE> ("triangle", func, max_level, c_thresh);
  run_shape<T8_ECLASS_QUAD> ("quad", func, max_level, c_thresh);
}

//=============================================================================
// Example 5: Hex (3D)
//=============================================================================

void
example_hex_3d ()
{
  root_out () << "\n=== 5. Hex: 3D adaptation ===\n";

  constexpr int U = 1;
  constexpr int P = 3;
  const int min_level = 0;
  const int max_level = 4;
  const double c_thresh = 1.0;

  t8_mra::multiscale<T8_ECLASS_HEX, U, P> mra (max_level, sc_MPI_COMM_WORLD);

  t8_cmesh_t cmesh = t8_cmesh_new_hypercube (T8_ECLASS_HEX, sc_MPI_COMM_WORLD, 0, 0, 0);
  auto *scheme = t8_scheme_new_default ();
  t8_cmesh_ref (cmesh);
  t8_scheme_ref (const_cast<t8_scheme *> (scheme));

  mra.initialize_data (cmesh, scheme, max_level, gaussian_bump_3d<U> ());
  print_grid_stats (mra, "Uniform level " + std::to_string (max_level));
  write_vtk_output (mra, "mra_output/05_hex_step0_uniform");

  mra.coarsen (min_level, max_level, t8_mra::hard_thresholding { .c_thresh = c_thresh });
  print_grid_stats (mra, "After coarsening");
  write_vtk_output (mra, "mra_output/05_hex_step1_coarsened");

  mra.refine (min_level, max_level, t8_mra::harten_prediction { .c_thresh = c_thresh });
  print_grid_stats (mra, "After refinement");
  write_vtk_output (mra, "mra_output/05_hex_step2_refined");

  mra.cleanup ();
  t8_cmesh_destroy (&cmesh);
  t8_scheme_unref (const_cast<t8_scheme **> (&scheme));

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

  constexpr int U = 2;
  constexpr int P = 3;
  const int min_level = 0;
  const int max_level = 7;
  const double c_thresh = 0.2;

  t8_mra::multiscale<T8_ECLASS_TRIANGLE, U, P> mra (max_level, sc_MPI_COMM_WORLD);

  t8_cmesh_t cmesh = t8_cmesh_new_hypercube (T8_ECLASS_TRIANGLE, sc_MPI_COMM_WORLD, 0, 0, 0);
  auto *scheme = t8_scheme_new_default ();
  t8_cmesh_ref (cmesh);
  t8_scheme_ref (const_cast<t8_scheme *> (scheme));

  mra.initialize_data_adaptive (cmesh, scheme, max_level, two_quarter_circles<U> ());
  print_grid_stats (mra, "Uniform level " + std::to_string (max_level));
  write_vtk_output (mra, "mra_output/06_two_components_step0_initial");

  mra.coarsen (min_level, max_level, t8_mra::hard_thresholding { .c_thresh = c_thresh });
  print_grid_stats (mra, "After coarsening");
  write_vtk_output (mra, "mra_output/06_two_components_step1_coarsened");

  mra.refine (min_level, max_level, t8_mra::harten_prediction { .c_thresh = c_thresh });
  print_grid_stats (mra, "After refinement");
  write_vtk_output (mra, "mra_output/06_two_components_step2_refined");

  root_out () << "  Color by u0 / u1 in ParaView: the grid follows both jumps.\n";

  mra.cleanup ();
  t8_cmesh_destroy (&cmesh);
  t8_scheme_unref (const_cast<t8_scheme **> (&scheme));
}

//=============================================================================
// Main
//=============================================================================

int
main (int argc, char **argv)
{
  int mpiret = sc_MPI_Init (&argc, &argv);
  SC_CHECK_MPI (mpiret);

  sc_init (sc_MPI_COMM_WORLD, 1, 1, nullptr, SC_LP_ESSENTIAL);
  t8_init (SC_LP_PRODUCTION);

  // Rank 0 creates the output directory; everyone writes into it
  int mpirank;
  sc_MPI_Comm_rank (sc_MPI_COMM_WORLD, &mpirank);
  if (mpirank == 0)
    std::filesystem::create_directory ("mra_output");
  sc_MPI_Barrier (sc_MPI_COMM_WORLD);
  if (!std::filesystem::exists ("mra_output"))
    t8_errorf ("Could not create directory");

  example_adaptation_cycle ();
  example_bottom_up ();
  example_custom_criterion ();
  example_triangle_vs_quad ();
  example_hex_3d ();
  example_two_components ();

  root_out () << "\nAll examples completed. Output in mra_output/ (open in ParaView).\n";

  sc_finalize ();

  mpiret = sc_MPI_Finalize ();
  SC_CHECK_MPI (mpiret);
  return 0;
}

#endif  // T8_ENABLE_MRA
