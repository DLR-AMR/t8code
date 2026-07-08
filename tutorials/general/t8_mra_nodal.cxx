/**
 * @file t8_mra_nodal.cxx
 * @brief MRA tutorial: load nodal DG data off an existing forest
 *
 * A nodal DG code stores its solution as values at nodal points per cell. The
 * MRA stores modal coefficients. initialize_data_nodal bridges the two: given a
 * set of reference nodes and a per-cell nodal-value provider, it converts each
 * leaf to modal coefficients (one Vandermonde solve per component), then the
 * usual MRA machinery (coarsening, detail thresholding) applies.
 *
 * Here the "nodal solver" is faked by sampling an analytic field at each cell's
 * nodes; in practice the provider reads the cell's data from your DG state.
 */

#include "t8.h"
#ifdef T8_ENABLE_MRA

#include "t8_mra/t8_mra.hxx"
#include "t8_cmesh/t8_cmesh.h"
#include "t8_cmesh/t8_cmesh_examples.h"
#include "t8_forest/t8_forest.h"

#include <array>
#include <cmath>
#include <filesystem>
#include <iostream>

static void
root_print (const std::string &s)
{
  int rank;
  sc_MPI_Comm_rank (sc_MPI_COMM_WORLD, &rank);
  if (rank == 0)
    std::cout << s;
}

/// Analytic field the fake nodal solver samples (smooth away from a jump ring).
static double
source_field (double x, double y)
{
  const double r = std::sqrt (x * x + y * y);
  return (r < 0.6) ? (1.0 + x * y) : (0.25 * x - 0.5 * y);
}

int
main (int argc, char **argv)
{
  int mpiret = sc_MPI_Init (&argc, &argv);
  SC_CHECK_MPI (mpiret);
  sc_init (sc_MPI_COMM_WORLD, 1, 1, nullptr, SC_LP_ESSENTIAL);
  t8_init (SC_LP_PRODUCTION);

  int mpirank;
  sc_MPI_Comm_rank (sc_MPI_COMM_WORLD, &mpirank);
  if (mpirank == 0)
    std::filesystem::create_directory ("mra_output");
  sc_MPI_Barrier (sc_MPI_COMM_WORLD);

  root_print ("\n=== MRA: load nodal DG data off a forest ===\n");

  constexpr int U = 1;
  constexpr int P = 3;
  constexpr int DOF = P * P;

  const int max_level = 6;
  const int min_level = 0;

  // Reference nodes: tensor of P equispaced points on [0,1]^2 (node j = jy*P+jx).
  std::array<std::array<double, 2>, DOF> nodes;
  for (int jy = 0; jy < P; ++jy)
    for (int jx = 0; jx < P; ++jx)
      nodes[jy * P + jx] = { jx / static_cast<double> (P - 1), jy / static_cast<double> (P - 1) };

  // Forest the nodal solver lives on (uniform here; the loader accepts any
  // committed, possibly adaptive forest just the same).
  t8_cmesh_t cmesh = t8_cmesh_new_hypercube (T8_ECLASS_QUAD, sc_MPI_COMM_WORLD, 0, 0, 0);
  auto *scheme = t8_scheme_new_default ();
  t8_forest_t forest = t8_forest_new_uniform (cmesh, scheme, max_level, 0, sc_MPI_COMM_WORLD);

  t8_mra::multiscale<T8_ECLASS_QUAD, U, P> mra (max_level, sc_MPI_COMM_WORLD);

  // Provider: sample the field at each node's physical coordinate. A cartesian
  // cell is an axis-aligned box, so vertex 0 is its min corner, vertex 3 its max.
  auto cell_nodal_values = [&] (int tree_idx, const t8_element_t *element) {
    double v0[3], v3[3];
    t8_forest_element_coordinate (forest, tree_idx, element, 0, v0);
    t8_forest_element_coordinate (forest, tree_idx, element, 3, v3);

    std::array<double, U * DOF> nodal;
    for (int j = 0; j < DOF; ++j) {
      const double x = v0[0] + nodes[j][0] * (v3[0] - v0[0]);
      const double y = v0[1] + nodes[j][1] * (v3[1] - v0[1]);
      nodal[j] = source_field (x, y);
    }
    return nodal;
  };

  mra.initialize_data_nodal (forest, nodes, cell_nodal_values);

  const auto num_loaded = t8_forest_get_global_num_leaf_elements (mra.get_forest ());
  root_print ("  Loaded nodal data on " + std::to_string (num_loaded) + " cells\n");

  t8_mra::write_forest_lagrange_vtk (mra, "mra_output/nodal_solution", P - 1);
  root_print ("  Wrote mra_output/nodal_solution (loaded nodal field)\n");

  // MRA now owns modal coeffs: coarsen where details are insignificant.
  mra.coarsen (min_level, max_level, t8_mra::hard_thresholding { .c_thresh = 1e-3 });
  mra.balance ();

  const auto num_coarsened = t8_forest_get_global_num_leaf_elements (mra.get_forest ());
  root_print ("  After MRA coarsening: " + std::to_string (num_coarsened) + " cells\n");

  t8_mra::write_forest_lagrange_vtk (mra, "mra_output/nodal_coarsened", P - 1);
  root_print ("  Wrote mra_output/nodal_coarsened (open in ParaView)\n");

  mra.cleanup ();
  t8_forest_unref (&forest);

  root_print ("\nDone. Output in mra_output/.\n");

  sc_finalize ();
  mpiret = sc_MPI_Finalize ();
  SC_CHECK_MPI (mpiret);
  return 0;
}

#endif  // T8_ENABLE_MRA
