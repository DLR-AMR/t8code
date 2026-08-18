/**
 * @file t8_mra_dataset.cxx
 * @brief MRA tutorial: adapt to a given structured dataset
 *
 * A structured_dataset is a func-compatible sampler, so a grid of samples
 * drops straight into initialize_data_adaptive: the MRA builds an adaptive grid
 * that resolves the dataset's features without ever building the uniform grid.
 * Here the source samples come from an analytic field, but any structured field
 * (measurement, simulation output, image) works the same way.
 */

#include "t8.h"
#ifdef T8_ENABLE_MRA

#include "t8_mra/t8_mra.hxx"
#include "t8_cmesh/t8_cmesh.h"
#include "t8_cmesh/t8_cmesh_examples.h"

#include <array>
#include <cmath>
#include <filesystem>
#include <iostream>
#include <vector>

static void
root_print (const std::string &s)
{
  int rank;
  sc_MPI_Comm_rank (sc_MPI_COMM_WORLD, &rank);
  if (rank == 0)
    std::cout << s;
}

/// Sample an analytic field on an N x N grid over [0,1]^2 to build a dataset.
/// Value layout matches structured_dataset: axis 0 (x) slowest, axis 1 (y)
/// fastest, one component per node.
static t8_mra::structured_dataset<2, 1>
build_source_dataset (std::size_t n)
{
  std::vector<double> values (n * n);
  const double h = 1.0 / static_cast<double> (n - 1);

  for (std::size_t i = 0; i < n; ++i)
    for (std::size_t j = 0; j < n; ++j) {
      const double x = i * h;
      const double y = j * h;
      const double r = x * x + y * y;
      values[i * n + j] = (r < 0.25) ? (x * y + x + 3.0) : (x * x * y - 2.0 * x * y * y + 3.0 * x);
    }

  return t8_mra::structured_dataset<2, 1> ({ n, n }, { 0.0, 0.0 }, { h, h }, std::move (values));
}

void
example_dataset ()
{
  root_print ("\n=== MRA: adaptive grid from a structured dataset ===\n");

  constexpr int U = 1;
  constexpr int P = 3;
  const int max_level = 7;
  const double c_thresh = 1.0;
  const std::size_t source_nodes = 129;  // 2^7 + 1: source resolves the level-7 grid

  const auto dataset = build_source_dataset (source_nodes);

  t8_mra::multiscale<T8_ECLASS_QUAD, U, P> mra (max_level, sc_MPI_COMM_WORLD);

  t8_cmesh_t cmesh = t8_cmesh_new_hypercube (T8_ECLASS_QUAD, sc_MPI_COMM_WORLD, 0, 0, 0);
  auto *scheme = t8_scheme_new_default ();
  t8_cmesh_ref (cmesh);
  t8_scheme_ref (const_cast<t8_scheme *> (scheme));

  mra.initialize_data_adaptive (cmesh, scheme, max_level, dataset, t8_mra::hard_thresholding { .c_thresh = c_thresh });
  mra.balance ();

  const auto num_adaptive = t8_forest_get_global_num_leaf_elements (mra.get_forest ());
  const auto num_uniform = static_cast<t8_gloidx_t> (std::pow (4, max_level));
  root_print ("  Adaptive grid: " + std::to_string (num_adaptive) + " elements (uniform level "
              + std::to_string (max_level) + ": " + std::to_string (num_uniform) + ")\n");

  t8_mra::write_forest_lagrange_vtk (mra, "mra_output/dataset_adaptive", P - 1);
  root_print ("  Wrote mra_output/dataset_adaptive (open in ParaView)\n");

  mra.cleanup ();
  t8_cmesh_destroy (&cmesh);
  t8_scheme_unref (const_cast<t8_scheme **> (&scheme));
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

  example_dataset ();

  root_print ("\nDone. Output in mra_output/.\n");

  sc_finalize ();

  mpiret = sc_MPI_Finalize ();
  SC_CHECK_MPI (mpiret);
  return 0;
}

#endif  // T8_ENABLE_MRA
