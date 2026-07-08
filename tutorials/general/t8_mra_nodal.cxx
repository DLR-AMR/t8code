/**
 * @file t8_mra_nodal.cxx
 * @brief MRA tutorial: adapt a nodal DG forest and get a nodal forest back
 *
 * A nodal DG code stores its solution as values at nodal points per cell; the
 * MRA stores modal coefficients. This tutorial takes a forest already carrying
 * the caller's nodal data (here two components), runs the MRA grid adaptation
 * on it, and returns a new forest on the adapted grid carrying the
 * reconstructed nodal data in the caller's own per-cell layout (attached as
 * forest user data). The returned forest is self-contained and outlives the
 * multiscale object.
 *
 * The bridge is two converters: initialize_data_nodal (nodal -> modal on load)
 * and export_data_nodal (modal -> nodal on the adapted grid). The nodal
 * solution is written to VTK before and after adaptation.
 */

#include "t8.h"
#ifdef T8_ENABLE_MRA

#include "t8_mra/t8_mra.hxx"
#include "t8_cmesh/t8_cmesh.h"
#include "t8_cmesh/t8_cmesh_examples.h"
#include "t8_forest/t8_forest.h"
#include "t8_forest/t8_forest_ghost.h"
#include "t8_forest/t8_forest_geometrical.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <filesystem>
#include <iostream>
#include <span>
#include <vector>

constexpr int U = 2;
constexpr int P = 3;
constexpr int DOF = P * P;

/// The caller's nodal DG cell state: values at the DOF nodes per component,
/// component-major (index u*DOF + j), matching the MRA nodal buffer layout.
struct nodal_cell
{
  std::array<double, U * DOF> values;
};

using node_set = std::array<std::array<double, 2>, DOF>;

static void
root_print (const std::string &s)
{
  int rank;
  sc_MPI_Comm_rank (sc_MPI_COMM_WORLD, &rank);
  if (rank == 0)
    std::cout << s;
}

/// Two-component analytic field the fake nodal solver samples. Both components
/// are smooth over large regions with localized sharp features, so the MRA
/// coarsens the smooth areas and keeps the grid fine only around the features:
/// component 0 a Gaussian bump, component 1 a jump ring on a gentle background.
static std::array<double, U>
source_field (double x, double y)
{
  const auto rc = std::hypot (x - 0.35, y - 0.35);
  const auto bump = std::exp (-120.0 * rc * rc);
  const auto ring = (std::hypot (x - 0.7, y - 0.7) < 0.2) ? 1.0 : 0.0;
  const auto background = 0.2 * std::sin (M_PI * x) * std::sin (M_PI * y);
  return { bump + background, ring + background };
}

/// Fake a nodal solution on `forest`: one nodal_cell per leaf (SFC order),
/// sampled from source_field at the node coordinates.
static std::vector<nodal_cell>
sample_nodal (t8_forest_t forest, const node_set &nodes)
{
  std::vector<nodal_cell> data;
  const auto num_trees = t8_forest_get_num_local_trees (forest);
  for (t8_locidx_t tree = 0; tree < num_trees; ++tree) {
    const auto num_leaves = t8_forest_get_tree_num_leaf_elements (forest, tree);
    for (t8_locidx_t e = 0; e < num_leaves; ++e) {
      const auto *element = t8_forest_get_leaf_element_in_tree (forest, tree, e);
      double v0[3], v3[3];  // axis-aligned box: vertex 0 min, vertex 3 max
      t8_forest_element_coordinate (forest, tree, element, 0, v0);
      t8_forest_element_coordinate (forest, tree, element, 3, v3);

      nodal_cell cell;
      for (int j = 0; j < DOF; ++j) {
        const auto x = v0[0] + nodes[j][0] * (v3[0] - v0[0]);
        const auto y = v0[1] + nodes[j][1] * (v3[1] - v0[1]);
        const auto f = source_field (x, y);

        for (int u = 0; u < U; ++u)
          cell.values[u * DOF + j] = f[u];
      }

      data.push_back (cell);
    }
  }

  return data;
}

/// Refine cells the solution actually needs: a tight ball on the bump and a
/// band along the jump ring. Leaves the original nodal forest nonuniform but
/// resolving the features (the loader accepts any committed forest).
static int
refine_near_features (t8_forest_t, t8_forest_t forest_from, t8_locidx_t which_tree, t8_eclass_t tree_class, t8_locidx_t,
                      const t8_scheme_c *scheme, int, int, t8_element_t *elements[])
{
  constexpr int feature_level = 7;
  if (scheme->element_get_level (tree_class, elements[0]) >= feature_level)
    return 0;

  double c[3];
  t8_forest_element_centroid (forest_from, which_tree, elements[0], c);
  const bool on_bump = std::hypot (c[0] - 0.35, c[1] - 0.35) < 0.12;
  const bool on_ring = std::abs (std::hypot (c[0] - 0.7, c[1] - 0.7) - 0.2) < 0.06;
  return (on_bump || on_ring) ? 1 : 0;
}

/// A nonuniform original forest: a uniform grid refined toward the features.
/// Consumes one ref of cmesh and scheme (like t8_forest_new_uniform).
static t8_forest_t
feature_refined_forest (t8_cmesh_t cmesh, const t8_scheme *scheme, sc_MPI_Comm comm)
{
  t8_forest_t forest = t8_forest_new_uniform (cmesh, scheme, 4, 0, comm);
  return t8_forest_new_adapt (forest, refine_near_features, 1, 0, nullptr);
}

/**
 * Run the MRA grid adaptation on a nodal DG forest and return the adapted grid
 * carrying the reconstructed nodal data. The returned forest has a face-ghost
 * layer, and the attached std::vector<nodal_cell> is ghost-complete: entries
 * [0, num_local) are local leaves, [num_local, num_local + num_ghost) the ghost
 * cells (filled by one ghost exchange), so it is ready for a nodal solver. The
 * forest owns the vector as user data (caller deletes it and unrefs the
 * forest). `forest_in` / `nodal_in` stay owned by the caller.
 */
static t8_forest_t
mra_coarsen_nodal (t8_forest_t forest_in, const std::vector<nodal_cell> &nodal_in, const node_set &nodes, int min_level,
                   int max_level)
{
  t8_mra::multiscale<T8_ECLASS_QUAD, U, P> mra (max_level, sc_MPI_COMM_WORLD);

  // Load: hand each leaf's nodal values to the MRA (SFC order matches nodal_in).
  size_t in_idx = 0;
  mra.initialize_data_nodal (forest_in, nodes, [&] (int, const t8_element_t *) { return nodal_in[in_idx++].values; });
  t8_mra::write_forest_lagrange_vtk (mra, "mra_output/nodal_before", P - 1);

  mra.coarsen (min_level, max_level, t8_mra::hard_thresholding { .c_thresh = 0.1 });
  mra.balance ();
  t8_mra::write_forest_lagrange_vtk (mra, "mra_output/nodal_after", P - 1);

  // Reconstruct nodal values on the adapted grid (local leaves, SFC order).
  t8_forest_t adapted = mra.get_forest ();
  const auto num_adapted = t8_forest_get_local_num_leaf_elements (adapted);
  auto *nodal_adapted = sc_array_new_count (sizeof (nodal_cell), num_adapted);
  auto *adapted_cells = reinterpret_cast<nodal_cell *> (nodal_adapted->array);
  size_t leaf = 0;

  mra.export_data_nodal (nodes, [&] (int, const t8_element_t *, std::span<const double> nodal) {
    std::copy (nodal.begin (), nodal.end (), adapted_cells[leaf++].values.begin ());
  });

  // Returned forest: the adapted grid load-balanced with a face-ghost layer.
  // set_partition takes ownership of `adapted`, so ref it to keep the MRA's.
  t8_forest_ref (adapted);
  t8_forest_t forest_out;
  t8_forest_init (&forest_out);
  t8_forest_set_ghost (forest_out, 1, T8_GHOST_FACES);
  t8_forest_set_partition (forest_out, adapted, 0);
  t8_forest_commit (forest_out);

  const auto num_local = t8_forest_get_local_num_leaf_elements (forest_out);
  const auto num_ghost = t8_forest_get_num_ghosts (forest_out);

  // Move the nodal data onto the new partition, then fill ghosts by exchange.
  auto *nodal_partitioned = sc_array_new_count (sizeof (nodal_cell), num_local);
  t8_forest_partition_data (adapted, forest_out, nodal_adapted, nodal_partitioned);

  auto *nodal_with_ghosts = sc_array_new_count (sizeof (nodal_cell), num_local + num_ghost);
  std::copy_n (reinterpret_cast<nodal_cell *> (nodal_partitioned->array), num_local,
               reinterpret_cast<nodal_cell *> (nodal_with_ghosts->array));
  t8_forest_ghost_exchange_data (forest_out, nodal_with_ghosts);

  auto *nodal_out = new std::vector<nodal_cell> (num_local + num_ghost);
  std::copy_n (reinterpret_cast<nodal_cell *> (nodal_with_ghosts->array), num_local + num_ghost, nodal_out->data ());

  sc_array_destroy (nodal_adapted);
  sc_array_destroy (nodal_partitioned);
  sc_array_destroy (nodal_with_ghosts);
  mra.cleanup ();

  t8_forest_set_user_data (forest_out, nodal_out);

  return forest_out;
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

  root_print ("\n=== MRA: adapt a two-component nodal DG forest ===\n");

  const int max_level = 7;
  const int min_level = 0;

  // Reference nodes: tensor of P equispaced points on [0,1]^2 (node j = jy*P+jx).
  node_set nodes;
  for (int jy = 0; jy < P; ++jy)
    for (int jx = 0; jx < P; ++jx)
      nodes[jy * P + jx] = { jx / static_cast<double> (P - 1), jy / static_cast<double> (P - 1) };

  // The caller's nodal DG state: a nonuniform forest plus one nodal_cell per leaf.
  t8_cmesh_t cmesh = t8_cmesh_new_hypercube (T8_ECLASS_QUAD, sc_MPI_COMM_WORLD, 0, 0, 0);
  auto *scheme = t8_scheme_new_default ();
  t8_forest_t forest_in = feature_refined_forest (cmesh, scheme, sc_MPI_COMM_WORLD);
  const auto nodal_in = sample_nodal (forest_in, nodes);
  root_print ("  Input nodal forest:  " + std::to_string (nodal_in.size ()) + " local cells (plotted: nodal_before)\n");

  // Adapt with the MRA; get the coarser grid back with nodal data attached.
  t8_forest_t forest_out = mra_coarsen_nodal (forest_in, nodal_in, nodes, min_level, max_level);
  auto *nodal_out = static_cast<std::vector<nodal_cell> *> (t8_forest_get_user_data (forest_out));
  const auto num_local = t8_forest_get_local_num_leaf_elements (forest_out);
  const auto num_ghost = t8_forest_get_num_ghosts (forest_out);
  root_print ("  Output nodal forest: " + std::to_string (num_local) + " local + " + std::to_string (num_ghost)
              + " ghost cells (plotted: nodal_after)\n");

  t8_forest_unref (&forest_in);

  // forest_out + *nodal_out is now a self-contained nodal DG state ready to hand
  // back to the solver. Release it once done.
  delete nodal_out;
  t8_forest_unref (&forest_out);

  root_print ("\nDone. Output in mra_output/.\n");

  sc_finalize ();
  mpiret = sc_MPI_Finalize ();
  SC_CHECK_MPI (mpiret);
  return 0;
}

#endif  // T8_ENABLE_MRA
