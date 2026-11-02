/*
  This file is part of t8code.
  t8code is a C library to manage a collection (a forest) of multiple
  connected adaptive space-trees of general element types in parallel.

  Copyright (C) 2026 the developers

  t8code is free software; you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation; either version 2 of the License, or
  (at your option) any later version.

  t8code is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with t8code; if not, write to the Free Software Foundation, Inc.,
  51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.
*/

/** \file t8_time_search_neighbors.cxx
 * Measure the runtime of the algorithms that traverse the leaf arrays of a forest:
 * the top-down search, the creation of the face ghost layer and the computation of leaf face neighbors.
 * The forest is a uniformly refined hypercube that is adaptively refined further, so that the leaves
 * have different levels.
 */

#include <t8.h>
#include <t8_cmesh/t8_cmesh.h>
#include <t8_cmesh/t8_cmesh_examples.h>
#include <t8_forest/t8_forest_general.h>
#include <t8_forest/t8_forest_search/t8_forest_search.hxx>
#include <t8_schemes/t8_default/t8_default.hxx>
#include <sc_flops.h>
#include <sc_statistics.h>
#include <sc_options.h>

/** The timed parts of the benchmark. */
enum t8_time_search_neighbors_stat {
  T8_TIME_SEARCH = 0,     /**< The top-down search over all elements. */
  T8_TIME_GHOST,          /**< Creating a forest with a face ghost layer. */
  T8_TIME_FACE_NEIGHBORS, /**< Computing the leaf face neighbors of all faces of all leaves. */
  T8_TIME_NUM_STATS       /**< The number of timed parts. */
};

/** Refine every element whose linear id is a multiple of three until the maximum level is reached.
 * The linear id does not depend on the partition, so the forest is the same for any number of processes.
 * The maximum level is passed as the user data of \a forest.
 */
static int
t8_time_adapt_every_third (t8_forest_t forest, [[maybe_unused]] t8_forest_t forest_from,
                           [[maybe_unused]] t8_locidx_t which_tree, const t8_eclass_t tree_class,
                           [[maybe_unused]] t8_locidx_t lelement_id, const t8_scheme_c *scheme,
                           [[maybe_unused]] const int is_family, [[maybe_unused]] const int num_elements,
                           t8_element_t *elements[])
{
  const int max_level = *static_cast<const int *> (t8_forest_get_user_data (forest));
  const int level = scheme->element_get_level (tree_class, elements[0]);
  return level < max_level && scheme->element_get_linear_id (tree_class, elements[0], level) % 3 == 0;
}

/** Create the adapted, partitioned forest with a face ghost layer.
 * \param [in] eclass     The element class of the hypercube.
 * \param [in] level      The level of the uniform refinement.
 * \param [in] max_level  The maximum level of the adaptive refinement. Must outlive the forest.
 * \return The committed forest.
 */
static t8_forest_t
t8_time_create_forest (const t8_eclass_t eclass, const int level, int *max_level)
{
  t8_cmesh_t cmesh;
  t8_cmesh_init (&cmesh);
  t8_cmesh_new_hypercube (&cmesh, eclass, sc_MPI_COMM_WORLD, 0, 0, 0);
  const t8_forest_t forest_uniform
    = t8_forest_new_uniform (cmesh, t8_scheme_new_default (), level, 0, sc_MPI_COMM_WORLD);

  t8_forest_t forest;
  t8_forest_init (&forest);
  t8_forest_set_user_data (forest, max_level);
  t8_forest_set_adapt (forest, forest_uniform, t8_time_adapt_every_third, 1);
  t8_forest_set_partition (forest, NULL, 0);
  t8_forest_set_ghost (forest, 1, T8_GHOST_FACES);
  t8_forest_commit (forest);
  return forest;
}

/** Search the whole forest top-down, visiting every element of the refinement hierarchy.
 * \param [in] forest The forest.
 */
static void
t8_time_search (const t8_forest_t forest)
{
  /* The generic lambda converts to the element callback, independent of the type of the leaf array. */
  t8_search<void> search ([] (auto &&...) { return true; }, forest);
  search.do_search ();
}

/** Create a copy of a forest with a face ghost layer.
 * \param [in] forest The forest.
 */
static void
t8_time_ghost (const t8_forest_t forest)
{
  t8_forest_ref (forest);
  t8_forest_t forest_ghost;
  t8_forest_init (&forest_ghost);
  t8_forest_set_copy (forest_ghost, forest);
  t8_forest_set_ghost (forest_ghost, 1, T8_GHOST_FACES);
  t8_forest_commit (forest_ghost);
  t8_forest_unref (&forest_ghost);
}

/** Compute the face neighbors of all faces of all leaves of a forest.
 * \param [in] forest The forest. Must have a face ghost layer.
 * \return The total number of face neighbors found.
 */
static t8_locidx_t
t8_time_face_neighbors (const t8_forest_t forest)
{
  const t8_scheme_c *scheme = t8_forest_get_scheme (forest);
  t8_locidx_t total_num_neighbors = 0;
  const t8_locidx_t num_local_trees = t8_forest_get_num_local_trees (forest);
  for (t8_locidx_t itree = 0; itree < num_local_trees; ++itree) {
    const t8_eclass_t tree_class = t8_forest_get_tree_class (forest, itree);
    const t8_locidx_t num_leaves = t8_forest_get_tree_num_leaf_elements (forest, itree);
    for (t8_locidx_t ileaf = 0; ileaf < num_leaves; ++ileaf) {
      const t8_element_t *leaf = t8_forest_get_leaf_element_in_tree (forest, itree, ileaf);
      const int num_faces = scheme->element_get_num_faces (tree_class, leaf);
      for (int iface = 0; iface < num_faces; ++iface) {
        const t8_element_t **neighbors;
        int *dual_faces;
        int num_neighbors;
        t8_locidx_t *neighbor_indices;
        t8_eclass_t neighbor_class;
        t8_forest_leaf_face_neighbors (forest, itree, leaf, &neighbors, iface, &dual_faces, &num_neighbors,
                                       &neighbor_indices, &neighbor_class);
        total_num_neighbors += num_neighbors;
        if (num_neighbors > 0) {
          T8_FREE (neighbors);
          T8_FREE (dual_faces);
          T8_FREE (neighbor_indices);
        }
      }
    }
  }
  return total_num_neighbors;
}

/** Run the benchmark and print the statistics.
 * \param [in] eclass       The element class of the hypercube.
 * \param [in] level        The level of the uniform refinement.
 * \param [in] adapt_levels The number of additional levels of the adaptive refinement.
 * \param [in] num_runs     The number of times each part is timed.
 */
static void
t8_time_search_neighbors (const t8_eclass_t eclass, const int level, const int adapt_levels, const int num_runs)
{
  int max_level = level + adapt_levels;
  t8_forest_t forest = t8_time_create_forest (eclass, level, &max_level);
  t8_global_essentialf ("Forest with %lli global leaves.\n",
                        static_cast<long long> (t8_forest_get_global_num_leaf_elements (forest)));

  sc_statinfo_t stats[T8_TIME_NUM_STATS];
  sc_stats_init (&stats[T8_TIME_SEARCH], "Search");
  sc_stats_init (&stats[T8_TIME_GHOST], "Ghost");
  sc_stats_init (&stats[T8_TIME_FACE_NEIGHBORS], "Face neighbors");

  long long num_neighbors = 0;
  for (int irun = 0; irun < num_runs; ++irun) {
    sc_flopinfo_t flops;
    sc_flopinfo_t snapshot;
    sc_flops_start (&flops);

    sc_flops_snap (&flops, &snapshot);
    t8_time_search (forest);
    sc_flops_shot (&flops, &snapshot);
    sc_stats_accumulate (&stats[T8_TIME_SEARCH], snapshot.iwtime);

    /* The face ghost layer is computed with a top-down search over the forest,
     * so creating it is timed as well. */
    sc_flops_snap (&flops, &snapshot);
    t8_time_ghost (forest);
    sc_flops_shot (&flops, &snapshot);
    sc_stats_accumulate (&stats[T8_TIME_GHOST], snapshot.iwtime);

    sc_flops_snap (&flops, &snapshot);
    num_neighbors = t8_time_face_neighbors (forest);
    sc_flops_shot (&flops, &snapshot);
    sc_stats_accumulate (&stats[T8_TIME_FACE_NEIGHBORS], snapshot.iwtime);
  }
  long long global_num_neighbors;
  const int mpiret
    = sc_MPI_Allreduce (&num_neighbors, &global_num_neighbors, 1, sc_MPI_LONG_LONG_INT, sc_MPI_SUM, sc_MPI_COMM_WORLD);
  SC_CHECK_MPI (mpiret);
  t8_global_essentialf ("Found %lli face neighbors.\n", global_num_neighbors);

  sc_stats_compute (sc_MPI_COMM_WORLD, T8_TIME_NUM_STATS, stats);
  sc_stats_print (t8_get_package_id (), SC_LP_ESSENTIAL, T8_TIME_NUM_STATS, stats, 1, 1);
  t8_forest_unref (&forest);
}

int
main (int argc, char **argv)
{
  int mpiret = sc_MPI_Init (&argc, &argv);
  SC_CHECK_MPI (mpiret);
  sc_init (sc_MPI_COMM_WORLD, 1, 1, NULL, SC_LP_ESSENTIAL);
  t8_init (SC_LP_ESSENTIAL);

  int help;
  int eclass_int;
  int level;
  int adapt_levels;
  int num_runs;
  sc_options_t *opt = sc_options_new (argv[0]);
  sc_options_add_switch (opt, 'h', "help", &help, "Display a short help message.");
  sc_options_add_int (opt, 'e', "elements", &eclass_int, T8_ECLASS_HEX,
                      "The element class of the hypercube.\n"
                      "\t\t\t\t\t2 - quad\n\t\t\t\t\t3 - triangle\n\t\t\t\t\t4 - hexahedron\n"
                      "\t\t\t\t\t5 - tetrahedron\n\t\t\t\t\t6 - prism\n\t\t\t\t\t7 - pyramid");
  sc_options_add_int (opt, 'l', "level", &level, 3, "The level of the uniform refinement.");
  sc_options_add_int (opt, 'a', "adapt-levels", &adapt_levels, 2,
                      "The number of additional levels of the adaptive refinement.");
  sc_options_add_int (opt, 'r', "runs", &num_runs, 5, "The number of times each part is timed.");

  const int parsed = sc_options_parse (t8_get_package_id (), SC_LP_ERROR, opt, argc, argv);
  if (help) {
    sc_options_print_usage (t8_get_package_id (), SC_LP_ERROR, opt, NULL);
  }
  else if (parsed >= 0 && T8_ECLASS_QUAD <= eclass_int && eclass_int <= T8_ECLASS_PYRAMID && level >= 0
           && adapt_levels >= 0 && num_runs > 0) {
    t8_time_search_neighbors (static_cast<t8_eclass_t> (eclass_int), level, adapt_levels, num_runs);
  }
  else {
    t8_global_productionf ("\n\tERROR: Wrong usage.\n\n");
    sc_options_print_usage (t8_get_package_id (), SC_LP_ERROR, opt, NULL);
  }

  sc_options_destroy (opt);
  sc_finalize ();
  mpiret = sc_MPI_Finalize ();
  SC_CHECK_MPI (mpiret);
  return 0;
}
