/*
  This file is part of t8code.
  t8code is a C library to manage a collection (a forest) of multiple
  connected adaptive space-trees of general element classes in parallel.

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

/**
 * \file This file adds a test for the forest profiling. For a simple example forest,
 *       it computes all profiling stats, accesses them and runs some sanity checks.
 */

// testing
#include <gtest/gtest.h>
#include <test/t8_gtest_custom_assertion.hxx>
#include <test/t8_gtest_adapt_callbacks.hxx>

// t8code
#include <t8.h>
#include <t8_cmesh/t8_cmesh.h>
#include <t8_cmesh/t8_cmesh_examples.h>
#include <t8_forest/t8_forest_types.h>
#include <t8_forest/t8_forest_profiling.h>

/**
 * Main function of testing suite. See doxygen file description for details on the testing strategy.
*/
TEST (t8_test_forest_profiling_test, test_forest_profiling)
{

  // (0) Create cmesh and initial forest.
  // ------------------------------------
  t8_cmesh_t cmesh;
  t8_cmesh_init (&cmesh);
  t8_cmesh_new_hypercube_hybrid (cmesh, sc_MPI_COMM_WORLD, 0);
  t8_forest_t uniform_forest = t8_forest_new_uniform (cmesh, t8_scheme_new_default (), 1, true, sc_MPI_COMM_WORLD);

  // (1) Define a second forest and activate profiling.
  // --------------------------------------------------
  t8_forest_t profiling_forest;
  t8_forest_init (&profiling_forest);
  t8_forest_set_profiling (profiling_forest, true);

  // (2) Operate on that forest.
  // ---------------------------
  // (2a) adapt, partition, balance, ghost, ghost-exchange data
  t8_forest_set_adapt (profiling_forest, uniform_forest, refine_every_nth_element_callback<5, 1>, false);
  t8_forest_set_partition (profiling_forest, uniform_forest, false);
  t8_forest_set_balance (profiling_forest, uniform_forest, true);
  t8_forest_set_ghost (profiling_forest, true, T8_GHOST_FACES);
  t8_forest_commit (profiling_forest);

  // (2b) Test ghost exchange data.
  sc_array_t element_data;
  t8_locidx_t num_data_entries
    = t8_forest_get_local_num_leaf_elements (profiling_forest) + t8_forest_get_num_ghosts (profiling_forest);
  sc_array_init_size (&element_data, sizeof (int), num_data_entries);
  for (t8_locidx_t ielem = 0; ielem < num_data_entries; ielem++) {
    *(int*) t8_sc_array_index_locidx (&element_data, ielem) = 1;
  }
  t8_forest_ghost_exchange_data (profiling_forest, &element_data);

  // (3) Compute and print profiling stats.
  // --------------------------------------
  t8_forest_compute_profile (profiling_forest);
  t8_forest_print_profile (profiling_forest);

  // (4) Access the specific stats computed.
  // ---------------------------------------
  // Stats
  const sc_statinfo_t* adapt_stats = t8_forest_profile_get_adapt_stats (profiling_forest);
  const sc_statinfo_t* partition_stats = t8_forest_profile_get_partition_stats (profiling_forest);
  const sc_statinfo_t* commit_stats = t8_forest_profile_get_commit_stats (profiling_forest);
  const sc_statinfo_t* balance_stats = t8_forest_profile_get_balance_stats (profiling_forest);
  const sc_statinfo_t* balance_rounds_stats = t8_forest_profile_get_balance_rounds_stats (profiling_forest);
  [[maybe unused]] const sc_statinfo_t* ghost_stats = t8_forest_profile_get_ghost_stats (profiling_forest);

  // Integers
  int balance_rounds, ghosts_sent, procs_sent;

  // Runtimes
  double time_adapt = t8_forest_profile_get_adapt_time (profiling_forest);
  double time_partition = t8_forest_profile_get_partition_time (profiling_forest, &procs_sent);
  double time_balance = t8_forest_profile_get_balance_time (profiling_forest, &balance_rounds);
  double time_cmesh_offsets = t8_forest_profile_get_cmesh_offsets_runtime (profiling_forest);
  double time_forest_offsets = t8_forest_profile_get_forest_offsets_runtime (profiling_forest);
  double time_first_descendant = t8_forest_profile_get_first_descendant_runtime (profiling_forest);
  double time_ghost = t8_forest_profile_get_ghost_time (profiling_forest, &ghosts_sent);
  double time_ghost_wait = t8_forest_profile_get_ghostexchange_waittime (profiling_forest);

  // (5) Simple sanity checks.
  // -------------------------
  // Preparation: get mpisize
  int mpisize;
  int mpiret = sc_MPI_Comm_size (sc_MPI_COMM_WORLD, &mpisize);
  SC_CHECK_MPI (mpiret);

  // Runtimes
  EXPECT_GT (time_adapt, 0.0);
  EXPECT_GT (time_partition, 0.0);
  EXPECT_GT (time_balance, 0.0);
  EXPECT_GT (time_cmesh_offsets, 0.0);
  EXPECT_GT (time_forest_offsets, 0.0);
  EXPECT_GT (time_first_descendant, 0.0);

  // Integers
  EXPECT_GT (balance_rounds, 0);
  EXPECT_GE (procs_sent, 0);
  EXPECT_LT (procs_sent, mpisize);

  // Tests only valid for parallel execution
  if (mpisize > 1) {
    EXPECT_GT (ghosts_sent, 0);
    EXPECT_GT (time_ghost, 0.0);
    EXPECT_GT (time_ghost_wait, 0.0);
  }

  // Stats
  EXPECT_GT (commit_stats->average, 0.0);

  // (6) Sanity checks comparing stats and times
  // -------------------------------------------
  EXPECT_GE (time_adapt, adapt_stats->min);
  EXPECT_LE (time_adapt, adapt_stats->max);
  EXPECT_GE (time_partition, partition_stats->min);
  EXPECT_LE (time_partition, partition_stats->max);
  EXPECT_GE (time_balance, balance_stats->min);
  EXPECT_LE (time_balance, balance_stats->max);
  EXPECT_GE (balance_rounds, balance_rounds_stats->min);
  EXPECT_LE (balance_rounds, balance_rounds_stats->max);
#if T8_ENABLE_MPI
  EXPECT_GE (time_ghost, ghost_stats->min);
  EXPECT_LE (time_ghost, ghost_stats->max);
#endif

  // (7) Memory clean up.
  // --------------------
  t8_forest_unref (&profiling_forest);
  sc_array_reset (&element_data);
}
