/*
  This file is part of t8code.
  t8code is a C library to manage a collection (a forest) of multiple
  connected adaptive space-trees of general element classes in parallel.

  Copyright (C) 2025 the developers

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
 * \file In this file, we test the functions t8_forest_set_partition_offset and
 * t8_forest_new_gather.
 *
 * They are closely related and allow to partition a forest according to custom
 * element offsets (t8_forest_set_partition_offset) or such that all elements
 * are gathered on one process (t8_forest_new_gather).
 *
 * The following steps are performed per example cmesh and scheme, with sanity
 * checks validating each of them:
 *
 * (1.) A uniform, partitioned base forest is created.
 *
 * (2.) A gathered forest is created, i.e., a copy of the base forest in which
 *      all elements and trees live on one rank.
 *
 * (3.) A forest with custom element offsets is created. We choose some example
 *      distribution such that each process gets a different number of elements.
 *
 * (4.) Lastly, as a sanity check, the forest with custom offsets is repartitioned
 *      back to the standard partitioning to verify the result matches the
 *      original base forest.
 *
 *  Note that "under the hood", t8_forest_new_gather is based on
 *  t8_forest_set_partition_offset, which is why we do not have to perform every
 *  check twice.
 */

// testing
#include <gtest/gtest.h>
#include <test/t8_gtest_schemes.hxx>
#include "test/t8_cmesh_generator/t8_cmesh_example_sets.hxx"
#include <test/t8_gtest_macros.hxx>
#include <test/t8_gtest_custom_assertion.hxx>

// t8code
#include <t8.h>
#include <t8_cmesh/t8_cmesh.h>
#include <t8_forest/t8_forest_general.h>
#include <t8_forest/t8_forest_partition.h>
#include <t8_forest/t8_forest_types.h>

/**
 * Test class for validating \ref t8_forest_set_partition_offset and \ref t8_forest_new_gather.
*/
class t8_test_set_partition_offset_test: public testing::TestWithParam<std::tuple<int, cmesh_example_base *>> {

 protected:
  /**
   * Setup routine of the test suite.
  */
  void
  SetUp () override
  {
    // Get scheme.
    const int scheme_id = std::get<0> (GetParam ());
    scheme = create_from_scheme_id (scheme_id);

    // Construct cmesh.
    cmesh = std::get<1> (GetParam ())->cmesh_create ();

    // Skip empty meshes.
    if (t8_cmesh_is_empty (cmesh)) {
      GTEST_SKIP ();
    }
  }

  /**
   * Tear down routine of the test suite.
   *
   * Makes sure to unref the cmesh and the scheme.
  */
  void
  TearDown () override
  {
    t8_cmesh_destroy (&cmesh);
    scheme->unref ();
  }

  // Member variables:
  const t8_scheme *scheme; /**< The currently tested scheme. */
  t8_cmesh_t cmesh;        /**< The currently tested cmesh. */
};

/**
 * Main function of testing suite. See doxygen file description for details on the testing strategy.
*/
TEST_P (t8_test_set_partition_offset_test, test_set_partition_offset)
{

  // -------------------------------------------
  // ----- (1.) Create uniform base forest -----
  // -------------------------------------------

  t8_global_productionf ("Create uniform base forest.\n");

  // Initial uniform refinement level
  const int level = 2;

  // Increase reference counters of cmesh and scheme to avoid reaching zero.
  t8_cmesh_ref (cmesh);
  scheme->ref ();

  // Create initial, uniform base forest.
  t8_forest_t base_forest = t8_forest_new_uniform (cmesh, scheme, level, 0, sc_MPI_COMM_WORLD);

  // ----------------------------------------
  // ----- (2.) Gather forest on rank 0 -----
  // ----------------------------------------

  t8_global_productionf ("Gather forest on rank 0.\n");

  // Gather forest on rank 0.
  t8_forest_ref (base_forest);
  t8_forest_t forest_gathered = t8_forest_new_gather (base_forest, 0);

  // Sanity checks: forest_gathered vs. base_forest
  EXPECT_EQ (forest_gathered->global_num_leaf_elements, base_forest->global_num_leaf_elements);
  EXPECT_EQ (forest_gathered->global_num_trees, base_forest->global_num_trees);

  // Get MPI rank
  int mpirank;
  sc_MPI_Comm_rank (sc_MPI_COMM_WORLD, &mpirank);

  // Check that all elements and trees have been gathered on rank 0.
  if (mpirank == 0) {
    EXPECT_EQ (forest_gathered->local_num_leaf_elements, forest_gathered->global_num_leaf_elements);
    EXPECT_EQ (t8_forest_get_num_local_trees (forest_gathered), forest_gathered->global_num_trees);
  }
  else {
    EXPECT_EQ (forest_gathered->local_num_leaf_elements, 0);
    EXPECT_EQ (t8_forest_get_num_local_trees (forest_gathered), 0);
  }

  // ----------------------------------------------
  // ----- (3.) Gather forest on another rank -----
  // ----------------------------------------------

  // Get number of processes
  int mpisize;
  sc_MPI_Comm_size (sc_MPI_COMM_WORLD, &mpisize);

  // Choose some nonzero rank to gather on.
  int other_gather_rank = mpisize / 2;

  t8_global_productionf ("Gather forest on rank .\n");

  // Gather forest on rank 0.
  t8_forest_ref (base_forest);
  t8_forest_t forest_gathered_other = t8_forest_new_gather (base_forest, other_gather_rank);

  // Sanity checks: forest_gathered vs. base_forest
  EXPECT_EQ (forest_gathered_other->global_num_leaf_elements, base_forest->global_num_leaf_elements);
  EXPECT_EQ (forest_gathered_other->global_num_trees, base_forest->global_num_trees);

  // Check that all elements and trees have been gathered on rank other_gather_rank.
  if (mpirank == other_gather_rank) {
    EXPECT_EQ (forest_gathered_other->local_num_leaf_elements, forest_gathered_other->global_num_leaf_elements);
    EXPECT_EQ (t8_forest_get_num_local_trees (forest_gathered_other), forest_gathered_other->global_num_trees);
  }
  else {
    EXPECT_EQ (forest_gathered_other->local_num_leaf_elements, 0);
    EXPECT_EQ (t8_forest_get_num_local_trees (forest_gathered_other), 0);
  }

  // ----------------------------------------------
  // ----- (4.) Test some manual partitioning -----
  // ----------------------------------------------

  t8_debugf ("Create forest with manual partitioning.\n");

  // Initialization
  t8_forest_t forest_manual_partition;
  t8_forest_init (&forest_manual_partition);

  // Set custom element offsets for partitioning
  // We set the first element to N * p*p / (P*P),
  // where N is the global number of elements, p the current process,
  // and P the number of processes.
  t8_gloidx_t my_element_offset = base_forest->global_num_leaf_elements * mpirank * mpirank / (mpisize * mpisize);
  t8_gloidx_t next_element_offset
    = base_forest->global_num_leaf_elements * (mpirank + 1) * (mpirank + 1) / (mpisize * mpisize);

  // Pass offset to forest.
  t8_forest_set_partition_offset (forest_manual_partition, my_element_offset);

  // Create forest with custom partitioning.
  t8_forest_set_partition (forest_manual_partition, base_forest, 0);
  t8_forest_ref (base_forest);
  t8_forest_commit (forest_manual_partition);

  // Sanity checks: forest_manual_partition vs. base_forest
  EXPECT_EQ (forest_manual_partition->global_num_leaf_elements, base_forest->global_num_leaf_elements);
  EXPECT_EQ (forest_manual_partition->global_num_trees, base_forest->global_num_trees);

  // Check element numbers
  EXPECT_EQ (forest_manual_partition->local_num_leaf_elements, next_element_offset - my_element_offset);

  // -----------------------------------------------
  // ----- (5.) Sanity check: Repartition back -----
  // -----------------------------------------------

  // Repartition the manually partitioned forest
  t8_forest_t forest_repartitioned;
  t8_forest_init (&forest_repartitioned);
  t8_forest_set_partition (forest_repartitioned, forest_manual_partition, 0);
  t8_forest_ref (forest_manual_partition);
  t8_forest_commit (forest_repartitioned);

  // It has to be the same as the base forest
  EXPECT_FOREST_EQ (base_forest, forest_repartitioned);

  // Destroy the forests.
  t8_forest_unref (&base_forest);
  t8_forest_unref (&forest_manual_partition);
  t8_forest_unref (&forest_repartitioned);
  t8_forest_unref (&forest_gathered);
  t8_forest_unref (&forest_gathered_other);
}

// Instantiate parameterized test to be run for all schemes and example cmeshes.
INSTANTIATE_TEST_SUITE_P (t8_gtest_set_partition_offset, t8_test_set_partition_offset_test,
                          testing::Combine (AllSchemeCollections, AllCmeshsParam), pretty_print_base_example_scheme);
