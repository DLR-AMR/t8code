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
 * \file In this file, we test the weighted-partitioning feature of t8code.
 *
 * It allows defining a weighting function that assigns a weight to each element
 * in the forest. The partitioning will then try to establish a (as far as possible)
 * ideal balancing of the weights, rather than the element numbers as "usual".
 *
 * The following steps are performed per example cmesh and scheme, with sanity
 * checks validating each of them:
 *
 * (1.) A uniform, partitioned base forest is created.
 *
 * (2.) In a first test, the weighting function is equal to one for all elements.
 *      This allows a simple sanity check: The resulting partition has to be the
 *      same as without assigning weights.
 *
 * (3.) Next, we test a non-uniform weight function. As an example, we defined the
 *      weight function such that an element's weight is simply its global element id.
 *      This allows to "manually" compute the expected element offsets to then compare 
 *      them.
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
 * Test class for validating t8code's weighted partitioning. See file description for details on the testing strategy.
*/
class t8_test_weighted_partitioning_test: public testing::TestWithParam<std::tuple<int, cmesh_example_base *>> {

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
  const t8_scheme *scheme;
  t8_cmesh_t cmesh;
};

/**
 * Main function of testing suite. See doxygen file description for details on the testing strategy.
*/
TEST_P (t8_test_weighted_partitioning_test, test_weighted_partitioning)
{

  // -------------------------------------------
  // ----- (1.) Create uniform base forest -----
  // -------------------------------------------

  t8_debugf ("Create uniform base forest.\n");

  // Initial uniform refinement level
  const int level = 2;

  // Increase reference counters of cmesh and scheme to avoid reaching zero.
  t8_cmesh_ref (cmesh);
  scheme->ref ();

  // Create initial, uniform base forest.
  t8_forest_t base_forest = t8_forest_new_uniform (cmesh, scheme, level, 0, sc_MPI_COMM_WORLD);

  // ---------------------------------------------
  // ----- (2.) Test 1: All weights the same -----
  // ---------------------------------------------
  t8_debugf ("Create forest with uniform partition weights.\n");

  // Initialize forest
  t8_forest_t forest_uniform_weights;
  t8_forest_init (&forest_uniform_weights);

  // Define uniform weight fcn.
  t8_weight_fcn_t *uniform_weight_fcn = [] (t8_forest_t, t8_locidx_t, t8_locidx_t) -> double { return 1; };

  // Prepare partitioning.
  t8_forest_set_partition (forest_uniform_weights, base_forest, 0);
  t8_forest_set_partition_weight_function (forest_uniform_weights, uniform_weight_fcn);

  // Commit forest (without destroying base_forest).
  t8_forest_ref (base_forest);
  t8_forest_commit (forest_uniform_weights);

  // Create partitioned forest without weights for comparison.
  t8_forest_t forest_no_weights;
  t8_forest_init (&forest_no_weights);
  t8_forest_set_partition (forest_no_weights, base_forest, 0);
  t8_forest_ref (base_forest);
  t8_forest_commit (forest_no_weights);

  // Make sure they are the same.
  EXPECT_FOREST_EQ (forest_uniform_weights, forest_no_weights);

  // ---------------------------------------------
  // ----- (3.) Test 2: Non-uniform weights ------
  // ---------------------------------------------
  t8_debugf ("Create forest with non-uniform partition weights.\n");

  // (3a.) Apply weighted partitioning:
  // ----------------------------------

  // Define non-uniform weight fcn: For this test, we pick the element number as element weight.
  t8_weight_fcn_t *weight_fcn = [] (t8_forest_t forest, t8_locidx_t ltree_id, t8_locidx_t ele_in_tree) -> double {
    const t8_gloidx_t gelem_id = t8_forest_get_first_local_leaf_element_id (forest)
                                 + t8_forest_get_tree_element_offset (forest, ltree_id) + ele_in_tree;
    return gelem_id;
  };

  // Prepare forest.
  t8_forest_t forest_weighted;
  t8_forest_init (&forest_weighted);
  t8_forest_set_partition (forest_weighted, base_forest, 0);
  t8_forest_set_partition_weight_function (forest_weighted, weight_fcn);

  // Commit forest (without destroying base_forest).
  t8_forest_ref (base_forest);
  t8_forest_commit (forest_weighted);

  // (3b.) Determine expected partitioning:
  // --------------------------------------

  // Preparation: MPI variables and global number of leaf elements.
  const int mpirank = base_forest->mpirank;
  const int mpisize = base_forest->mpisize;
  const t8_gloidx_t num_global_leaves = t8_forest_get_global_num_leaf_elements (base_forest);

  // Total sum of the weights for this case (Gauss sum formula for num_global_leaves - 1)
  const double total_weight_sum = num_global_leaves * (num_global_leaves - 1) / 2;

  // "Manually" compute the element offsets expected for this distribution.
  // Note: For simplicity, every process computes all offsets here, so no MPI communication is required.
  std::vector<t8_gloidx_t> element_offsets (mpisize + 1);
  double sum = 0.0;
  int on_rank = 0;
  element_offsets[0] = 0;
  for (t8_gloidx_t ielem = 0; ielem < num_global_leaves; ielem++) {

    // Add weight (which matches element number here)
    sum += ielem;

    // If we are above ideal bound for current process, add offset to vector and go to next rank.
    if (sum > total_weight_sum * (on_rank + 1) / mpisize) {
      on_rank++;
      std::fill (element_offsets.begin () + on_rank, element_offsets.end (), ielem);
      t8_global_productionf ("Computed offset for rank %i: %li \n", on_rank, element_offsets[on_rank]);
    }
  }
  // The last entry is always the global number of elements.
  element_offsets[mpisize] = num_global_leaves;

  // Determine number of local elements.
  t8_locidx_t num_local_leaves = element_offsets[mpirank + 1] - element_offsets[mpirank];

  // (3c.) Comparison:
  // -----------------

  // Make sure the number of local leaves match.
  EXPECT_EQ (t8_forest_get_local_num_leaf_elements (forest_weighted), num_local_leaves);

  // ---------------------------
  // ----- Memory clean-up -----
  // ---------------------------

  // Destroy the forests.
  t8_forest_unref (&base_forest);
  t8_forest_unref (&forest_uniform_weights);
  t8_forest_unref (&forest_no_weights);
  t8_forest_unref (&forest_weighted);
}

// Instantiate parameterized test to be run for all schemes and example cmeshes.
INSTANTIATE_TEST_SUITE_P (t8_gtest_weighted_partitioning, t8_test_weighted_partitioning_test,
                          testing::Combine (AllSchemeCollections, AllCmeshsParam), pretty_print_base_example_scheme);
