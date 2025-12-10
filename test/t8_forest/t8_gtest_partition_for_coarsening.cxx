/*
  This file is part of t8code.
  t8code is a C library to manage a collection (a forest) of multiple
  connected adaptive space-trees of general element classes in parallel.

  Copyright (C) 2015 the developers

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
 * \file In this file, we test the partition-for-coarsening (PFC) functionality
 * of t8code. Its main motivation is to avoid splitting families at process
 * boundaries. Such a split family is problematic in that it is not recognized
 * by the adaptation callback, making coarsening dependent on the number of
 * processes etc.
 *
 * This goal is achieved by slightly softening the constraint of equal load
 * distribution: If the PFC flag is set in t8_forest_set_partition, the
 * standard partitioning is followed by a correction step that checks for all
 * process boundaries whether a family is split; if so, the boundary is slightly
 * moved to make sure the family is on one process.
 *
 * The test validating the PFC functionality takes the following steps:
 *
 * (1.) Create a uniform forest from the cmesh.
 *
 * (2.) Adapt the forest based on some arbitrary criterion.
 *
 * (3.) Apply the PFC partitioning
 *
 * (4.) Coarsen all families that are on one process (should be all families).
 *
 * (5.) Perform checks:
 *
 *    (5a.) Repartition the forest from step (4.) without PFC option.
 *
 *    (5b.) Gather forest from step (2.) on one rank and then coarsen all families
 *      (since this whole forest lives on one rank, we know no family is split);
 *      after that, this coarsened forest is repartitioned without PFC option.
 *
 *    (5c.) The forests resulting from (5a.) and (5b.) have to match!
 */

// testing
#include <gtest/gtest.h>
#include <test/t8_gtest_schemes.hxx>
#include <test/t8_gtest_macros.hxx>
#include "test/t8_cmesh_generator/t8_cmesh_example_sets.hxx"

// t8code
#include <t8.h>
#include <t8_forest/t8_forest_general.h>
#include <t8_forest/t8_forest_geometrical.h>
#include <t8_forest/t8_forest_partition.h>
#include <t8_forest/t8_forest_io.h>
#include <t8_forest/t8_forest_types.h>

/**
 * Callback to perform some arbitrary mesh adaptation within the test suite.
 *
 * Note: The argument list has to be the same as for \ref t8_forest_adapt_t, even
 *       if most arguments are unused. For their meaning, please refer to \ref t8_forest_adapt_t.
 *       (The following doxygen documentation is just to make sure it is technically documented.)
 *
 * \param[in] forest        "forest" argument of \ref t8_forest_adapt_t
 * \param[in] forest_from   "forest_from" argument of \ref t8_forest_adapt_t
 * \param[in] which_tree    "which_tree" argument of \ref t8_forest_adapt_t
 * \param[in] tree_class    "tree_class" argument of \ref t8_forest_adapt_t
 * \param[in] lelement_id   "lelement_id" argument of \ref t8_forest_adapt_t
 * \param[in] scheme        "scheme" argument of \ref t8_forest_adapt_t
 * \param[in] is_family     "is_family" argument of \ref t8_forest_adapt_t
 * \param[in] num_elements  "num_elements" argument of \ref t8_forest_adapt_t
 * \param[in] elements      "elements" argument of \ref t8_forest_adapt_t
 *
 * \return 1 if the element will be refined, 0 otherwise.
*/
int
refine_some_callback ([[maybe_unused]] t8_forest_t forest, [[maybe_unused]] t8_forest_t forest_from,
                      [[maybe_unused]] t8_locidx_t which_tree, [[maybe_unused]] t8_eclass_t tree_class,
                      [[maybe_unused]] t8_locidx_t lelement_id, [[maybe_unused]] const t8_scheme *scheme,
                      [[maybe_unused]] const int is_family, [[maybe_unused]] const int num_elements,
                      [[maybe_unused]] t8_element_t *elements[])
{
  // Refine some elements.
  if (lelement_id % (forest_from->mpirank + 1) == 1) {
    return 1;
  }
  return 0;
}

/**
 * Callback to coarsen all families (that are on the same process).
 *
 * This callback is crucial for testing the PFC functionality.
 *
 * Note: The argument list has to be the same as for \ref t8_forest_adapt_t, even
 *       if most arguments are unused. For their meaning, please refer to \ref t8_forest_adapt_t.
 *       (The following doxygen documentation is just to make sure it is technically documented.)
 * \param[in] forest        "forest" argument of \ref t8_forest_adapt_t
 * \param[in] forest_from   "forest_from" argument of \ref t8_forest_adapt_t
 * \param[in] which_tree    "which_tree" argument of \ref t8_forest_adapt_t
 * \param[in] tree_class    "tree_class" argument of \ref t8_forest_adapt_t
 * \param[in] lelement_id   "lelement_id" argument of \ref t8_forest_adapt_t
 * \param[in] scheme        "scheme" argument of \ref t8_forest_adapt_t
 * \param[in] is_family     "is_family" argument of \ref t8_forest_adapt_t
 * \param[in] num_elements  "num_elements" argument of \ref t8_forest_adapt_t
 * \param[in] elements      "elements" argument of \ref t8_forest_adapt_t
 *
 * \return -1 if the elements are a family, 0 otherwise.
 */
int
coarsen_all_callback ([[maybe_unused]] t8_forest_t forest, [[maybe_unused]] t8_forest_t forest_from,
                      [[maybe_unused]] t8_locidx_t which_tree, [[maybe_unused]] t8_eclass_t tree_class,
                      [[maybe_unused]] t8_locidx_t lelement_id, [[maybe_unused]] const t8_scheme *scheme,
                      const int is_family, [[maybe_unused]] const int num_elements,
                      [[maybe_unused]] t8_element_t *elements[])
{
  if (is_family) {
    return -1;  // coarsen every family
  }
  return 0;
}

/**
 * Class to test the partition-for-coarsening functionality.
*/
class t8_test_partition_for_coarsening_test: public testing::TestWithParam<std::tuple<int, cmesh_example_base *>> {

 protected:
  /** During SetUp, set the scheme and the eclass based on the current testing parameters.*/
  void
  SetUp () override
  {
    // Get scheme.
    const int scheme_id = std::get<0> (GetParam ());
    scheme = create_from_scheme_id (scheme_id);

    // Construct cmesh and store name.
    cmesh = std::get<1> (GetParam ())->cmesh_create ();
    cmesh_name = std::get<1> (GetParam ())->name;

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

  /**
   * Helper function to create vtk output for debugging purposes.
   * Only writes output if the option is manually activated by changing
   * write_vtk_output_flag from "false" to "true" in the code below.
   *
   * \param[in] forest      the forest to be written to vtk
   * \param[in] forest_name its name used for the output file name
   *
   */
  void
  write_forest_to_vtk_if_flag_set (t8_forest_t forest, std::string forest_name)
  {

    // Manually change flag to "true" here to write vtk output.
    bool write_vtk_output_flag = false;

    // Skip debug vtk output for bigmeshes.
    if (cmesh_name.find (std::string ("bigmesh")) != std::string::npos) {
      return;
    }

    // If flag is true, write vtk output.
    if (write_vtk_output_flag) {
      std::string fileName = cmesh_name;
      fileName = (fileName.append ("_")).append (forest_name.c_str ());
      t8_forest_write_vtk (forest, fileName.c_str ());
    }
  }

  // Member variables: The currently tested scheme and eclass.
  const t8_scheme *scheme; /**< The currently tested scheme. */
  t8_cmesh_t cmesh;        /**< The currently tested cmesh. */
  std::string cmesh_name;  /**< Name of the currently tested cmesh.*/
};

// The test's main function.
TEST_P (t8_test_partition_for_coarsening_test, test_partition_for_coarsening)
{

  t8_global_productionf ("##############################################################################\n");
  t8_global_productionf ("######## Testing mesh: %s\n", cmesh_name.c_str ());
  t8_global_productionf ("##############################################################################\n");
  sc_MPI_Barrier (sc_MPI_COMM_WORLD);

  // -------------------------------------------
  // ----- (1.) Create uniform base forest -----
  // -------------------------------------------
  t8_global_productionf ("Create uniform forest.\n");

  // Initial uniform refinement level
  const int level = 2;

  // Increase reference counters of cmesh and scheme to avoid reaching zero.
  t8_cmesh_ref (cmesh);
  scheme->ref ();

  // Create initial, uniform base forest.
  t8_forest_t uniform_forest = t8_forest_new_uniform (cmesh, scheme, level, 0, sc_MPI_COMM_WORLD);

  // If the associated flag is set, write forest to vtk.
  write_forest_to_vtk_if_flag_set (uniform_forest, "uniform_forest");

  // -----------------------------------------------------------
  // ----- (2.) Create adapted forest with some refinement -----
  // -----------------------------------------------------------
  t8_global_productionf ("Adapt uniform forest.\n");

  // Create adapted base forest.
  t8_forest_t adapted_base_forest = t8_forest_new_adapt (uniform_forest, refine_some_callback, 0, 0, nullptr);

  // If the associated flag is set, write forest to vtk.
  write_forest_to_vtk_if_flag_set (adapted_base_forest, "adapted_base_forest");

  // ---------------------------------------
  // ----- (3.) Apply PFC partitioning -----
  // ---------------------------------------
  t8_global_productionf ("Partition forest with PFC.\n");

  // Initialize forest.
  t8_forest_t pfc_forest;
  t8_forest_init (&pfc_forest);

  // Set partitioning mode and source.
  // Note: The third argument being nonzero activates the PFC correction.
  t8_forest_set_partition (pfc_forest, adapted_base_forest, 1);

  // Increase reference counter of base forest to avoid destruction.
  t8_forest_ref (adapted_base_forest);

  // Commit PFC-partitioned forest.
  t8_forest_commit (pfc_forest);

  // If the associated flag is set, write forest to vtk.
  write_forest_to_vtk_if_flag_set (pfc_forest, "pfc_forest");

  // -----------------------------------------------------
  // ----- (4.) Coarsen all families on same process -----
  // -----------------------------------------------------
  t8_global_productionf ("Coarsen all families in PFC forest.\n");

  // Create coarsened PFC forest.
  t8_forest_t coarsened_pfc_forest = t8_forest_new_adapt (pfc_forest, coarsen_all_callback, 0, 0, nullptr);

  // If the associated flag is set, write forest to vtk.
  write_forest_to_vtk_if_flag_set (coarsened_pfc_forest, "coarsened_pfc_forest");

  // ---------------------------------------
  // ----- (5.) Checks for correctness -----
  // ---------------------------------------
  t8_global_productionf ("Start preparing checks.\n");

  // (5a.) Repartition coarsened PFC forest:
  // ---------------------------------------
  t8_global_productionf ("Repartition coarsened forest without PFC option.\n");

  // Initialize forest
  t8_forest_t repartitioned_coarse_pfc_forest;
  t8_forest_init (&repartitioned_coarse_pfc_forest);

  // Set partition mode and source.
  t8_forest_set_partition (repartitioned_coarse_pfc_forest, coarsened_pfc_forest, 0);

  // Commit the forest.
  t8_forest_commit (repartitioned_coarse_pfc_forest);

  // (5b.) Take different route to allow comparison:
  // -----------------------------------------------
  // Steps:
  // * Gather base base forest on one rank
  // * Coarsen all families
  // * Repartition without PFC
  t8_global_productionf ("Create the same forest via different route.\n");

  // Gather base forest on rank zero.
  t8_forest_t forest_gathered = t8_forest_new_gather (adapted_base_forest, 0);

  // Coarsen all families.
  t8_forest_t coarse_gathered_forest = t8_forest_new_adapt (forest_gathered, coarsen_all_callback, 0, 0, nullptr);

  // Partition without PFC option.
  t8_forest_t repartitioned_coarse_gathered_forest;
  t8_forest_init (&repartitioned_coarse_gathered_forest);
  t8_forest_set_partition (repartitioned_coarse_gathered_forest, coarse_gathered_forest, 0);
  t8_forest_commit (repartitioned_coarse_gathered_forest);

  // (5c.) Verify that the two forests are equal:
  // --------------------------------------------
  t8_global_productionf ("Verify that the two coarsened and repartitioned forests match!\n");
  EXPECT_TRUE (t8_forest_is_equal (repartitioned_coarse_pfc_forest, repartitioned_coarse_gathered_forest));

  // --------------------------------
  // ----- (6.) Clean up memory -----
  // --------------------------------
  t8_forest_unref (&repartitioned_coarse_pfc_forest);
  t8_forest_unref (&repartitioned_coarse_gathered_forest);
}

// Instantiate parameterized test to be run for all schemes.
INSTANTIATE_TEST_SUITE_P (t8_gtest_partition_for_coarsening, t8_test_partition_for_coarsening_test,
                          testing::Combine (AllSchemeCollections, AllCmeshsParam), pretty_print_base_example_scheme);
