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

/** \file t8_gtest_subelement.cxx
 * Unit test for the subelement scheme. This test checks the complete pipeline of hanging node resolution for
 *´hybrid 2D meshes. Currently, we only test that the functions required for visualization work correctly.
 */

#include <gtest/gtest.h>
#include <test/t8_gtest_adapt_callbacks.hxx>

#include <t8.h>
#include <t8_cmesh/t8_cmesh.h>
#include <t8_cmesh/t8_cmesh_examples.h>
#include <t8_forest/t8_forest_general.h>
#include <t8_forest/t8_forest_subelement.hxx>
#include <t8_schemes/t8_subelement/t8_subelement.hxx>

/** Check that the hanging node resolution for 2D hybrid meshes works. At the moment we only check the functionality 
* needed for visualization (so e.g. no connectivity). 
*/
TEST (t8_gtest_subelement, hybrid_hanging_nodes_visualization)
{ 
  /* Setup: Build hypercube cmesh and uniform forest with the subelement scheme. */
  const int level = 3;
  t8_cmesh_t cmesh;
  t8_cmesh_init (&cmesh);
  t8_cmesh_new_2D_hypercube_hybrid (cmesh, sc_MPI_COMM_WORLD);

  t8_forest_t forest = t8_forest_new_uniform (cmesh, t8_scheme_new_subelement (), level, 0, sc_MPI_COMM_WORLD);

  /* Initial uniform forest should not have any subelements. */
  EXPECT_FALSE (t8_forest_has_global_subelements (forest));

  /* Adapt the forest (refining every second element). */
  forest = t8_forest_new_adapt (forest, t8_test_adapt_even_global_id, 0, 0, NULL);

  /* Before resolving hanging nodes, subelements should not yet be introduced. */
  EXPECT_FALSE (t8_forest_has_global_subelements (forest));
  const t8_gloidx_t num_leaves_adapted = t8_forest_get_global_num_leaf_elements (forest);

  /* Remove hanging nodes by inserting subelements. The forest is already balanced as we only adapted once. */
  forest = t8_forest_remove_hanging_nodes (forest);
  EXPECT_TRUE(t8_forest_is_committed (forest));

  /* Hanging node resolution must introduce subelements into the forest. */
  EXPECT_TRUE (t8_forest_has_global_subelements (forest));

  /* Adding transition subelements must increase (or equal) the total leaf count. */
  const t8_gloidx_t num_leaves_sub = t8_forest_get_global_num_leaf_elements (forest);
  EXPECT_GT (num_leaves_sub, num_leaves_adapted);

  /* Repartition the forest containing subelements (exercises MPI_Pack / MPI_Unpack). */
  t8_forest_t forest_partitioned;
  t8_forest_init (&forest_partitioned);
  t8_forest_set_partition (forest_partitioned, forest, 0);
  t8_forest_commit (forest_partitioned);

  /* Subelements and leaf count must remain consistent after repartitioning. */
  EXPECT_TRUE (t8_forest_has_global_subelements  (forest_partitioned));
  EXPECT_EQ (t8_forest_get_global_num_leaf_elements (forest_partitioned), num_leaves_sub);

  /* Discard subelements from the partitioned forest. */
  forest = t8_forest_discard_subelements (forest_partitioned);
  /* Subelements should now be completely removed. */
  EXPECT_FALSE (t8_forest_has_global_subelements  (forest));
  /* Discarding subelements should restore the pre-resolution leaf count. */
  EXPECT_EQ (t8_forest_get_global_num_leaf_elements (forest), num_leaves_adapted);

  /* Clean up forest */
  t8_forest_unref (&forest);
}