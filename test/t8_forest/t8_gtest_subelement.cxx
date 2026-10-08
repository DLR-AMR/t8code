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
#include <test/t8_gtest_custom_assertion.hxx>

#include <t8.h>
#include <t8_cmesh/t8_cmesh.h>
#include <t8_cmesh/t8_cmesh_examples.h>
#include <t8_forest/t8_forest_io.h>
#include <t8_forest/t8_forest_general.h>
#include <t8_forest/t8_forest_subelement.hxx>
#include <t8_schemes/t8_subelement/t8_subelement.hxx>
#include <t8_forest/t8_forest_geometrical.hxx>
#include <t8_schemes/t8_scheme.hxx>
#include <t8_types/t8_vec.hxx>
#include <cmath>
#include <ranges>
#include <vector>

/** Check that the hanging node resolution for 2D hybrid meshes works. At the moment we only check the functionality 
* needed for visualization (so e.g. no connectivity). 
*/
TEST (t8_gtest_subelement, hybrid_hanging_nodes_visualization)
{
  /* Setup: Build hypercube cmesh and uniform forest with the subelement scheme. */
  const int level = 2;
  t8_cmesh_t cmesh;
  t8_cmesh_init (&cmesh);
  t8_cmesh_new_2D_hypercube_hybrid (cmesh, sc_MPI_COMM_WORLD);

  t8_forest_t forest = t8_forest_new_uniform (cmesh, t8_scheme_new_subelement (), level, 0, sc_MPI_COMM_WORLD);

  /* Initial uniform forest should not have any subelements. */
  EXPECT_FALSE (t8_forest_has_subelements (forest));

  /* Adapt the forest (refining every second element). */
  forest = t8_forest_new_adapt (forest, refine_every_nth_element_callback<2>, 0, 0, NULL);

  /* Before resolving hanging nodes, subelements should not yet be introduced. */
  EXPECT_FALSE (t8_forest_has_subelements (forest));
  EXPECT_FALSE (t8_forest_is_conforming (forest));

  // Check that discarding without subelements just does nothing.
  // Introduce second forest for comparisons.
  t8_forest_ref (forest);
  auto forest_compare = t8_forest_discard_subelements (forest);
  EXPECT_FALSE (t8_forest_has_subelements (forest_compare));
  EXPECT_FOREST_EQ (forest, forest_compare);

  /* Remove hanging nodes by inserting subelements. The forest is already balanced as we only adapted once. */
  forest = t8_forest_remove_hanging_nodes (forest);
  EXPECT_TRUE (t8_forest_is_committed (forest));

  /* Hanging node resolution must introduce subelements into the forest. */
  EXPECT_TRUE (t8_forest_has_subelements (forest));

  /* Adding transition subelements must increase the total leaf count. */
  const t8_gloidx_t num_leaves_sub = t8_forest_get_global_num_leaf_elements (forest);
  const t8_gloidx_t num_leaves_adapted = t8_forest_get_global_num_leaf_elements (forest_compare);
  EXPECT_GT (num_leaves_sub, num_leaves_adapted);

  /* Repartition the forest containing subelements (exercises MPI_Pack / MPI_Unpack). */
  t8_forest_t forest_partitioned;
  t8_forest_init (&forest_partitioned);
  t8_forest_set_partition (forest_partitioned, forest, true);
  t8_forest_commit (forest_partitioned);

  /* Subelements and leaf count must remain consistent after repartitioning. */
  EXPECT_TRUE (t8_forest_has_subelements (forest_partitioned));
  EXPECT_EQ (t8_forest_get_global_num_leaf_elements (forest_partitioned), num_leaves_sub);

#if T8_ENABLE_DEBUG
  /* Write vtk file in debug mode. This checks that all functions are available that are required for visualization. */
  t8_forest_write_vtk (forest_partitioned, "test_subelements");
#endif

  /* Discard subelements from the partitioned forest. */
  forest = t8_forest_discard_subelements (forest_partitioned);
  /* Subelements should now be completely removed. */
  EXPECT_FALSE (t8_forest_has_subelements (forest));
  /* Discarding subelements should restore the pre-resolution leaf count. */
  EXPECT_EQ (t8_forest_get_global_num_leaf_elements (forest), num_leaves_adapted);

  /* Clean up. */
  t8_forest_unref (&forest);
  t8_forest_unref (&forest_compare);
}

/* Check that converting a batch of reference coordinates gives the same result as converting each point on its own,
 * for regular elements as well as for subelements. */
TEST (t8_gtest_subelement, batch_reference_coords)
{
  /* Setup: Build a forest with subelements as above. */
  const int level = 2;
  t8_cmesh_t cmesh;
  t8_cmesh_init (&cmesh);
  t8_cmesh_new_2D_hypercube_hybrid (cmesh, sc_MPI_COMM_WORLD);
  t8_forest_t forest = t8_forest_new_uniform (cmesh, t8_scheme_new_subelement (), level, 0, sc_MPI_COMM_WORLD);
  forest = t8_forest_new_adapt (forest, refine_every_nth_element_callback<2>, 0, 0, NULL);
  forest = t8_forest_remove_hanging_nodes (forest);
  ASSERT_TRUE (t8_forest_has_subelements (forest));

  /** Points in the reference space of triangles and quads. The elements are 2D and ignore the third component,
   * hence we set it to an arbitrary value. */
  const std::vector<t8_3D_vec> ref_coords = { { 0, 0, 0.3 }, { 1, 0, 0.3 }, { 1, 1, 0.3 }, { 0.75, 0.25, 0.3 } };
  const t8_scheme *scheme = t8_forest_get_scheme (forest);
  for (t8_locidx_t itree = 0; itree < t8_forest_get_num_local_trees (forest); ++itree) {
    const t8_eclass_t tree_class = t8_forest_get_tree_class (forest, itree);
    for (t8_locidx_t ielement = 0; ielement < t8_forest_get_tree_num_leaf_elements (forest, itree); ++ielement) {
      const t8_element_t *element = t8_forest_get_leaf_element_in_tree (forest, itree, ielement);
      /** The points in the reference space of the tree, converted as one batch. */
      std::vector<t8_3D_vec> tree_ref_coords (ref_coords.size ());
      scheme->element_get_reference_coords (tree_class, element, ref_coords, tree_ref_coords);
      /** The points in the domain, evaluated as one batch. */
      std::vector<t8_3D_vec> domain_coords (ref_coords.size ());
      t8_forest_element_from_ref_coords (forest, itree, element, ref_coords, domain_coords);
      for (size_t ipoint = 0; ipoint < ref_coords.size (); ++ipoint) {
        const t8_3D_vec single_tree_ref_coords
          = scheme->element_get_reference_coords (tree_class, element, ref_coords[ipoint]);
        /* Only the first two components are part of the result, since the trees are 2D. */
        EXPECT_TRUE (t8_eq (tree_ref_coords[ipoint] | std::views::take (2),
                            single_tree_ref_coords | std::views::take (2), T8_PRECISION_EPS));
#if T8_ENABLE_DEBUG
        /* In debug mode, the third component is NaN. */
        EXPECT_TRUE (std::isnan (tree_ref_coords[ipoint][2]));
#endif
        const t8_3D_vec single_domain_coords
          = t8_forest_element_from_ref_coords (forest, itree, element, ref_coords[ipoint]);
        EXPECT_TRUE (t8_eq (domain_coords[ipoint], single_domain_coords, T8_PRECISION_EPS));
      }
    }
  }
  t8_forest_unref (&forest);
}
