/*
  This file is part of t8code.
  t8code is a C library to manage a collection (a forest) of multiple
  connected adaptive space-trees of general element classes in parallel.

  Copyright (C) 2023 the developers

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

#include <sc/src/sc_functions.h>
#include <gtest/gtest.h>
#include <t8_eclass.h>
#include <t8_schemes/t8_default/t8_default.hxx>
#include <t8_schemes/t8_default/t8_default_pyramid/t8_dpyramid_bits.h>
#include <t8_cmesh/t8_cmesh_examples.h>
#include <t8_forest/t8_forest_general.h>
#include <t8_forest/t8_forest_geometrical.h>
#include <test/t8_gtest_macros.hxx>

/**
 * This file tests the face normal computation of elements.
 */

class class_forest_face_normal: public testing::TestWithParam<std::tuple<t8_eclass_t, int>> {
 protected:
  void
  SetUp () override
  {
    eclass = std::get<0> (GetParam ());
    level = std::get<1> (GetParam ());
    scheme = t8_scheme_new_default ();
    t8_cmesh_t cmesh = t8_cmesh_new_hypercube (eclass, sc_MPI_COMM_WORLD, 0, 0, 0);
    const int do_face_ghost = 1;
    forest = t8_forest_new_uniform (cmesh, scheme, level, do_face_ghost, sc_MPI_COMM_WORLD);
  }
  void
  TearDown () override
  {
    t8_forest_unref (&forest);
  }
  t8_forest_t forest;
  t8_scheme *scheme;
  t8_eclass_t eclass;
  int level;
};

TEST_P (class_forest_face_normal, back_and_forth)
{
  /** Iterate over all elements of a uniformly refined forest. For all faceneighbors elements, if they are local,
    * check if their facenormal is the negative of the corresponding facenormal of the neighbor elements.
    */

  const t8_scheme *scheme = t8_forest_get_scheme (forest);
  const t8_locidx_t local_num_trees = t8_forest_get_num_local_trees (forest);
  /* Iterate over all elements. */
  for (t8_locidx_t itree = 0; itree < local_num_trees; itree++) {
    const t8_locidx_t tree_elements = t8_forest_get_tree_num_elements (forest, itree);
    const t8_eclass_t tree_class = t8_forest_get_tree_class (forest, itree);
    ASSERT_EQ (eclass, tree_class);
    for (t8_locidx_t ielement = 0; ielement < tree_elements; ielement++) {
      const t8_element_t *element = t8_forest_get_element_in_tree (forest, itree, ielement);
      const int num_faces = scheme->element_get_num_faces (tree_class, element);
      for (int iface = 0; iface < num_faces; iface++) {
        /* Compute facenormal */
        double face_normal[3] = { 0, 0, 0 };
        t8_forest_element_face_normal (forest, itree, element, iface, face_normal);

        /* Get all face neighbors */

        t8_element_t **neighbors;
        int num_neighbors;
        const int forest_is_balanced = 1;
        t8_eclass_t neigh_eclass;
        int *dual_faces;
        t8_locidx_t *neigh_ids;

        t8_gloidx_t gneightree;
        t8_forest_leaf_face_neighbors_ext (forest, itree, element, &neighbors, iface, &dual_faces, &num_neighbors,
                                           &neigh_ids, &neigh_eclass, forest_is_balanced, &gneightree, NULL);

        /* Iterate and compute their facenormal */
        for (int ineigh = 0; ineigh < num_neighbors; ineigh++) {
          t8_locidx_t lneightree = t8_forest_get_local_or_ghost_id (forest, gneightree);
          double neigh_face_normal[3];
          t8_forest_element_face_normal (forest, lneightree, neighbors[ineigh], dual_faces[ineigh], neigh_face_normal);
          EXPECT_NEAR (face_normal[0], -neigh_face_normal[0], T8_PRECISION_SQRT_EPS);
          EXPECT_NEAR (face_normal[1], -neigh_face_normal[1], T8_PRECISION_SQRT_EPS);
          EXPECT_NEAR (face_normal[2], -neigh_face_normal[2], T8_PRECISION_SQRT_EPS);
        }

        if (num_neighbors > 0) {
          scheme->element_destroy (neigh_eclass, num_neighbors, neighbors);
          T8_FREE (neigh_ids);
          T8_FREE (neighbors);
          T8_FREE (dual_faces);
        }
      }
    }
  }
}

INSTANTIATE_TEST_SUITE_P (t8_gtest_forest_face_normal, class_forest_face_normal,
                          testing::Combine (AllEclasses, testing::Range (0, 2)));
