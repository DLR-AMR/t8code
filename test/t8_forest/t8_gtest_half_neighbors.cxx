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

#include <gtest/gtest.h>
#include <test/t8_gtest_custom_assertion.hxx>

#include <t8_eclass.h>
#include <t8_cmesh.h>
#include <t8_forest/t8_forest_general.h>
#include <t8_forest/t8_forest_types.h>
#include <t8_schemes/t8_default/t8_default_cxx.hxx>
#include <t8_cmesh/t8_cmesh_offset.h>
#include <t8_cmesh/t8_cmesh_examples.h>
#include <t8_forest/t8_forest_partition.h>
#include <t8_forest/t8_forest_private.h>
#include <test/t8_gtest_macros.hxx>

class forest_half_neighbors: public testing::TestWithParam<std::tuple<t8_eclass, int>> {
 protected:
  void
  SetUp () override
  {
    eclass = std::get<0> (GetParam ());
    cmesh_type = std::get<1> (GetParam ());

    default_scheme = t8_scheme_new_default_cxx ();
    /* Construct a coarse mesh of one tree */
    cmesh = t8_cmesh_new_from_class (eclass, sc_MPI_COMM_WORLD);
  }

  t8_eclass_t eclass;
  int cmesh_type;
  t8_cmesh_t cmesh;
  t8_scheme_cxx_t *default_scheme;
  t8_element_t *neighbor;
};

#if 0
/* Depending on an integer i create a different cmesh.
 * i = 0: cmesh_new_class
 * i = 1: cmesh_new_hypercube
 * i = 2: cmesh_new_bigmesh (100 trees)
 * else:  cmesh_new_class
 */
TEST_P (forest_halt_neighbors, test_create_cmesh)
{
  switch (cmesh_type) {
  case 0:
    return t8_cmesh_new_from_class (eclass, sc_MPI_COMM_WORLD);
  case 1:
    return t8_cmesh_new_hypercube (eclass, sc_MPI_COMM_WORLD, 0, 0, 0);
  case 2:
    return t8_cmesh_new_bigmesh (eclass, 100, sc_MPI_COMM_WORLD);
  default:
    return t8_cmesh_new_from_class (eclass, sc_MPI_COMM_WORLD);
  }
}
#endif

TEST_P (forest_half_neighbors, test_half_neighbors)
{

  const int level = 3;
  sc_array_t owners;
  int dual_face;

  t8_debugf ("Testing half neighbors with eclass %s, cmesh type %i.\n", t8_eclass_to_string[eclass], cmesh_type);

  /* initialize the array of owners to store ints */
  sc_array_init (&owners, sizeof (int));
  /* Build a uniform forest */
  t8_forest_t forest = t8_forest_new_uniform (cmesh, default_scheme, level, 0, sc_MPI_COMM_WORLD);
  t8_eclass_scheme_c *ts = t8_forest_get_eclass_scheme (forest, eclass);
  /* iterate over all elements */
  for (t8_locidx_t itree = 0; itree < t8_forest_get_num_local_trees (forest); itree++) {
    for (t8_locidx_t ielement = 0; ielement < t8_forest_get_tree_num_elements (forest, itree); ielement++) {
      const t8_element_t *element = t8_forest_get_element_in_tree (forest, itree, ielement);
      /* iterate over the faces */
      for (int face = 0; face < ts->t8_element_num_faces (element); face++) {
        /* Get the eclass of the face neighbor and get the scheme */
        const t8_eclass_t neigh_class = t8_forest_element_neighbor_eclass (forest, itree, element, face);
        t8_eclass_scheme_c *neigh_scheme = t8_forest_get_eclass_scheme (forest, neigh_class);
        const int num_face_neighs = ts->t8_element_num_face_children (element, face);
        t8_element_t **half_neighbors = T8_ALLOC (t8_element_t *, num_face_neighs);
        ts->t8_element_new (num_face_neighs, half_neighbors);
        t8_forest_element_half_face_neighbors (forest, itree, element, half_neighbors, neigh_scheme, face,
                                               num_face_neighs, NULL);
        /* allocate memory for element's neighbor and construct it */
        neigh_scheme->t8_element_new (1, &neighbor);
        const t8_locidx_t neigh_tree
          = t8_forest_element_face_neighbor (forest, itree, element, neighbor, neigh_scheme, face, &dual_face);
        if (neigh_tree > 0) {
          /* We now check whether the face children of neighbor are the half neighbors. */
          T8_ASSERT (num_face_neighs == neigh_scheme->t8_element_num_face_children (neighbor, dual_face));
          t8_element_t **neighbor_face_children = T8_ALLOC (t8_element_t *, num_face_neighs);
          neigh_scheme->t8_element_new (num_face_neighs, neighbor_face_children);
          int *child_ids = T8_ALLOC (int, num_face_neighs);
          neigh_scheme->t8_element_children_at_face (neighbor, dual_face, neighbor_face_children, num_face_neighs,
                                                     child_ids);
          /* Check that the children at face of the neighbor are the half neighbors of the element */
          for (int ineigh = 0; ineigh < num_face_neighs; ineigh++) {
            EXPECT_ELEM_EQ (neigh_scheme, neighbor_face_children[ineigh], half_neighbors[ineigh])
              << "ineigh = " << ineigh << " face = " << face;
          }
          neigh_scheme->t8_element_destroy (num_face_neighs, neighbor_face_children);
          T8_FREE (child_ids);
          T8_FREE (neighbor_face_children);
        }
        neigh_scheme->t8_element_destroy (1, &neighbor);
        neigh_scheme->t8_element_destroy (num_face_neighs, half_neighbors);
        T8_FREE (half_neighbors);
      }
    }
  }
  t8_forest_unref (&forest);
  sc_array_reset (&owners);
}

INSTANTIATE_TEST_SUITE_P (t8_gtest_half_neighbors, forest_half_neighbors,
                          testing::Combine (AllEclasses, testing::Range (0, 3)));
