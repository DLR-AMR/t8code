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
#include <test/t8_gtest_memory_macros.hxx>

#include <t8_eclass/t8_eclass.h>
#include <t8_cmesh/t8_cmesh.h>
#include <t8_forest/t8_forest_general.h>
#include <t8_forest/t8_forest_types.h>
#include <t8_schemes/t8_default/t8_default.hxx>
#include <t8_cmesh/t8_cmesh_internal/t8_cmesh_offset.h>
#include <t8_cmesh/t8_cmesh_internal/t8_cmesh_trees.h>
#include <t8_cmesh/t8_cmesh_examples.h>
#include <t8_forest/t8_forest_partition.h>
#include <t8_forest/t8_forest_private.h>
#include <test/t8_gtest_schemes.hxx>

struct element_neighbor_eclass: public testing::TestWithParam<int>
{
 protected:
  void
  SetUp () override
  {
    const int level = 1;
    const int scheme_id = GetParam ();
    scheme = create_from_scheme_id (scheme_id);
    
    // Construct a hybrid coarse mesh
    t8_cmesh_t cmesh = t8_cmesh_new_full_hybrid (sc_MPI_COMM_WORLD);

    // Build a uniform forest
    forest = t8_forest_new_uniform (cmesh, scheme, level, 0, sc_MPI_COMM_WORLD);

    // Go through all the trees in the forest and store their eclass
    tree_eclass.resize( t8_forest_get_num_local_trees (forest) );
    for (t8_locidx_t itree = 0; itree < t8_forest_get_num_local_trees (forest); itree++) {
      tree_eclass[itree] = t8_forest_get_tree_class (forest, itree);
      std::cerr << "tree " << itree << "  eclass " << tree_eclass[itree] << std::endl;
    }
  }

  void
  TearDown () override
  {
    t8_forest_unref (&forest);
  }

  t8_forest_t forest;
  const t8_scheme *scheme;
  std::vector<t8_eclass_t> tree_eclass;
};

TEST_P (element_neighbor_eclass, test_half_neighbors)
{
  // TODO: print scheme?
  //t8_debugf ("Testing element neighbor eclass with eclass %s.\n", t8_eclass_to_string[eclass]);

  // Iterate over all trees
  // TODO: ghost trees?
  for (t8_locidx_t itree = 0; itree < t8_forest_get_num_local_trees (forest); itree++) {
    const t8_eclass_t tree_class = t8_forest_get_eclass (forest, itree);

    // Iterate over all elements
    for (t8_locidx_t ielement = 0; ielement < t8_forest_get_tree_num_leaf_elements (forest, itree); ielement++) {
      const t8_element_t *element = t8_forest_get_leaf_element_in_tree (forest, itree, ielement);

      // Iterate over all faces
      for (int iface = 0; iface < scheme->element_get_num_faces (tree_class, element); iface++) {

        // Get neighbor tree id
        t8_cmesh_t cmesh = t8_forest_get_cmesh (forest);
        t8_locidx_t const ltreeid_in_cmesh = t8_forest_ltreeid_to_cmesh_ltreeid (forest, itree);
        t8_locidx_t const cmesh_dual_itree = t8_cmesh_get_face_neighbor (cmesh, ltreeid_in_cmesh, iface, nullptr, nullptr);
        t8_locidx_t const forest_dual_itree = t8_forest_cmesh_ltreeid_to_ltreeid(forest, cmesh_dual_itree);

        // Get the eclass of the face neighbor
        const t8_eclass_t neigh_class = t8_forest_element_neighbor_eclass (forest, itree, element, iface);

        if (forest_dual_itree > -1) {
          // Compare
          EXPECT_EQ (neigh_class, tree_eclass[forest_dual_itree]);
        }
        else {
          // TODO: ?
        }
      }
    }
  }
}

INSTANTIATE_TEST_SUITE_P (t8_gtest_element_neighbor_eclass, element_neighbor_eclass, AllSchemeCollections);
