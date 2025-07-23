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
#include <t8_eclass.h>
#include <t8_schemes/t8_default/t8_default.hxx>
#include <t8_forest/t8_forest_general.h>
#include <t8_forest/t8_forest_ghost/t8_forest_ghost.h>
#include <t8_forest/t8_forest_private.h>
#include <t8_cmesh.h>
#include "test/t8_cmesh_generator/t8_cmesh_example_sets.hxx"
#include <test/t8_gtest_macros.hxx>
#include <test/t8_gtest_schemes.hxx>

/* This test program tests the forest ghost exchange routine.
 * Given a forest for which the ghost layer was created and an array
 * storing data for the local elements and the ghost elements, ghost_exchange
 * communicates the data of the local elements to the ghost entries of the
 * processes for which these elements are ghost.
 * We test the ghost exchange routine for several forests on different
 * coarse meshes.
 * One test is an integer entry '42' for each element,
 * in a second test, we store the element's linear id in the data array.
 */

class forest_ghost_exchange: public testing::TestWithParam<std::tuple<int, cmesh_example_base *>> {
 protected:
  void
  SetUp () override
  {
    const int scheme_id = std::get<0> (GetParam ());
    scheme = create_from_scheme_id (scheme_id);
    /* Construct a cmesh */
    cmesh = std::get<1> (GetParam ())->cmesh_create ();
    if (t8_cmesh_is_empty (cmesh)) {
      /* empty cmeshes are currently not supported */
      GTEST_SKIP ();
    }
  }
  void
  TearDown () override
  {
    t8_cmesh_destroy (&cmesh);
    scheme->unref ();
  }
  const t8_scheme *scheme;
  t8_cmesh_t cmesh;
};

static int
t8_test_exchange_adapt (t8_forest_t forest, [[maybe_unused]] t8_forest_t forest_from,
                        [[maybe_unused]] t8_locidx_t which_tree, const t8_eclass_t tree_class,
                        [[maybe_unused]] t8_locidx_t lelement_id, const t8_scheme *scheme,
                        [[maybe_unused]] const int is_family, [[maybe_unused]] const int num_elements,
                        t8_element_t *elements[])
{
  /* refine every second element up to the maximum level */
  const int level = scheme->element_get_level (tree_class, elements[0]);
  const t8_linearidx_t eid = scheme->element_get_linear_id (tree_class, elements[0], level);
  const int maxlevel = *(int *) t8_forest_get_user_data (forest);

  if (eid % 2 && level < maxlevel) {
    return 1;
  }
  return 0;
}

/* Construct a data array of uin64_t for all elements and all ghosts,
 * fill the element's entries with their linear id, perform the ghost exchange and
 * check whether the ghost's entries are their linear id.
 */
static void
t8_test_ghost_exchange_data_id (t8_forest_t forest)
{
  const t8_scheme *scheme = t8_forest_get_scheme (forest);
  size_t array_pos = 0;
  sc_array_t element_data;

  t8_locidx_t num_elements = t8_forest_get_local_num_leaf_elements (forest);
  t8_locidx_t num_ghosts = t8_forest_get_num_ghosts (forest);
  /* Allocate a uin64_t as data for each element and each ghost */
  sc_array_init_size (&element_data, sizeof (t8_linearidx_t), num_elements + num_ghosts);

  /* Fill the local element entries with their linear id */
  for (t8_locidx_t itree = 0; itree < t8_forest_get_num_local_trees (forest); itree++) {
    const t8_eclass_t tree_class = t8_forest_get_tree_class (forest, itree);
    for (t8_locidx_t ielem = 0; ielem < t8_forest_get_tree_num_leaf_elements (forest, itree); ielem++) {
      /* Get a pointer to this element */
      const t8_element_t *elem = t8_forest_get_leaf_element_in_tree (forest, itree, ielem);
      const int level = scheme->element_get_level (tree_class, elem);
      /* Compute the linear id of this element */
      const t8_linearidx_t elem_id = scheme->element_get_linear_id (tree_class, elem, level);
      /* Store this id at the element's index in the array */
      *(t8_linearidx_t *) sc_array_index (&element_data, array_pos) = elem_id;
      array_pos++;
    }
  }

  /* Perform the data exchange */
  t8_forest_ghost_exchange_data (forest, &element_data);

  /* We now iterate over all ghost elements and check whether the correct
   * id was received */
  for (t8_locidx_t itree = 0; itree < t8_forest_get_num_ghost_trees (forest); itree++) {
    const t8_eclass_t tree_class = t8_forest_ghost_get_tree_class (forest, itree);
    for (t8_locidx_t ielem = 0; ielem < t8_forest_ghost_tree_num_leaf_elements (forest, itree); ielem++) {
      /* Get a pointer to this ghost */
      const t8_element_t *elem = t8_forest_ghost_get_leaf_element (forest, itree, ielem);
      /* Compute its ghost_id */
      const t8_linearidx_t ghost_id
        = scheme->element_get_linear_id (tree_class, elem, scheme->element_get_level (tree_class, elem));
      /* Compare this id with the entry in the element_data array */
      const t8_linearidx_t ghost_entry = *(t8_linearidx_t *) sc_array_index (&element_data, array_pos);
      ASSERT_EQ (ghost_id, ghost_entry) << "Error when exchanging ghost data. Received wrong element id.\n";
      /* Since array pos ended with the last element in the loop above, we can
       * continue counting for the ghost elements */
      array_pos++;
    }
  }
  /* clean-up */
  sc_array_reset (&element_data);
}

/* Construct a data array of ints for all elements and all ghosts,
 * fill the element's entries with '42', perform the ghost exchange and
 * check whether the ghost's entries are '42'.
 */
static void
t8_test_ghost_exchange_data_int (t8_forest_t forest)
{
  sc_array_t element_data;

  t8_locidx_t num_elements = t8_forest_get_local_num_leaf_elements (forest);
  t8_locidx_t num_ghosts = t8_forest_get_num_ghosts (forest);
  /* Allocate an integer as data for each element and each ghost */
  sc_array_init_size (&element_data, sizeof (int), num_elements + num_ghosts);

  /* Fill the local element entries with the integer 42 */
  for (t8_locidx_t ielem = 0; ielem < num_elements; ielem++) {
    *(int *) t8_sc_array_index_locidx (&element_data, ielem) = 42;
  }
  /* Perform the ghost data exchange */
  t8_forest_ghost_exchange_data (forest, &element_data);

  /* Check for the ghosts that we received the correct data */
  for (t8_locidx_t ielem = 0; ielem < num_ghosts; ielem++) {
    /* Get the integer for this ghost */
    int ghost_int = *(int *) t8_sc_array_index_locidx (&element_data, num_elements + ielem);
    ASSERT_EQ (ghost_int, 42) << "Error when exchanging ghost data. Received wrong data.\n";
  }
  /* clean-up */
  sc_array_reset (&element_data);
}

TEST_P (forest_ghost_exchange, test_ghost_exchange)
{

  /* Compute the minimum level, such that the forest is nonempty */
  int min_level = t8_forest_min_nonempty_level (cmesh, scheme);
  /* we start with an empty level */
  min_level = SC_MAX (min_level - 1, 0);
#if T8_TEST_LEVEL_INT >= 2
  const int max_level = min_level + 2;
#else
  const int max_level = min_level + 3;
#endif
  for (int level = min_level; level < max_level; level++) {
    /* ref the scheme since we reuse it */
    scheme->ref ();
    /* ref the cmesh since we reuse it */
    t8_cmesh_ref (cmesh);
    /* Create a uniformly refined forest */
    t8_forest_t forest = t8_forest_new_uniform (cmesh, scheme, level, 1, sc_MPI_COMM_WORLD);
    /* exchange ghost data */
    t8_test_ghost_exchange_data_int (forest);
    t8_test_ghost_exchange_data_id (forest);
    /* Adapt the forest and exchange data again */
    int maxlevel = level + 2;
    t8_forest_t forest_adapt = t8_forest_new_adapt (forest, t8_test_exchange_adapt, 1, 1, &maxlevel);
    t8_test_ghost_exchange_data_int (forest_adapt);
    t8_test_ghost_exchange_data_id (forest_adapt);
    t8_forest_unref (&forest_adapt);
  }
}

INSTANTIATE_TEST_SUITE_P (t8_gtest_ghost_exchange, forest_ghost_exchange,
                          testing::Combine (AllSchemeCollections, AllCmeshsParam), pretty_print_base_example_scheme);
