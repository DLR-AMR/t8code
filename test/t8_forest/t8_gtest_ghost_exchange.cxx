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
#include <t8_eclass/t8_eclass.h>
#include <t8_schemes/t8_default/t8_default.hxx>
#include <t8_forest/t8_forest_general.h>
#include <t8_forest/t8_forest_ghost/t8_forest_ghost.h>
#include <t8_forest/t8_forest_ghost/t8_forest_ghost.hxx>
#include <t8_data/t8_vector_handler.hxx>
#include <t8_data/t8_data_handler_specs.hxx>
#include <t8_forest/t8_forest_private.h>
#include <t8_cmesh/t8_cmesh.h>
#include <t8_cmesh_generator/t8_cmesh_example_sets.hxx>
#include <t8_gtest_macros.hxx>
#include <t8_gtest_schemes.hxx>
#include <utility>
#include <vector>

/* This test program tests the forest ghost exchange routine.
 * Given a forest for which the ghost layer was created and an array
 * storing data for the local elements and the ghost elements, ghost_exchange
 * communicates the data of the local elements to the ghost entries of the
 * processes for which these elements are ghost.
 * We test the ghost exchange routine for several forests on different
 * coarse meshes.
 * One test is an integer entry '42' for each element,
 * in a second test, we store the element's linear id in the data array.
 * A third test exchanges the linear id and the level of the elements through the C++ interface with vector
 * handlers, as bytewise items and as items that a passed data handler serializes.
 */

struct forest_ghost_exchange: public testing::TestWithParam<std::tuple<int, cmesh_example_base *>>
{
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

/** Compute the linear id and the level of every local element of a forest, followed by those of its ghosts.
 * \param [in] forest  A committed forest.
 * \return             For each local element and each ghost, its linear id and its level.
 */
static std::vector<std::pair<t8_linearidx_t, int>>
t8_test_ghost_exchange_ids_levels (t8_forest_t forest)
{
  const t8_scheme *scheme = t8_forest_get_scheme (forest);
  std::vector<std::pair<t8_linearidx_t, int>> ids_levels;
  for (t8_locidx_t itree = 0; itree < t8_forest_get_num_local_trees (forest); itree++) {
    const t8_eclass_t tree_class = t8_forest_get_tree_class (forest, itree);
    for (t8_locidx_t ielem = 0; ielem < t8_forest_get_tree_num_leaf_elements (forest, itree); ielem++) {
      const t8_element_t *elem = t8_forest_get_leaf_element_in_tree (forest, itree, ielem);
      const int level = scheme->element_get_level (tree_class, elem);
      ids_levels.emplace_back (scheme->element_get_linear_id (tree_class, elem, level), level);
    }
  }
  for (t8_locidx_t itree = 0; itree < t8_forest_get_num_ghost_trees (forest); itree++) {
    const t8_eclass_t tree_class = t8_forest_ghost_get_tree_class (forest, itree);
    for (t8_locidx_t ielem = 0; ielem < t8_forest_ghost_tree_num_leaf_elements (forest, itree); ielem++) {
      const t8_element_t *elem = t8_forest_ghost_get_leaf_element (forest, itree, ielem);
      const int level = scheme->element_get_level (tree_class, elem);
      ids_levels.emplace_back (scheme->element_get_linear_id (tree_class, elem, level), level);
    }
  }
  return ids_levels;
}

/** The element data of the vector handler tests: an element's linear id and level. */
struct t8_test_ghost_id_level
{
  t8_linearidx_t id; /**< The linear id of the element. */
  int level;         /**< The level of the element. */
};

/* Fill the entries of the local elements with their linear id and level, exchange the ghost data with the C++
 * interface and check that every ghost entry holds the linear id and level of its ghost. This is done with bytewise
 * items and with cells that a passed data handler serializes. Each handler is exchanged twice, the second time
 * starting with the ghost entries of the first exchange.
 */
static void
t8_test_ghost_exchange_data_handler (t8_forest_t forest)
{
  const std::vector<std::pair<t8_linearidx_t, int>> ids_levels = t8_test_ghost_exchange_ids_levels (forest);
  const size_t num_elements = t8_forest_get_local_num_leaf_elements (forest);
  ASSERT_EQ (ids_levels.size (), num_elements + t8_forest_get_num_ghosts (forest));

  std::vector<t8_test_ghost_id_level> items;
  std::vector<t8_test_cell> cells;
  for (size_t ielem = 0; ielem < num_elements; ++ielem) {
    const auto [id, level] = ids_levels[ielem];
    items.push_back ({ id, level });
    cells.push_back ({ { static_cast<double> (id), static_cast<double> (level) } });
  }
  t8_vector_handler item_handler (std::move (items));
  t8_vector_handler cell_handler (std::move (cells), t8_test_cell_handler { 2 });

  for (int iexchange = 0; iexchange < 2; ++iexchange) {
    t8_forest_ghost_exchange_data (forest, item_handler);
    t8_forest_ghost_exchange_data (forest, cell_handler);
    ASSERT_EQ (item_handler.get_data ().size (), ids_levels.size ());
    ASSERT_EQ (cell_handler.get_data ().size (), ids_levels.size ());
    for (size_t ientry = 0; ientry < ids_levels.size (); ++ientry) {
      const auto [id, level] = ids_levels[ientry];
      const t8_test_ghost_id_level &item = item_handler.get_data ()[ientry];
      ASSERT_TRUE (item.id == id && item.level == level) << "Wrong bytewise entry " << ientry;
      const std::vector<double> expected_values = { static_cast<double> (id), static_cast<double> (level) };
      ASSERT_EQ (cell_handler.get_data ()[ientry].values, expected_values) << "Wrong serialized entry " << ientry;
    }
  }
}

/* Exchange the ghost data of a forest without ghost layer with the C++ interface: The handler only keeps the entries
 * of the local elements and drops the trailing ones.
 */
static void
t8_test_ghost_exchange_data_handler_no_ghosts (t8_forest_t forest)
{
  ASSERT_EQ (t8_forest_get_num_ghosts (forest), 0);
  const size_t num_elements = t8_forest_get_local_num_leaf_elements (forest);
  std::vector<int> data (num_elements + 3, 42);
  t8_vector_handler handler (std::move (data));
  t8_forest_ghost_exchange_data (forest, handler);
  EXPECT_EQ (handler.get_data (), std::vector<int> (num_elements, 42));
}

TEST_P (forest_ghost_exchange, test_ghost_exchange)
{

  /* Compute the minimum level, such that the forest is nonempty */
  int min_level = t8_forest_min_nonempty_level (cmesh, scheme);
  /* we start with an empty level */
  min_level = SC_MAX (min_level - 1, 0);
#if T8_TEST_LEVEL_INT >= 2
  const int max_level = min_level + 1;
#else
  const int max_level = min_level + 2;
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
    t8_test_ghost_exchange_data_handler (forest);
    /* Adapt the forest and exchange data again */
    int adapt_maxlevel = level + 2;
    t8_forest_t forest_adapt = t8_forest_new_adapt (forest, t8_test_exchange_adapt, 1, 1, &adapt_maxlevel);
    t8_test_ghost_exchange_data_int (forest_adapt);
    t8_test_ghost_exchange_data_id (forest_adapt);
    t8_test_ghost_exchange_data_handler (forest_adapt);
    t8_forest_unref (&forest_adapt);

    /* A forest without ghost layer */
    scheme->ref ();
    t8_cmesh_ref (cmesh);
    t8_forest_t forest_no_ghosts = t8_forest_new_uniform (cmesh, scheme, level, 0, sc_MPI_COMM_WORLD);
    t8_test_ghost_exchange_data_handler_no_ghosts (forest_no_ghosts);
    t8_forest_unref (&forest_no_ghosts);
  }
}

INSTANTIATE_TEST_SUITE_P (t8_gtest_ghost_exchange, forest_ghost_exchange,
                          testing::Combine (AllSchemeCollections, AllCmeshsParam), pretty_print_base_example_scheme);
