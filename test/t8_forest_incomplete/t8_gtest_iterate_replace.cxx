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
#include <t8.h>
#include <t8_cmesh/t8_cmesh_examples.h>
#include "test/t8_cmesh_generator/t8_cmesh_example_sets.hxx"
#include <t8_forest/t8_forest.h>
#include <t8_forest/t8_forest_iterate.h>
#include <t8_schemes/t8_default/t8_default_cxx.hxx>
#include <test/t8_gtest_macros.hxx>

/* In this test, we first adapt a forest and store every callback return value.
 * In the next step, we call t8_forest_iterate_replace. Instead of interpolating
 * the stored data, we check in every callback call inside of
 * t8_forest_iterate_replace if it is passed the correct values.
 */

class forest_iterate: public testing::TestWithParam<cmesh_example_base *> {
 protected:
  void
  SetUp () override
  {
    cmesh = GetParam ()->cmesh_create ();
    if (t8_cmesh_is_empty (cmesh)) {
      /* empty cmeshes are currently not supported */
      GTEST_SKIP ();
    }
    forest = t8_forest_new_uniform (cmesh, t8_scheme_new_default_cxx (), 4, 0, sc_MPI_COMM_WORLD);
  }
  void
  TearDown () override
  {
    if (t8_cmesh_is_empty (cmesh)) {
      t8_cmesh_unref (&cmesh);
    }
    else {
      t8_forest_unref (&forest);
    }
  }
  t8_cmesh_t cmesh;
  t8_forest_t forest;
};

/** This structure contains an array with all return values of all
 * callback function calls in the adaptation process of a forest.
 */
struct t8_return_data
{
  int *callbacks;
};

/** Inside the callback of iterate_replace we compare \a refine 
 * with the according return value of the callback of forest_adapt.
 * If true, we check the parameter \a num_outgoing, \a first_outgoing
 * \a num_incoming and \a first_incoming for correctness. */
void
t8_forest_replace (t8_forest_t forest_old, t8_forest_t forest_new, t8_locidx_t which_tree, t8_eclass_scheme_c *ts,
                   int refine, int num_outgoing, t8_locidx_t first_outgoing, int num_incoming,
                   t8_locidx_t first_incoming)
{
  /* Note, the new forest contains the callback returns of the old forest */
  struct t8_return_data *adapt_data = (struct t8_return_data *) t8_forest_get_user_data (forest_new);
  T8_ASSERT (adapt_data != NULL);

  /* Local element index of the old and new forest. */
  t8_locidx_t elidx_old = first_outgoing;
  t8_locidx_t elidx_new = first_incoming;
  for (t8_locidx_t tidx = 0; tidx < which_tree; tidx++) {
    elidx_new += t8_forest_get_tree_num_elements (forest_new, tidx);
    elidx_old += t8_forest_get_tree_num_elements (forest_old, tidx);
  }

  ASSERT_EQ (adapt_data->callbacks[elidx_old], refine);

  /* Element remained untouched. */
  if (refine == 0) {
    ASSERT_EQ (num_outgoing, 1);
    ASSERT_EQ (num_incoming, 1);
  }
  /* Element/family got coarsened. */
  if (refine == -1) {
    ASSERT_EQ (num_incoming, 1);

    /* Begin check family */
    const t8_element_t *parent = t8_forest_get_element_in_tree (forest_new, which_tree, first_incoming);
    t8_element_t *parent_compare;
    ts->t8_element_new (1, &parent_compare);
    int family_size = 1;
    t8_locidx_t tree_num_elements_old = t8_forest_get_tree_num_elements (forest_old, which_tree);
    for (t8_locidx_t elidx = 1;
         elidx < ts->t8_element_num_children (parent) && elidx + first_outgoing < tree_num_elements_old; elidx++) {
      const t8_element_t *child = t8_forest_get_element_in_tree (forest_old, which_tree, first_outgoing + elidx);
      ts->t8_element_parent (child, parent_compare);
      if (ts->t8_element_equal (parent, parent_compare)) {
        family_size++;
      }
    }
    ts->t8_element_destroy (1, &parent_compare);
    ASSERT_EQ (num_outgoing, family_size);
    /* End check family */

    /* If element got coarsen, only the first element
     * should be called in the callback of forest_adapt. */
    for (t8_locidx_t i = 1; i < num_outgoing; i++) {
      ASSERT_EQ (adapt_data->callbacks[elidx_old + i], -3);
    }
  }
  /* Element got removed. */
  if (refine == -2) {
    ASSERT_EQ (num_outgoing, 1);
    ASSERT_EQ (num_incoming, 0);
    ASSERT_EQ (first_incoming, -1);
  }
  /* Element got refined. */
  if (refine == 1) {
    ASSERT_EQ (num_outgoing, 1);
    const t8_element_t *element = t8_forest_get_element_in_tree (forest_old, which_tree, first_outgoing);
    const t8_locidx_t family_size = ts->t8_element_num_children (element);
    ASSERT_EQ (num_incoming, family_size);
  }
}

/** For each locale element: Remove, coarsen, leave untouched, or refine it depending on its index.
 *      if \a lelement_id mod 12 < 3  -> remove element
 * else if \a lelement_id mod 12 < 6  -> coarse element
 * else if \a lelement_id mod 12 < 9  -> leave element untouched
 * else if \a lelement_id mod 12 < 12 -> refine element
*/
int
t8_adapt_callback (t8_forest_t forest, t8_forest_t forest_from, t8_locidx_t which_tree, t8_locidx_t lelement_id,
                   t8_eclass_scheme_c *ts, const int is_family, const int num_elements, t8_element_t *elements[])
{
  struct t8_return_data *return_data = (struct t8_return_data *) t8_forest_get_user_data (forest);
  T8_ASSERT (return_data != NULL);

  int return_val = lelement_id % 12;
  if (return_val < 3) {
    return_val = -2;
    if (t8_forest_get_eclass (forest_from, which_tree) == T8_ECLASS_VERTEX) {
      return_val = 0;
    }
  }
  else if (return_val < 6) {
    if (is_family) {
      return_val = -1;
    }
    else {
      return_val = 0;
    }
  }
  else if (return_val < 9) {
    return_val = 0;
  }
  else {
    return_val = 1;
  }
  T8_ASSERT (-3 < return_val);
  T8_ASSERT (return_val < 2);

  /* Get the local index of current element in the local forest. */
  t8_locidx_t lelement_id_forest = lelement_id;
  for (t8_locidx_t tidx = 0; tidx < which_tree; tidx++) {
    lelement_id_forest += t8_forest_get_tree_num_elements (forest_from, tidx);
  }
  /* Store the return value. */
  return_data->callbacks[lelement_id_forest] = return_val;
  return return_val;
}

t8_forest_t
t8_adapt_forest (t8_forest_t forest_from, t8_forest_adapt_t adapt_fn, int do_adapt, int do_partition, void *user_data)
{
  t8_forest_t forest_new;

  t8_forest_init (&forest_new);
  if (do_adapt) {
    t8_forest_set_adapt (forest_new, forest_from, adapt_fn, 0);
    if (do_partition) {
      t8_forest_set_partition (forest_new, NULL, 0);
    }
  }
  else if (do_partition) {
    t8_forest_set_partition (forest_new, forest_from, 0);
  }
  if (user_data != NULL) {
    t8_forest_set_user_data (forest_new, user_data);
  }
  t8_forest_commit (forest_new);

  return forest_new;
}

TEST_P (forest_iterate, test_iterate_replace)
{
#if T8_ENABLE_LESS_TESTS
  const int runs = 1;
#else
  const int runs = 2;
#endif
  for (int run = 0; run < runs; run++) {
    t8_locidx_t num_elements = t8_forest_get_local_num_elements (forest);
    int *adapt_callbacks = T8_ALLOC (int, num_elements);

    for (t8_locidx_t elidx = 0; elidx < num_elements; elidx++) {
      adapt_callbacks[elidx] = -3;
    }
    struct t8_return_data data
    {
      adapt_callbacks
    };

    t8_forest_ref (forest);
    t8_forest_t forest_adapt;
    forest_adapt = t8_adapt_forest (forest, t8_adapt_callback, 1, 0, &data);

    t8_forest_iterate_replace (forest_adapt, forest, t8_forest_replace);
    t8_forest_unref (&forest);

    /* Partition the forest. This is useful,
     * if we run the test a second time with the adapted forest. */
    forest_adapt = t8_adapt_forest (forest_adapt, NULL, 0, 1, NULL);

    T8_FREE (adapt_callbacks);
    forest = forest_adapt;
  }
}

INSTANTIATE_TEST_SUITE_P (t8_gtest_iterate_replace, forest_iterate, AllCmeshsParam, pretty_print_base_example);
