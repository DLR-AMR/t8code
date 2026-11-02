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

/** \file t8_gtest_adapt_batched.cxx
 * Test the batched adapt callback \ref t8_forest_adapt_batched_t.
 * A batched callback and a per-element callback \ref t8_forest_adapt_t that decide with the same rule for each
 * element have to produce equal forests. The per-element callback coarsens a family if all its members vote for
 * coarsening, which is the rule of the batched adaptation.
 */

#include <gtest/gtest.h>
#include <t8_cmesh/t8_cmesh.h>
#include <t8_cmesh/t8_cmesh_examples.h>
#include <t8_forest/t8_forest_general.h>
#include <t8_forest/t8_forest_partition.h>
#include <t8_forest/t8_forest_private.h>
#include <t8_forest/t8_forest_types.h>
#include <t8_schemes/t8_default/t8_default.hxx>
#include <t8_schemes/t8_subelement/t8_subelement.hxx>
#include <t8_cmesh_generator/t8_cmesh_example_sets.hxx>
#include <t8_gtest_macros.hxx>
#include <t8_gtest_schemes.hxx>

#include <cstdint>

/** The actions the test callbacks may request, passed as user data of the adapted forest. */
struct batched_test_settings
{
  bool remove;      /**< True if the callbacks remove elements, false otherwise. */
  bool subelements; /**< True if the callbacks request subelements, false otherwise. */
};

/** Scramble the bits of \a value, so that consecutive values give unrelated results.
 * \param [in] value The value to scramble.
 * \return The scrambled value.
 */
static uint64_t
t8_test_mix (uint64_t value)
{
  value ^= value >> 33;
  value *= 0xff51afd7ed558ccdULL;
  value ^= value >> 33;
  value *= 0xc4ceb9fe1a85ec53ULL;
  value ^= value >> 33;
  return value;
}

/** The vote of a single element, the same for both callbacks.
 * All siblings vote for coarsening if their parent is selected, so that many families are coarsened. Single elements
 * vote for coarsening, too, which must not coarsen a family whose other members do not vote for it.
 * \param [in] forest      The forest under construction, its user data are the \ref batched_test_settings.
 * \param [in] forest_from The forest that is adapted.
 * \param [in] which_tree  The local tree of \a element.
 * \param [in] tree_class  The eclass of \a which_tree.
 * \param [in] scheme      The scheme of the forest.
 * \param [in] element     The element.
 * \return The adapt action of \a element, see \ref t8_forest_adapt_t.
 */
static int
t8_test_vote (const t8_forest_t forest, const t8_forest_t forest_from, const t8_locidx_t which_tree,
              const t8_eclass_t tree_class, const t8_scheme *scheme, const t8_element_t *element)
{
  const batched_test_settings *settings = static_cast<const batched_test_settings *> (t8_forest_get_user_data (forest));
  const int level = scheme->element_get_level (tree_class, element);
  const uint64_t tree_hash = t8_test_mix ((uint64_t) t8_forest_global_tree_id (forest_from, which_tree));
  /* The linear id at level - 1 is the id of the parent, which is the same for all siblings. */
  if (level > 0 && t8_test_mix (tree_hash + scheme->element_get_linear_id (tree_class, element, level - 1)) % 3 == 0) {
    return -1;
  }
  const uint64_t hash
    = t8_test_mix (tree_hash + scheme->element_get_linear_id (tree_class, element, level) * 64 + level);
  switch (hash % 8) {
  case 0:
    return 1;
  case 1:
    return -1;
  case 2:
    /* Never remove a first child, so that no tree loses all its elements. */
    return settings->remove && scheme->element_get_child_id (tree_class, element) != 0 ? -2 : 0;
  case 3:
    /* Subelement type + 1, the types 1 to 6 exist for quads and triangles. */
    return settings->subelements ? 2 + (int) ((hash / 8) % 6) : 0;
  default:
    return 0;
  }
}

/** Per-element adapt callback: coarsen a family if all its members vote for it, otherwise use the vote of the first
 * element, where a vote for coarsening means keep.
 * \see t8_forest_adapt_t */
static int
t8_test_adapt_element (t8_forest_t forest, t8_forest_t forest_from, t8_locidx_t which_tree,
                       const t8_eclass_t tree_class, [[maybe_unused]] t8_locidx_t lelement_id,
                       const t8_scheme_c *scheme, const int is_family, const int num_elements, t8_element_t *elements[])
{
  if (is_family) {
    bool all_coarsen = true;
    for (int ielement = 0; ielement < num_elements && all_coarsen; ielement++) {
      all_coarsen = t8_test_vote (forest, forest_from, which_tree, tree_class, scheme, elements[ielement]) == -1;
    }
    if (all_coarsen) {
      return -1;
    }
  }
  const int action = t8_test_vote (forest, forest_from, which_tree, tree_class, scheme, elements[0]);
  return action == -1 ? 0 : action;
}

/** Batched adapt callback: the vote of each element.
 * \see t8_forest_adapt_batched_t */
static void
t8_test_adapt_batched (t8_forest_t forest, t8_forest_t forest_from, t8_locidx_t which_tree,
                       const t8_eclass_t tree_class, const t8_scheme_c *scheme, const t8_element_array_t *elements,
                       int *actions)
{
  const size_t num_elements = t8_element_array_get_count (elements);
  for (size_t ielement = 0; ielement < num_elements; ielement++) {
    actions[ielement] = t8_test_vote (forest, forest_from, which_tree, tree_class, scheme,
                                      t8_element_array_index_locidx (elements, (t8_locidx_t) ielement));
  }
}

/** Adapt \a forest_from with one of the test callbacks. \a forest_from is not unreferenced.
 * \param [in] forest_from       The forest to adapt.
 * \param [in] batched           True to use \ref t8_test_adapt_batched, false to use \ref t8_test_adapt_element.
 * \param [in] settings          The actions the callbacks may request.
 * \param [in] balance_partition True to balance and partition in the same commit, false otherwise.
 * \return The adapted forest.
 */
static t8_forest_t
t8_test_adapt (t8_forest_t forest_from, const bool batched, batched_test_settings *settings,
               const bool balance_partition)
{
  t8_forest_t forest;
  t8_forest_init (&forest);
  t8_forest_ref (forest_from);
  t8_forest_set_user_data (forest, settings);
  if (batched) {
    t8_forest_set_adapt_batched (forest, forest_from, t8_test_adapt_batched);
  }
  else {
    t8_forest_set_adapt (forest, forest_from, t8_test_adapt_element, 0);
  }
  if (balance_partition) {
#if T8_TEST_LEVEL_INT < 2
    t8_forest_set_balance (forest, NULL, 0);
#endif
    t8_forest_set_partition (forest, NULL, 0);
  }
  t8_forest_commit (forest);
  return forest;
}

/** Adapt \a forest_from with both callbacks and check that the results are equal. \a forest_from is not
 * unreferenced.
 * \param [in] forest_from       The forest to adapt.
 * \param [in] settings          The actions the callbacks may request.
 * \param [in] balance_partition True to balance and partition in the same commit, false otherwise.
 */
static void
t8_test_expect_equal_adapt (t8_forest_t forest_from, batched_test_settings settings, const bool balance_partition)
{
  t8_forest_t forest_element = t8_test_adapt (forest_from, false, &settings, balance_partition);
  t8_forest_t forest_batched = t8_test_adapt (forest_from, true, &settings, balance_partition);
  EXPECT_TRUE (t8_forest_is_equal (forest_element, forest_batched));
  EXPECT_EQ (forest_element->incomplete_trees, forest_batched->incomplete_trees);
  t8_forest_unref (&forest_element);
  t8_forest_unref (&forest_batched);
}

struct forest_adapt_batched: public testing::TestWithParam<std::tuple<int, cmesh_example_base *>>
{
 protected:
  void
  SetUp () override
  {
    scheme = create_from_scheme_id (std::get<0> (GetParam ()));
    cmesh = std::get<1> (GetParam ())->cmesh_create ();
    if (t8_cmesh_is_empty (cmesh)) {
      /* A uniform forest needs a nonempty cmesh. */
      GTEST_SKIP ();
    }
    /* Start one level above the first nonempty level, so that families can be coarsened. */
    const int level = t8_forest_min_nonempty_level (cmesh, scheme) + 1;
    t8_cmesh_ref (cmesh);
    scheme->ref ();
    forest = t8_forest_new_uniform (cmesh, scheme, level, 0, sc_MPI_COMM_WORLD);
  }

  void
  TearDown () override
  {
    if (forest != nullptr) {
      t8_forest_unref (&forest);
    }
    scheme->unref ();
    t8_cmesh_destroy (&cmesh);
  }

  t8_cmesh_t cmesh;             /**< The coarse mesh of the test case. */
  const t8_scheme *scheme;      /**< The scheme of the test case. */
  t8_forest_t forest = nullptr; /**< The uniform forest to adapt, nullptr if the test is skipped. */
};

/* Refine, keep and coarsen. */
TEST_P (forest_adapt_batched, refine_coarsen)
{
  t8_test_expect_equal_adapt (forest, { false, false }, false);
}

/* Additionally remove elements, which makes the forest incomplete. */
TEST_P (forest_adapt_batched, remove)
{
  t8_test_expect_equal_adapt (forest, { true, false }, false);
}

/* Adapt a forest with holes. */
TEST_P (forest_adapt_batched, holes)
{
  batched_test_settings remove_settings = { true, false };
  t8_forest_t forest_holes = t8_test_adapt (forest, false, &remove_settings, false);
  t8_test_expect_equal_adapt (forest_holes, { false, false }, false);
  t8_forest_unref (&forest_holes);
}

/* Adapt, balance and partition in one commit, which adapts an intermediate forest. */
TEST_P (forest_adapt_batched, balance_partition)
{
  t8_test_expect_equal_adapt (forest, { false, false }, true);
}

INSTANTIATE_TEST_SUITE_P (t8_gtest_adapt_batched, forest_adapt_batched,
                          testing::Combine (AllSchemeCollections, AllCmeshsParam), pretty_print_base_example_scheme);

/* Request subelements with a subelement scheme on a hybrid mesh of quads and triangles. */
TEST (t8_gtest_adapt_batched, subelements)
{
  t8_cmesh_t cmesh;
  t8_cmesh_init (&cmesh);
  t8_cmesh_new_hypercube_hybrid_2D (cmesh, sc_MPI_COMM_WORLD);
  t8_forest_t forest = t8_forest_new_uniform (cmesh, t8_scheme_new_subelement (), 2, 0, sc_MPI_COMM_WORLD);
  t8_test_expect_equal_adapt (forest, { false, true }, false);
  t8_forest_unref (&forest);
}
