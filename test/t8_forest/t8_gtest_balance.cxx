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

/** \file t8_gtest_balance.cxx
* Performs some tests of the balance functionality. 
*/

#include <gtest/gtest.h>
#include <test/t8_gtest_macros.hxx>
#include <test/t8_gtest_schemes.hxx>
#include <test/t8_gtest_custom_assertion.hxx>

#include <t8_eclass.h>
#include <t8_cmesh/t8_cmesh.h>
#include <t8_cmesh/t8_cmesh_examples.h>
#include <t8_forest/t8_forest_general.h>
#include <t8_geometry/t8_geometry_implementations/t8_geometry_linear.hxx>
#include <t8_schemes/t8_default/t8_default.hxx>
#include <t8_forest/t8_forest_balance.h>

#include <array>
#include <vector>
#include <algorithm>

struct gtest_balance: public testing::TestWithParam<std::tuple<std::tuple<int, t8_eclass_t>, int, int>>
{
 public:
  static const int kNumTrees = 4;

 protected:
  void
  SetUp () override
  {
    const int scheme_id = std::get<0> (std::get<0> (GetParam ()));
    scheme = create_from_scheme_id (scheme_id);
    eclass = std::get<1> (std::get<0> (GetParam ()));
    ilevel = std::get<1> (GetParam ());
    ido_periodic = std::get<2> (GetParam ());
  }
  t8_eclass_t eclass;
  const t8_scheme *scheme;
  int ilevel;
  int ido_periodic;
};

/**
 * \brief This test confirms that the function 't8_forest_is_balanced' recognizes uniform forests as balanced.
 */
TEST_P (gtest_balance, confirm_is_balanced_check_for_uniform_forests)
{
  if (eclass == t8_eclass_t::T8_ECLASS_PYRAMID && ido_periodic == 1) {
    scheme->unref ();
    GTEST_SKIP_ ("The pyramid cube mesh cannot be periodic.");
  }
  t8_cmesh_t cmesh = t8_cmesh_new_hypercube (eclass, sc_MPI_COMM_WORLD, 0, 0, ido_periodic);
  t8_forest_t forest = t8_forest_new_uniform (cmesh, scheme, ilevel, 0, sc_MPI_COMM_WORLD);

  EXPECT_EQ (t8_forest_is_balanced (forest), 1);

  t8_forest_unref (&forest);
}

struct gtest_balance_adapt_data
{
  std::vector<t8_gloidx_t> trees_to_refine;
  int max_refinement_level;
};

static int
t8_gtest_balance_refine_certain_trees (t8_forest_t forest, t8_forest_t forest_from, t8_locidx_t which_tree,
                                       const t8_eclass_t tree_class, [[maybe_unused]] t8_locidx_t lelement_id,
                                       const t8_scheme *scheme, [[maybe_unused]] const int is_family,
                                       [[maybe_unused]] const int num_elements, t8_element_t *elements[])
{
  gtest_balance_adapt_data *adapt_data = static_cast<gtest_balance_adapt_data *> (t8_forest_get_user_data (forest));

  const t8_gloidx_t gtree_id = t8_forest_global_tree_id (forest_from, which_tree);

  if (std::find (adapt_data->trees_to_refine.begin (), adapt_data->trees_to_refine.end (), gtree_id)
        != adapt_data->trees_to_refine.end ()
      && scheme->element_get_level (tree_class, elements[0]) < adapt_data->max_refinement_level) {
    return 1;
  }
  else {
    return 0;
  }
}

/**
 * \brief This function generates a simple forest which is not balanced.
 * \param [in] trees_to_refine A vector holding the IDs of the trees ought to be refined
 * \param [in] additional_refinement Indicates how much more refinement takes place within the supplied trees \var trees_to_refine
 * 
 * \note The resulting forest consists of four trees (\see t8_cmesh_new_hypercube_pad).
 * The elements of a tree are refined to the level '2 + \a additional_refinement'.
 * 
 * \note The forest being created looks similar to (in this case: \a additional_refinement = 0):
 *             Cmesh:             ->        Adapted Forest:
 *     __ __ __ __ __ __ __ __           __ __ __ __ __ __ __ __ 
 *    |           |           |         |           |           |
 *    |  Tree: 2  |  Tree: 3  |         |           |           |
 *    |           |           |         |           |           |
 *    |__ __ __ __|__ __ __ __|         |__ __ __ __|__ __ __ __|
 *    |           |           |         |__|__|__|__|           |
 *    |  Tree: 0  |  Tree: 1  |         |__|__|__|__|           |
 *    |           |           |         |__|__|__|__|           |
 *    |__ __ __ __|__ __ __ __|         |__|__|__|__|__ __ __ __|
 * 
 * \return The adapted forest; as shown above
 */
static t8_forest_t
t8_gtest_obtain_forest_for_balance_tests (const std::vector<t8_gloidx_t> &trees_to_refine, const t8_scheme *scheme,
                                          int additional_refinement = 0)
{
  const double boundary_coords[12] = { 0, 0, 0, 1, 0, 0, 0, 1, 0, 1, 1, 0 };

  t8_cmesh_t cmesh = t8_cmesh_new_hypercube_pad (T8_ECLASS_QUAD, sc_MPI_COMM_WORLD, boundary_coords, 2, 2, 1, 0);

  t8_forest_t forest;
  t8_forest_init (&forest);
  t8_forest_set_cmesh (forest, cmesh, sc_MPI_COMM_WORLD);
  t8_forest_set_scheme (forest, scheme);
  t8_forest_commit (forest);

  gtest_balance_adapt_data adapt_data;
  adapt_data.trees_to_refine = trees_to_refine;
  adapt_data.max_refinement_level = 2 + additional_refinement;

  return t8_forest_new_adapt (forest, t8_gtest_balance_refine_certain_trees, 1, 0, &adapt_data);
}

/**
 * \brief This function checks whether each tree only holds elements that are on the refinement level given by \a expected_elem_level_per_tree
 * 
 * \param [in] balanced_forest A forest consisting of gtest_balance::kNumTrees trees
 * \param [in] expected_elem_level_per_tree An array holding a refinement level for each tree id
 * \return true If each element with a tree corresponds to the given refinement level supplied by the array \var expected_elem_level_per_tree
 * \return false If not every element complies to the given refinement level per tree
 */
static bool
t8_gtest_check_custom_balanced_forest (t8_forest_t balanced_forest,
                                       const std::array<int, gtest_balance::kNumTrees> &expected_elem_level_per_tree)
{
  const t8_locidx_t num_local_trees = t8_forest_get_num_local_trees (balanced_forest);

  for (t8_locidx_t tree_id = 0; tree_id < num_local_trees; ++tree_id) {
    const t8_locidx_t num_tree_local_elems = t8_forest_get_tree_num_leaf_elements (balanced_forest, tree_id);

    const t8_gloidx_t gtree_id = t8_forest_global_tree_id (balanced_forest, tree_id);
    const t8_eclass_t tree_class = t8_forest_get_tree_class (balanced_forest, tree_id);
    const t8_scheme *scheme = t8_forest_get_scheme (balanced_forest);

    for (t8_locidx_t elem_id = 0; elem_id < num_tree_local_elems; ++elem_id) {
      const t8_element_t *element = t8_forest_get_leaf_element_in_tree (balanced_forest, tree_id, elem_id);

      const int elem_level = scheme->element_get_level (tree_class, element);

      if (elem_level != expected_elem_level_per_tree[gtree_id]) {
        return false;
      }
    }
  }

  return true;
}

/**
 * \brief The test receives a custom forest which is not balanced and checks whether the balanced version of the forest complies
 * with a maximum level difference of +/-1 between neighboring elements.
 * 
 * \note The resulting forest should look like:
 *         Adapted Forest:         ->       Balanced Forest:
 *     __ __ __ __ __ __ __ __           __ __ __ __ __ __ __ __ 
 *    |           |           |         |     |     |           |
 *    |           |           |         |__ __|__ __|           |
 *    |           |           |         |     |     |           |
 *    |__ __ __ __|__ __ __ __|         |__ __|__ __|__ __ __ __|
 *    |__|__|__|__|           |         |__|__|__|__|     |     |
 *    |__|__|__|__|           |         |__|__|__|__|__ __|__ __|
 *    |__|__|__|__|           |         |__|__|__|__|     |     |
 *    |__|__|__|__|__ __ __ __|         |__|__|__|__|__ __|__ __|
 * 
 */
TEST_P (gtest_balance, balance_adapted_forest_no_repartition)
{
  std::vector<t8_gloidx_t> trees_to_refine { 0 };
  t8_forest_t forest = t8_gtest_obtain_forest_for_balance_tests (trees_to_refine, scheme, 0);

  const int flag_no_repartition = 1;

  t8_forest_t balanced_forest;
  t8_forest_init (&balanced_forest);
  t8_forest_set_balance (balanced_forest, forest, flag_no_repartition);
  t8_forest_commit (balanced_forest);

  std::array<int, gtest_balance::kNumTrees> expected_elem_level_per_tree { 2, 1, 1, 0 };

  EXPECT_TRUE (t8_gtest_check_custom_balanced_forest (balanced_forest, expected_elem_level_per_tree));

  t8_forest_unref (&balanced_forest);
}

/**
 * \brief Tests whether an already balanced forest remains unchanged after another balance iteration.
 */
TEST_P (gtest_balance, balance_consistency_test)
{
  const int additional_refinement = 2;
  std::vector<t8_gloidx_t> trees_to_refine { 1, 2 };

  t8_forest_t forest = t8_gtest_obtain_forest_for_balance_tests (trees_to_refine, scheme, additional_refinement);

  int flag_no_repartition = 0;
  t8_forest_t balanced_forest;
  t8_forest_init (&balanced_forest);
  t8_forest_set_balance (balanced_forest, forest, flag_no_repartition);
  t8_forest_commit (balanced_forest);

  t8_forest_ref (balanced_forest);

  flag_no_repartition = 1;
  t8_forest_t already_balanced_forest;
  t8_forest_init (&already_balanced_forest);
  t8_forest_set_balance (already_balanced_forest, balanced_forest, flag_no_repartition);
  t8_forest_commit (already_balanced_forest);

  EXPECT_EQ (t8_forest_is_equal (balanced_forest, already_balanced_forest), 1);

  t8_forest_unref (&balanced_forest);
  t8_forest_unref (&already_balanced_forest);
}
#if T8_TEST_LEVEL_INT >= 2
const int maxlvl = 3;
#elif T8_TEST_LEVEL_INT >= 1
const int maxlvl = 4;
#else
const int maxlvl = 5;
#endif

INSTANTIATE_TEST_SUITE_P (t8_gtest_balance, gtest_balance,
                          testing::Combine (AllSchemes, testing::Range (0, maxlvl), testing::Range (0, 2)));
