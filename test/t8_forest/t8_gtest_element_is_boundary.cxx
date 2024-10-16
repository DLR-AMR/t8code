/*
  This file is part of t8code.
  t8code is a C library to manage a collection (a forest) of multiple
  connected adaptive space-trees of general element classes in parallel.

  Copyright (C) 2024 the developers

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
#include <t8_cmesh.h>
#include <t8_forest/t8_forest_general.h>
#include <t8_schemes/t8_default/t8_default.hxx>
#include <t8_forest/t8_forest_private.h>
#include "test/t8_cmesh_generator/t8_cmesh_example_sets.hxx"
#include <test/t8_gtest_macros.hxx>

/* In this test we check the t8_forest_element_is_boundary function.
 * Iterating over all cmesh test cases, we create a uniform and an adaptive forest.
 * For each forest, we check that for each element t8_forest_element_is_boundary returns true
 * if and only if the element is at the domain boundary.
 */

/* Maximum uniform level for forest. */
#define T8_IS_BOUNDARY_MAX_LVL 3

/* Adapt a forest such that always the first child of a
 * family is refined and no other elements. This results in a highly
 * imbalanced forest. */
static int
t8_test_adapt_first_child (t8_forest_t forest, t8_forest_t forest_from, t8_locidx_t which_tree, t8_locidx_t lelement_id,
                           t8_eclass_scheme_c *ts, const int is_family, const int num_elements,
                           t8_element_t *elements[])
{
  int level = ts->t8_element_level (elements[0]);

  /* we set a maximum refinement level as forest user data */
  int maxlevel = *(int *) t8_forest_get_user_data (forest);
  if (level >= maxlevel) {
    /* Do not refine after the maxlevel */
    return 0;
  }
  int child_id = ts->t8_element_child_id (elements[0]);
  if (child_id == 1) {
    return 1;
  }
  return 0;
}

/* In a quad forest remove the first and fourth child of each element.
 * This will reside in a forest where each element and each face is a boundary element.
 * */
static int
t8_test_adapt_quad_remove_first_and_fourth_child (t8_forest_t forest, t8_forest_t forest_from, t8_locidx_t which_tree,
                                                  t8_locidx_t lelement_id, t8_eclass_scheme_c *ts, const int is_family,
                                                  const int num_elements, t8_element_t *elements[])
{
  /* For the purpose of this test, this function should only get called with quad elements. */
  SC_CHECK_ABORT (ts->t8_element_shape (elements[0]) == T8_ECLASS_QUAD,
                  "Special test adapt function must only be used on quad elements.\n");
  int child_id = ts->t8_element_child_id (elements[0]);
  /* Remove child_id 0 and 3, do not change any other element. */
  if (child_id == 0 || child_id == 3) {
    return -2;
  }
  return 0;
}

class element_is_boundary: public testing::TestWithParam<std::tuple<int, cmesh_example_base *>> {
 protected:
  void
  SetUp () override
  {
    /* Construct a cmesh */
    const int level = std::get<0> (GetParam ());
    cmesh = std::get<1> (GetParam ())->cmesh_create ();
    if (t8_cmesh_is_empty (cmesh)) {
      /* forest_commit does not support empty cmeshes, we skip this case */
      GTEST_SKIP ();
    }
    /* Build the default scheme (TODO: Test this with all schemes) */
    scheme = t8_scheme_new_default_cxx ();
    forest = t8_forest_new_uniform (cmesh, scheme, level, 0, sc_MPI_COMM_WORLD);
    t8_forest_ref (forest);
    int maxlevel = 7;
    const int recursive_adapt = 1;
    forest_adapt = t8_forest_new_adapt (forest, t8_test_adapt_first_child, recursive_adapt, 0, &maxlevel);
  }

  void
  TearDown () override
  {
    if (t8_cmesh_is_empty (cmesh)) {
      t8_cmesh_destroy (&cmesh);
    }
    else {
      t8_forest_unref (&forest);
      t8_forest_unref (&forest_adapt);
    }
  }

  t8_cmesh_t cmesh;
  t8_forest_t forest, forest_adapt;
  t8_scheme_cxx_t *scheme;
};

void
t8_test_element_is_boundary_for_forest (t8_forest_t forest, t8_cmesh_t cmesh,
                                        const int each_element_face_is_expected_boundary)
{
  const t8_locidx_t num_local_trees = t8_forest_get_num_local_trees (forest);

  for (t8_locidx_t itree = 0; itree < num_local_trees; ++itree) {
    const t8_locidx_t num_elements_in_tree = t8_forest_get_tree_num_elements (forest, itree);
    const t8_eclass_t tree_class = t8_forest_get_tree_class (forest, itree);
    const t8_eclass_scheme_c *scheme = t8_forest_get_eclass_scheme (forest, tree_class);
    /* Iterate over all the tree's leaf elements, check whether the leaf
     * is correctly identified by t8_forest_element_is_boundary,
     * build its parent and its first child (if they exist), and verify
     * that t8_forest_element_is_boundary returns false. */
    for (t8_locidx_t ielement = 0; ielement < num_elements_in_tree; ++ielement) {
      const t8_element_t *leaf_element = t8_forest_get_element_in_tree (forest, itree, ielement);
      const int num_element_faces = scheme->t8_element_num_faces (leaf_element);
      for (int iface = 0; iface < num_element_faces; ++iface) {
        /* Iterate over all faces */
        int face_is_at_boundary = 0;
        if (!each_element_face_is_expected_boundary) {
          /* Manually check whether this element is at the boundary */
          if (scheme->t8_element_is_root_boundary (leaf_element, iface)) {
            /* This face is at a tree boundary, so it might be at the domain boundary. */
            /* We check whether the tree is at the domain boundary. */
            const int tree_face = scheme->t8_element_tree_face (leaf_element, iface);
            const t8_locidx_t cmesh_local_tree = t8_forest_ltreeid_to_cmesh_ltreeid (forest, itree);
            if (t8_cmesh_tree_face_is_boundary (cmesh, cmesh_local_tree, tree_face)) {
              /* This element is at the domain boundary. */
              face_is_at_boundary = 1;
            }
          }
        }
        else {
          /* each element and face is an expected boundary */
          face_is_at_boundary = 1;
        }
        EXPECT_EQ (t8_forest_leaf_is_boundary (forest, itree, leaf_element, iface), face_is_at_boundary);
      }
    }
  }
}

TEST_P (element_is_boundary, element_is_boundary)
{
  t8_test_element_is_boundary_for_forest (forest, cmesh, 0);
}

TEST_P (element_is_boundary, element_is_boundary_adapt)
{
  t8_test_element_is_boundary_for_forest (forest_adapt, cmesh, 0);
}

TEST (element_is_boundary, quad_forest_with_holes)
{
  /* This test is deactivated as long as Ghost does not work in combination
   * with deleted elements.
   * Once this is fixed, reactivate this test.
   * See https://github.com/DLR-AMR/t8code/issues/825
   */
  GTEST_SKIP ();
  /* Create a 10 x 5 2D brick cmesh, periodic in x direction. */
  t8_cmesh_t cmesh = t8_cmesh_new_brick_2d (10, 5, 1, 0, sc_MPI_COMM_WORLD);

  t8_scheme_cxx_t *scheme = t8_scheme_new_default_cxx ();
  t8_forest_t forest = t8_forest_new_uniform (cmesh, scheme, T8_IS_BOUNDARY_MAX_LVL, 0, sc_MPI_COMM_WORLD);
  t8_forest_t forest_adapt = t8_forest_new_adapt (forest, t8_test_adapt_quad_remove_first_and_fourth_child, 0, 1, NULL);

  t8_forest_write_vtk (forest_adapt, "test_quad_w_holes");

  t8_test_element_is_boundary_for_forest (forest_adapt, cmesh, 1);
  t8_forest_unref (&forest_adapt);
}

/* This test class creates a single tree cmesh for each eclass and
 * builds a uniform level 0 forest on it.
 * For the single element in that forest all faces lie on the boundary. */
class element_is_boundary_known_boundary: public testing::TestWithParam<t8_eclass_t> {
 protected:
  void
  SetUp () override
  {
    eclass = GetParam ();
    cmesh = t8_cmesh_new_from_class (eclass, sc_MPI_COMM_WORLD);
    t8_scheme_cxx_t *scheme = t8_scheme_new_default_cxx ();
    forest = t8_forest_new_uniform (cmesh, scheme, 0, 0, sc_MPI_COMM_WORLD);
  }

  void
  TearDown () override
  {
    t8_forest_unref (&forest);
  }

  t8_cmesh_t cmesh;
  t8_forest_t forest;
  t8_eclass_t eclass;
};

/* For a level 0 forest of a single class cmesh all faces should
 * be at the boundary. */
TEST_P (element_is_boundary_known_boundary, level_0)
{
  const t8_locidx_t num_elements = t8_forest_get_local_num_elements (forest);
  if (num_elements > 0) {
    T8_ASSERT (num_elements == 1);
    const t8_element_t *element = t8_forest_get_element_in_tree (forest, 0, 0);
    const t8_eclass_scheme_c *scheme = t8_forest_get_eclass_scheme (forest, eclass);
    const int num_faces = scheme->t8_element_num_faces (element);
    for (int iface = 0; iface < num_faces; ++iface) {
      EXPECT_TRUE (t8_forest_leaf_is_boundary (forest, 0, element, iface));
    }
  }
}

/* Define a lambda to beautify gtest output for tuples <level, cmesh>.
 * This will set the correct level and cmesh name as part of the test case name. */
auto pretty_print_level_and_cmesh_params
  = [] (const testing::TestParamInfo<std::tuple<int, cmesh_example_base *>> &info) {
      std::string name = std::string ("Level_") + std::to_string (std::get<0> (info.param));
      std::string cmesh_name;
      std::get<1> (info.param)->param_to_string (cmesh_name);
      name += std::string ("_") + cmesh_name;
      return name;
    };

INSTANTIATE_TEST_SUITE_P (t8_gtest_element_is_boundary, element_is_boundary,
                          testing::Combine (testing::Range (0, T8_IS_BOUNDARY_MAX_LVL), AllCmeshsParam),
                          pretty_print_level_and_cmesh_params);

INSTANTIATE_TEST_SUITE_P (t8_gtest_element_is_boundary, element_is_boundary_known_boundary, AllEclasses, print_eclass);
