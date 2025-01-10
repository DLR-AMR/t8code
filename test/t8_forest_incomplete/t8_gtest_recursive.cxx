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
#include <t8_forest/t8_forest.h>
#include <test/t8_gtest_schemes.hxx>
#include <test/t8_gtest_macros.hxx>

/* In this test, we recursively constructs a mesh containing only the first 
 * and last elements of a family. Furthermore, for every two elements 
 * e_a and e_b that do not belong to the same family, it holds that 
 * if linear_id(e_a) < linear_id(a_b) then level(e_a) > level(e_b).
 * Thus, recursive coarsening should result in a tree containing only the
 * root element.
 * 
 * Note, that each rank has its own local/global tree. No trees are shared.
 */

class recursive_tree: public testing::TestWithParam<int> {
 protected:
  void
  SetUp () override
  {
    scheme_id = GetParam ();
    scheme = t8_scheme_all_schemes ();
    tree_class = scheme->get_eclass_scheme_eclass (static_cast<t8_eclass_t>(scheme_id));
    if (tree_class == T8_ECLASS_ZERO) {
      GTEST_SKIP ();
    }
    sc_MPI_Comm_size (sc_MPI_COMM_WORLD, &MPI_size);

    /* Construct a cmesh such that each process will get one rooted tree */
    cmesh = t8_cmesh_new_bigmesh (tree_class, MPI_size, sc_MPI_COMM_WORLD);

    scheme->ref ();
    t8_cmesh_ref (cmesh);

    /* The forest to be adapted. */
    forest = t8_forest_new_uniform (cmesh, scheme, 1, 0, sc_MPI_COMM_WORLD);
    /* The forest contains only root elements and serves as a comparison. */
    forest_base = t8_forest_new_uniform (cmesh, scheme, 0, 0, sc_MPI_COMM_WORLD);
  }
  void
  TearDown () override
  {
    if (tree_class != T8_ECLASS_ZERO) {
      t8_forest_unref (&forest);
      t8_forest_unref (&forest_base);
    }
    else {
      scheme->unref ();
    }
  }
  int MPI_size;
  int scheme_id;
  t8_eclass_t tree_class;
  t8_scheme *scheme;
  t8_cmesh_t cmesh;
  t8_forest_t forest;
  t8_forest_t forest_base;
};

/** Remove every element except last and first of a family. */
static int
t8_adapt_remove_but_last_first (t8_forest_t forest, t8_forest_t forest_from, t8_locidx_t which_tree,
                                const t8_eclass_t tree_class, t8_locidx_t lelement_id, const t8_scheme *scheme,
                                const int is_family, const int num_elements, t8_element_t *elements[])
{
  const int num_children = scheme->element_get_num_children (tree_class, elements[0]);
  const int child_id = scheme->element_get_child_id (tree_class, elements[0]);
  if (num_children - 1 != child_id && 0 != child_id) {
    return -2;
  }
  return 0;
}

/** Refine the first element of a family. */
static int
t8_adapt_refine_first (t8_forest_t forest, t8_forest_t forest_from, t8_locidx_t which_tree,
                       const t8_eclass_t tree_class, t8_locidx_t lelement_id, const t8_scheme *scheme,
                       const int is_family, const int num_elements, t8_element_t *elements[])
{
  const int level = scheme->element_get_level (tree_class, elements[0]);
  const int level_max = scheme->get_maxlevel (tree_class);
  const int child_id = scheme->element_get_child_id (tree_class, elements[0]);
  if (child_id == 0 && level < (int) (0.2 * level_max)) {
    return 1;
  }
  return 0;
}

/** Refine every element. */
static int
t8_adapt_refine_all (t8_forest_t forest, t8_forest_t forest_from, t8_locidx_t which_tree, const t8_eclass_t tree_class,
                     t8_locidx_t lelement_id, const t8_scheme *scheme, const int is_family, const int num_elements,
                     t8_element_t *elements[])
{
  return 1;
}

/** Coarse every family. */
static int
t8_adapt_coarse_all (t8_forest_t forest, t8_forest_t forest_from, t8_locidx_t which_tree, const t8_eclass_t tree_class,
                     t8_locidx_t lelement_id, const t8_scheme *scheme, const int is_family, const int num_elements,
                     t8_element_t *elements[])
{
  if (is_family) {
    return -1;
  }
  return 0;
}

static t8_forest_t
t8_adapt_forest (t8_forest_t forest_from, t8_forest_adapt_t adapt_fn, int recursive)
{
  t8_forest_t forest_new;

  t8_forest_init (&forest_new);
  t8_forest_set_adapt (forest_new, forest_from, adapt_fn, recursive);
  t8_forest_commit (forest_new);

  return forest_new;
}

TEST_P (recursive_tree, test_recursive)
{
  forest = t8_adapt_forest (forest, t8_adapt_refine_first, 1);
  forest = t8_adapt_forest (forest, t8_adapt_remove_but_last_first, 0);
  forest = t8_adapt_forest (forest, t8_adapt_refine_all, 0);
  forest = t8_adapt_forest (forest, t8_adapt_remove_but_last_first, 1);
  forest = t8_adapt_forest (forest, t8_adapt_coarse_all, 1);

  /* The adaptet forest should only contain root elements as forest_base */
  ASSERT_TRUE (t8_forest_is_equal (forest, forest_base));
}

INSTANTIATE_TEST_SUITE_P (t8_gtest_recursive, recursive_tree, AllSchemes);
