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
#include <t8_forest/t8_forest_ghost.h>
#include <t8_forest/t8_forest_private.h>
#include <t8_cmesh.h>
#include "test/t8_cmesh_generator/t8_cmesh_example_sets.hxx"
#include <test/t8_gtest_macros.hxx>

/* This test program tests the forest ghost layer.
 * We adapt a forest and create its ghost layer. Afterwards, we
 * parse through all ghost elements and test whether the owner of an
 * element is in face the owner that is stored in the ghost layer.
  */

class forest_ghost_owner: public testing::TestWithParam<cmesh_example_base *> {
 protected:
  void
  SetUp () override
  {

    scheme = t8_scheme_new_default_cxx ();
    /* Construct a cmesh */
    cmesh = GetParam ()->cmesh_create ();
    if (t8_cmesh_is_empty (cmesh)) {
      /* empty cmeshes are currently not supported */
      GTEST_SKIP ();
    }
  }
  void
  TearDown () override
  {
    t8_cmesh_destroy (&cmesh);
    t8_schemexx_unref (&scheme);
  }
  t8_cmesh_t cmesh;
  t8_scheme *scheme;
};

static int
t8_test_gao_adapt (t8_forest_t forest, t8_forest_t forest_from, t8_locidx_t which_tree, t8_locidx_t lelement_id,
                   t8_scheme *ts, const int is_family, const int num_elements, t8_element_t *elements[])
{
  /* refine every second element up to the maximum level */
  int level = ts->t8_element_level (elements[0]);
  t8_linearidx_t eid = ts->t8_element_get_linear_id (elements[0], level);
  int maxlevel = *(int *) t8_forest_get_user_data (forest);

  if (eid % 2 && level < maxlevel) {
    return 1;
  }
  return 0;
}

static void
t8_test_gao_check (t8_forest_t forest)
{
  int remote = -1;
  int next_remote = -1;
  int num_remotes;
  int pos = 0;

  t8_locidx_t num_ghost_trees = t8_forest_ghost_num_trees (forest);
  int *remotes = t8_forest_ghost_get_remotes (forest, &num_remotes);

  /* remote stores the remote process of the current ghost element */
  if (num_remotes > 0) {
    remote = remotes[0];
  }
  if (num_remotes > 1) {
    next_remote = remotes[1];
  }
  /* loop over ghost trees */
  for (t8_locidx_t itree = 0, lelement_id = 0; itree < num_ghost_trees; itree++) {
    t8_locidx_t num_elems_in_tree = t8_forest_ghost_tree_num_elements (forest, itree);
    /* Get the global id and element class of this tree */
    t8_gloidx_t gtreeid = t8_forest_ghost_get_global_treeid (forest, itree);
    t8_eclass_t eclass = t8_forest_ghost_get_tree_class (forest, itree);
    for (t8_locidx_t ielem = 0; ielem < num_elems_in_tree; ielem++, lelement_id++) {
      if (pos + 1 < num_remotes && t8_forest_ghost_remote_first_elem (forest, next_remote) <= lelement_id) {
        /* This element belongs to the next remote rank */
        remote = next_remote;
        pos++;
        if (pos + 1 < num_remotes) {
          next_remote = remotes[pos + 1];
        }
      }

      t8_element_t *ghost_element = t8_forest_ghost_get_element (forest, itree, ielem);
      int is_owner = t8_forest_element_check_owner (forest, ghost_element, gtreeid, eclass, remote, 0);
      ASSERT_TRUE (is_owner) << "Owner check for ghost element failed.\n";
      int owner = t8_forest_element_find_owner (forest, gtreeid, ghost_element, eclass);
      ASSERT_EQ (owner, remote) << "Found wrong owner for ghost element.\n";
    }
  }
}

TEST_P (forest_ghost_owner, test_ghost_owner)
{

  /* Compute the minimum level, such that the forest is nonempty */
  int min_level = t8_forest_min_nonempty_level (cmesh, scheme);
  /* start with an empty level */
  min_level = SC_MAX (0, min_level - 1);
  t8_debugf ("Testing ghost exchange with start level %i\n", min_level);
  for (int level = min_level; level < min_level + 3; level++) {
    /* ref the scheme since we reuse it */
    t8_schemexx_ref (scheme);
    /* ref the cmesh since we reuse it */
    t8_cmesh_ref (cmesh);
    /* Create a uniformly refined forest */
    t8_forest_t forest = t8_forest_new_uniform (cmesh, scheme, level, 1, sc_MPI_COMM_WORLD);
    /* Check the owners of the ghost elements */
    t8_test_gao_check (forest);
    /* Adapt the forest and exchange data again */
    int maxlevel = level + 2;
    t8_forest_t forest_adapt = t8_forest_new_adapt (forest, t8_test_gao_adapt, 1, 1, &maxlevel);
    /* Check the owners of the ghost elements */
    t8_test_gao_check (forest_adapt);
    t8_forest_unref (&forest_adapt);
  }
}

INSTANTIATE_TEST_SUITE_P (t8_gtest_ghost_and_owner, forest_ghost_owner, AllCmeshsParam, pretty_print_base_example);
