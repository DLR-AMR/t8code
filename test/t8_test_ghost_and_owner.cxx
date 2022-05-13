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

#include <t8_eclass.h>
#include <t8_schemes/t8_default/t8_default_cxx.hxx>
#include <t8_forest.h>
#include <t8_forest/t8_forest_ghost.h>
#include <t8_forest/t8_forest_private.h>
#include <t8_cmesh.h>
#include "t8_cmesh/t8_cmesh_testcases.h"

/* TODO: when this test works for all cmeshes remove if statement in test_cmesh_ghost_and_owner_all () */

/* This test program tests the forest ghost layer.
 * We adapt a forest and create its ghost layer. Afterwards, we
 * parse through all ghost elements and test whether the owner of an
 * element is in face the owner that is stored in the ghost layer.
  */

static int
t8_test_gao_adapt (t8_forest_t forest, t8_forest_t forest_from,
                   t8_locidx_t which_tree, t8_locidx_t lelement_id,
                   t8_eclass_scheme_c *ts, const int is_family,
                   const int num_elements, t8_element_t *elements[])
{
  t8_linearidx_t      eid;
  int                 level, maxlevel;

  /* refine every second element up to the maximum level */
  level = ts->t8_element_level (elements[0]);
  eid = ts->t8_element_get_linear_id (elements[0], level);
  maxlevel = *(int *) t8_forest_get_user_data (forest);

  if (eid % 2 && level < maxlevel) {
    return 1;
  }
  return 0;
}

static void
t8_test_gao_check (t8_forest_t forest)
{
  t8_locidx_t         num_ghost_trees, num_elems_in_tree;
  t8_locidx_t         itree, ielem, lelement_id;
  t8_element_t       *ghost_element;
  t8_gloidx_t         gtreeid;
  t8_eclass_t         eclass;
  int                 remote = -1, next_remote = -1, num_remotes, *remotes;
  int                 pos = 0, is_owner, owner;

  num_ghost_trees = t8_forest_ghost_num_trees (forest);
  remotes = t8_forest_ghost_get_remotes (forest, &num_remotes);
  pos = 0;
  /* remote stores the remote process of the current ghost element */
  if (num_remotes > 0) {
    remote = remotes[0];
  }
  if (num_remotes > 1) {
    next_remote = remotes[1];
  }
  /* loop over ghost trees */
  for (itree = 0, lelement_id = 0; itree < num_ghost_trees; itree++) {
    num_elems_in_tree = t8_forest_ghost_tree_num_elements (forest, itree);
    /* Get the global id and element class of this tree */
    gtreeid = t8_forest_ghost_get_global_treeid (forest, itree);
    eclass = t8_forest_ghost_get_tree_class (forest, itree);
    for (ielem = 0; ielem < num_elems_in_tree; ielem++, lelement_id++) {
      if (pos + 1 < num_remotes
          && t8_forest_ghost_remote_first_elem (forest,
                                                next_remote) <= lelement_id) {
        /* This element belongs to the next remote rank */
        remote = next_remote;
        pos++;
        if (pos + 1 < num_remotes) {
          next_remote = remotes[pos + 1];
        }
      }

      ghost_element = t8_forest_ghost_get_element (forest, itree, ielem);
      is_owner =
        t8_forest_element_check_owner (forest, ghost_element, gtreeid, eclass,
                                       remote, 0);
      SC_CHECK_ABORT (is_owner, "Owner check for ghost element failed.\n");
      owner =
        t8_forest_element_find_owner (forest, gtreeid, ghost_element, eclass);
      SC_CHECK_ABORT (owner == remote,
                      "Found wrong owner for ghost element.\n");
    }
  }
}

static void
t8_test_ghost_owner (int cmesh_id)
{
  int                 level, min_level, maxlevel;
  t8_cmesh_t          cmesh;
  t8_forest_t         forest, forest_adapt;
  t8_scheme_cxx_t    *scheme;

  scheme = t8_scheme_new_default_cxx ();
  /* Construct a cmesh */
  cmesh = t8_test_create_cmesh (cmesh_id);
  /* Compute the minimum level, such that the forest is nonempty */
  min_level = t8_forest_min_nonempty_level (cmesh, scheme);
  /* start with an empty level */
  min_level = SC_MAX (0, min_level - 1);
  t8_global_productionf
    ("Testing ghost exchange with start level %i\n", min_level);
  for (level = min_level; level < min_level + 3; level++) {
    /* ref the scheme since we reuse it */
    t8_scheme_cxx_ref (scheme);
    /* ref the cmesh since we reuse it */
    t8_cmesh_ref (cmesh);
    /* Create a uniformly refined forest */
    forest = t8_forest_new_uniform (cmesh, scheme, level, 1,
                                    sc_MPI_COMM_WORLD);
    /* Check the owners of the ghost elements */
    t8_test_gao_check (forest);
    /* Adapt the forest and exchange data again */
    maxlevel = level + 2;
    forest_adapt =
      t8_forest_new_adapt (forest, t8_test_gao_adapt, 1, 1, &maxlevel);
    /* Check the owners of the ghost elements */
    t8_test_gao_check (forest_adapt);
    t8_forest_unref (&forest_adapt);
  }
  t8_cmesh_destroy (&cmesh);

  t8_scheme_cxx_unref (&scheme);
}

/* The function test_cmesh_ghost_exchange_all () runs the ghost_exchange test for all cmeshes we want to test.
 * We run over all testcases using t8_get_all_testcases() to know how many to check. 
 */
static void
test_cmesh_ghost_owner_all ()
{
  /* Test all cmeshes over all different inputs we get through their id */
  for (int cmesh_id = 0; cmesh_id < t8_get_number_of_all_testcases ();
       cmesh_id++) {
    /* This if statement is necessary to make the test work by avoiding specific cmeshes which do not work yet for this test.
     * When the issues are gone, remove the if statement. */
    if (cmesh_id != 89 && (cmesh_id < 237 || cmesh_id > 256)) {
      t8_test_ghost_owner (cmesh_id);
    }
  }
}

int
main (int argc, char **argv)
{
  int                 mpiret;
  sc_MPI_Comm         mpic;

  mpiret = sc_MPI_Init (&argc, &argv);
  SC_CHECK_MPI (mpiret);

  mpic = sc_MPI_COMM_WORLD;
  sc_init (mpic, 1, 1, NULL, SC_LP_PRODUCTION);
  p4est_init (NULL, SC_LP_ESSENTIAL);
  t8_init (SC_LP_DEFAULT);

  test_cmesh_ghost_owner_all ();
  t8_debugf ("Test successful\n");

  sc_finalize ();

  mpiret = sc_MPI_Finalize ();
  SC_CHECK_MPI (mpiret);

  return 0;
}
