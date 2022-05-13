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

/* TODO: when this test works for all cmeshes remove if statement in test_cmesh_ghost_exchange_all () */

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

static int
t8_test_exchange_adapt (t8_forest_t forest, t8_forest_t forest_from,
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

/* Construct a data array of uin64_t for all elements and all ghosts,
 * fill the element's entries with their linear id, perform the ghost exchange and
 * check whether the ghost's entries are their linear id.
 */
static void
t8_test_ghost_exchange_data_id (t8_forest_t forest)
{
  t8_eclass_scheme_c *ts;

  t8_locidx_t         num_elements, ielem, num_ghosts, itree;
  t8_linearidx_t      ghost_id, elem_id, ghost_entry;
  t8_element_t       *elem;
  size_t              array_pos = 0;
  sc_array_t          element_data;

  num_elements = t8_forest_get_local_num_elements (forest);
  num_ghosts = t8_forest_get_num_ghosts (forest);
  /* Allocate a uin64_t as data for each element and each ghost */
  sc_array_init_size (&element_data, sizeof (t8_linearidx_t),
                      num_elements + num_ghosts);

  /* Fill the local element entries with their linear id */
  for (itree = 0; itree < t8_forest_get_num_local_trees (forest); itree++) {
    /* Get the eclass scheme for this tree */
    ts = t8_forest_get_eclass_scheme (forest,
                                      t8_forest_get_tree_class (forest,
                                                                itree));
    for (ielem = 0; ielem < t8_forest_get_tree_num_elements (forest, itree);
         ielem++) {
      /* Get a pointer to this element */
      elem = t8_forest_get_element_in_tree (forest, itree, ielem);
      /* Compute the linear id of this element */
      elem_id = ts->t8_element_get_linear_id (elem,
                                              ts->t8_element_level (elem));
      /* Store this id at the element's index in the array */
      *(t8_linearidx_t *) sc_array_index (&element_data, array_pos) = elem_id;
      array_pos++;
    }
  }

  /* Perform the data exchange */
  t8_forest_ghost_exchange_data (forest, &element_data);

  /* We now iterate over all ghost elements and check whether the correct
   * id was received */
  for (itree = 0; itree < t8_forest_get_num_ghost_trees (forest); itree++) {
    /* Get the eclass scheme of this ghost tree */
    ts =
      t8_forest_get_eclass_scheme (forest,
                                   t8_forest_ghost_get_tree_class (forest,
                                                                   itree));
    for (ielem = 0; ielem < t8_forest_ghost_tree_num_elements (forest, itree);
         ielem++) {
      /* Get a pointer to this ghost */
      elem = t8_forest_ghost_get_element (forest, itree, ielem);
      /* Compute its ghost_id */
      ghost_id =
        ts->t8_element_get_linear_id (elem, ts->t8_element_level (elem));
      /* Compare this id with the entry in the element_data array */
      ghost_entry =
        *(t8_linearidx_t *) sc_array_index (&element_data, array_pos);
      SC_CHECK_ABORT (ghost_id == ghost_entry,
                      "Error when exchanging ghost data. Received wrong element id.\n");
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
  sc_array_t          element_data;
  t8_locidx_t         num_elements, ielem, num_ghosts;
  int                 ghost_int;

  num_elements = t8_forest_get_local_num_elements (forest);
  num_ghosts = t8_forest_get_num_ghosts (forest);
  /* Allocate an integer as data for each element and each ghost */
  sc_array_init_size (&element_data, sizeof (int), num_elements + num_ghosts);

  /* Fill the local element entries with the integer 42 */
  for (ielem = 0; ielem < num_elements; ielem++) {
    *(int *) t8_sc_array_index_locidx (&element_data, ielem) = 42;
  }
  /* Perform the ghost data exchange */
  t8_forest_ghost_exchange_data (forest, &element_data);

  /* Check for the ghosts that we received the correct data */
  for (ielem = 0; ielem < num_ghosts; ielem++) {
    /* Get the integer for this ghost */
    ghost_int =
      *(int *) t8_sc_array_index_locidx (&element_data, num_elements + ielem);
    SC_CHECK_ABORT (ghost_int == 42,
                    "Error when exchanging ghost data. Received wrong data.\n");
  }
  /* clean-up */
  sc_array_reset (&element_data);
}

static void
t8_test_ghost_exchange (int cmesh_id)
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
  /* we start with an empty level */
  min_level = SC_MAX (min_level - 1, 0);
  t8_global_productionf
    ("Testing ghost exchange start level %i. cmesh_id = %i\n", min_level,
     cmesh_id);
  for (level = min_level; level < min_level + 3; level++) {
    /* ref the scheme since we reuse it */
    t8_scheme_cxx_ref (scheme);
    /* ref the cmesh since we reuse it */
    t8_cmesh_ref (cmesh);
    /* Create a uniformly refined forest */
    forest = t8_forest_new_uniform (cmesh, scheme, level, 1,
                                    sc_MPI_COMM_WORLD);
    /* exchange ghost data */
    t8_test_ghost_exchange_data_int (forest);
    t8_test_ghost_exchange_data_id (forest);
    /* Adapt the forest and exchange data again */
    maxlevel = level + 2;
    forest_adapt =
      t8_forest_new_adapt (forest, t8_test_exchange_adapt, 1, 1, &maxlevel);
    t8_test_ghost_exchange_data_int (forest_adapt);
    t8_test_ghost_exchange_data_id (forest_adapt);
    t8_forest_unref (&forest_adapt);
  }
  t8_cmesh_destroy (&cmesh);
  t8_scheme_cxx_unref (&scheme);
}

/* The function test_cmesh_ghost_exchange_all () runs the ghost_exchange test for all cmeshes we want to test.
 * We run over all testcases using t8_get_all_testcases() to know how many to check. 
 */
static void
test_cmesh_ghost_exchange_all ()
{
  /* Test all cmeshes over all different inputs we get through their id */
  for (int cmesh_id = 0; cmesh_id < t8_get_number_of_all_testcases ();
       cmesh_id++) {
    /* This if statement is necessary to make the test work by avoiding specific cmeshes which do not work yet for this test.
     * When the issues are gone, remove the if statement. */
    if (cmesh_id != 89 && (cmesh_id < 237 || cmesh_id > 256)) {
      t8_test_ghost_exchange (cmesh_id);
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

  test_cmesh_ghost_exchange_all ();

  sc_finalize ();

  mpiret = sc_MPI_Finalize ();
  SC_CHECK_MPI (mpiret);

  return 0;
}
