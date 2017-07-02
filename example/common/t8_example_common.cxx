/*
  This file is part of t8code.
  t8code is a C library to manage a collection (a forest) of multiple
  connected adaptive space-trees of general element types in parallel.

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

/** \file t8_example_common.cxx Provide routines that are used in more than one
  * example. */

#include <sc_refcount.h>
#include <t8_forest/t8_forest_adapt.h>
#include <t8_element_cxx.hxx>
#include <t8_forest.h>
#include <example/common/t8_example_common.h>

T8_EXTERN_C_BEGIN ();

/* Adapt a forest such that always the second child of the first
 * tree is refined and no other elements. This results in a highly
 * imbalanced forest.
 * The user data of forest must an integer set to the maximum refinement level.
 */
int
t8_common_adapt_balance (t8_forest_t forest, t8_forest_t forest_from,
                         t8_locidx_t which_tree, t8_eclass_scheme_c * ts,
                         int num_elements, t8_element_t * elements[])
{
  int                 level;
  int                 maxlevel, child_id;
  T8_ASSERT (num_elements == 1 || num_elements ==
             ts->t8_element_num_children (elements[0]));
  level = ts->t8_element_level (elements[0]);

  /* we set a maximum refinement level as forest user data */
  maxlevel = *(int *) t8_forest_get_user_data (forest);
  if (level >= maxlevel) {
    /* Do not refine after the maxlevel */
    return 0;
  }
  child_id = ts->t8_element_child_id (elements[0]);
  /* refine the last child of even trees */
  if ((which_tree + t8_forest_get_first_local_tree_id (forest_from)) % 2 == 0
      && child_id == ts->t8_element_num_children (elements[0]) - 1) {
    return 1;
  }
  return 0;
}

#if 0
static int
t8_basic_adapt (t8_forest_t forest, t8_locidx_t which_tree,
                t8_eclass_scheme_c * ts,
                int num_elements, t8_element_t * elements[])
{
  int                 level, mpirank, mpiret;
  T8_ASSERT (num_elements == 1 || num_elements ==
             ts->t8_element_num_children (elements[0]));
  level = ts->t8_element_level (elements[0]);
#if 0
  if (num_elements > 1) {
    /* Do coarsen here */
    if (level > 0)
      return -1;
    return 0;
  }
#endif
  mpiret = sc_MPI_Comm_rank (sc_MPI_COMM_WORLD, &mpirank);
  SC_CHECK_MPI (mpiret);
  if (level < 5)
    /* refine randomly if level is smaller 4 */
    return (unsigned) ((mpirank + 1) * rand ()) % 2;
  return 0;
}
#endif

T8_EXTERN_C_END ();
