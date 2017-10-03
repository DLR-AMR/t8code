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
                         t8_locidx_t which_tree, t8_locidx_t lelement_id,
                         t8_eclass_scheme_c * ts, int num_elements,
                         t8_element_t * elements[])
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

/* Get the coordinates of the anchor node of an element */
void
t8_common_midpoint (t8_forest_t forest, t8_locidx_t which_tree,
                    t8_eclass_scheme_c * ts, t8_element_t * element,
                    double elem_midpoint_f[3])
{
  double             *tree_vertices;

  tree_vertices = t8_cmesh_get_tree_vertices (t8_forest_get_cmesh (forest),
                                              t8_forest_ltreeid_to_cmesh_ltreeid
                                              (forest, which_tree));

  t8_forest_element_centroid (forest, which_tree, element, tree_vertices,
                              elem_midpoint_f);
}

/** Adapt a forest along a given level-set function.
 * The user data of forest must be a pointer to a \a t8_example_level_set_struct_t.
 * An element in the forest is refined, if it is in a band of \a band_with many
 * \a max_level elements around the zero level-set Gamma = { x | L(x) = 0}
 */
/* TODO: Currently the band_width control is not working yet. */
int
t8_common_adapt_level_set (t8_forest_t forest,
                           t8_forest_t forest_from,
                           t8_locidx_t which_tree,
                           t8_eclass_scheme_c * ts,
                           int num_elements, t8_element_t * elements[])
{
  t8_example_level_set_struct_t *data;

  data = (t8_example_level_set_struct_t *) t8_forest_get_user_data (forest);
  t8_example_level_set_fn L;
  int                 level, min_level, max_level;
  double              elem_midpoint[3];
  double              value;

  T8_ASSERT (num_elements == 1 || num_elements ==
             ts->t8_element_num_children (elements[0]));
  level = ts->t8_element_level (elements[0]);

  /* Get the minimum and maximum x-coordinate from the user data pointer of forest */
  data = (t8_example_level_set_struct_t *) t8_forest_get_user_data (forest);
  min_level = data->min_level;
  max_level = data->max_level;
  L = data->L;
  /* Compute the coordinates of the anchor node X. */
  t8_common_midpoint (forest_from, which_tree, ts,
                      elements[0], elem_midpoint);

  /* Compute L(X) */
  value =
    L (elem_midpoint[0], elem_midpoint[1], elem_midpoint[2], data->udata);
#if 1
  if (value >= -data->band_width / (5 * level)
      && value < data->band_width / (5 * level) && level < max_level) {
#else
  if (value >= -0.1 && value < 0.1 && level < max_level) {
#endif
    /* The element is in the band that should be refined. */
    return 1;
  }
  else if (num_elements > 1 && level > min_level) {
    /* If element lies out of the refinement region and a family was given
     * as argument, we coarsen to level base level */
    return -1;
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
