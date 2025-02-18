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
#include <t8_schemes/t8_scheme.hxx>
#include <t8_forest/t8_forest_general.h>
#include <t8_forest/t8_forest_geometrical.h>
#include <example/common/t8_example_common.h>

T8_EXTERN_C_BEGIN ();

/* Adapt a forest such that always the second child of the first
 * tree is refined and no other elements. This results in a highly
 * imbalanced forest.
 * The user data of forest must an integer set to the maximum refinement level.
 */
int
t8_common_adapt_balance (t8_forest_t forest, t8_forest_t forest_from, t8_locidx_t which_tree,
                         const t8_eclass_t tree_class, [[maybe_unused]] t8_locidx_t lelement_id, const t8_scheme *scheme,
                         [[maybe_unused]] const int is_family, [[maybe_unused]] const int num_elements, t8_element_t *elements[])
{
  int level;
  int maxlevel, child_id;
  T8_ASSERT (!is_family || num_elements == scheme->element_get_num_children (tree_class, elements[0]));
  level = scheme->element_get_level (tree_class, elements[0]);

  /* we set a maximum refinement level as forest user data */
  maxlevel = *(int *) t8_forest_get_user_data (forest);
  if (level >= maxlevel) {
    /* Do not refine after the maxlevel */
    return 0;
  }
  child_id = scheme->element_get_child_id (tree_class, elements[0]);
  /* refine the last child of even trees */
  if ((which_tree + t8_forest_get_first_local_tree_id (forest_from)) % 2 == 0
      && child_id == scheme->element_get_num_children (tree_class, elements[0]) - 1) {
    return 1;
  }
  return 0;
}

int
t8_common_within_levelset (t8_forest_t forest, const t8_locidx_t ltreeid, const t8_element_t *element,
                           t8_example_level_set_fn levelset, double band_width, double t, void *udata)
{
  double elem_midpoint[3], elem_diam;
  double value;
  const t8_eclass_t tree_class = t8_forest_get_eclass (forest, ltreeid);
  const t8_scheme *scheme = t8_forest_get_scheme (forest);

  T8_ASSERT (band_width >= 0);
  if (band_width == 0) {
    /* If bandwidth = 0, we only refine the elements that are intersected by the zero level-set */
    const int num_corners = scheme->element_get_num_corners (tree_class, element);
    int sign = 1, icorner;
    double coords[3];

    /* Compute LS function at first corner */
    t8_forest_element_coordinate (forest, ltreeid, element, 0, coords);
    /* compute the level-set function at this corner */
    value = levelset (coords, t, udata);
    /* sign = 1 if value > 0, -1 if value < 0, 0 if value = 0 */
    sign = value > 0 ? 1 : -(value < 0);
    /* iterate over all corners */
    for (icorner = 1; icorner < num_corners; icorner++) {
      t8_forest_element_coordinate (forest, ltreeid, element, icorner, coords);
      /* compute the level-set function at this corner */
      value = levelset (coords, t, udata);
      if ((value > 0 && sign <= 0) || (value == 0 && sign != 0) || (value < 0 && sign >= 0)) {
        /* The sign of the LS function changes across the element, we refine it */
        return 1;
      }
    }
    return 0;
  }

  /* Compute the coordinates of the anchor node X. */
  t8_forest_element_centroid (forest, ltreeid, element, elem_midpoint);
  /* Compute the element's diameter */
  elem_diam = t8_forest_element_diam (forest, ltreeid, element);
  /* Compute L(X) */
  value = levelset (elem_midpoint, t, udata);

  if (fabs (value) < band_width * elem_diam) {
    /* The element is in the band that should be refined. */
    return 1;
  }
  return 0;
}

/** Adapt a forest along a given level-set function.
 * The user data of forest must be a pointer to a \a t8_example_level_set_struct_t.
 * An element in the forest is refined, if it is in a band of \a band_with many
 * \a max_level elements around the zero level-set Gamma = { x | L(x) = 0}
 */
/* TODO: Currently the band_width control is not working yet. */
int
t8_common_adapt_level_set (t8_forest_t forest, t8_forest_t forest_from, t8_locidx_t which_tree,
                           const t8_eclass_t tree_class, [[maybe_unused]] t8_locidx_t lelement_id, const t8_scheme *scheme,
                           const int is_family, [[maybe_unused]] const int num_elements, t8_element_t *elements[])
{
  t8_example_level_set_struct_t *data;
  int within_band;
  int level;

  T8_ASSERT (!is_family || num_elements == scheme->element_get_num_children (tree_class, elements[0]));

  data = (t8_example_level_set_struct_t *) t8_forest_get_user_data (forest);
  level = scheme->element_get_level (tree_class, elements[0]);

  /* Get the minimum and maximum x-coordinate from the user data pointer of forest */
  data = (t8_example_level_set_struct_t *) t8_forest_get_user_data (forest);

  /* If maxlevel is exceeded then coarsen */
  if (level > data->max_level && is_family) {
    return -1;
  }
  /* Refine at least until min level */
  if (level < data->min_level) {
    return 1;
  }
  within_band = t8_common_within_levelset (forest_from, which_tree, elements[0], data->L, data->band_width / 2, data->t,
                                           data->udata);
  if (within_band && level < data->max_level) {
    /* The element can be refined and lies inside the refinement region */
    return 1;
  }
  else if (is_family && level > data->min_level && !within_band) {
    /* If element lies out of the refinement region and a family was given
     * as argument, we coarsen to level base level */
    return -1;
  }
  return 0;
}

T8_EXTERN_C_END ();
