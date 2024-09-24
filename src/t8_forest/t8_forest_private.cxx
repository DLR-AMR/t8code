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

#include <t8_forest/t8_forest_types.h>
#include <t8_forest/t8_forest_private.h>
#include <t8_forest/t8_forest_general.h>
#include <t8_element.hxx>

T8_EXTERN_C_BEGIN ();

const t8_element_t*
t8_forest_get_tree_element (t8_tree_t tree, t8_locidx_t elem_in_tree)
{
  T8_ASSERT (tree != NULL);
  T8_ASSERT (0 <= elem_in_tree && elem_in_tree < t8_forest_get_tree_element_count (tree));
  return t8_element_array_index_locidx (&tree->elements, elem_in_tree);
}

t8_element_t*
t8_forest_get_tree_element_mutable (t8_tree_t tree, t8_locidx_t elem_in_tree)
{
  return (t8_element_t*) t8_forest_get_tree_element (tree, elem_in_tree);
}

const t8_element_array_t*
t8_forest_get_tree_element_array (const t8_forest_t forest, t8_locidx_t ltreeid)
{
  T8_ASSERT (t8_forest_is_committed (forest));
  T8_ASSERT (0 <= ltreeid && ltreeid < t8_forest_get_num_local_trees (forest));

  return &t8_forest_get_tree (forest, ltreeid)->elements;
}

t8_element_array_t*
t8_forest_get_tree_element_array_mutable (const t8_forest_t forest, t8_locidx_t ltreeid)
{
  return (t8_element_array_t*) t8_forest_get_tree_element_array (forest, ltreeid);
}

t8_locidx_t
t8_forest_bin_search_lower (const t8_element_array_t* elements, const t8_linearidx_t element_id, const int maxlevel)
{
  t8_linearidx_t query_id;
  t8_locidx_t low, high, guess;

  const t8_eclass_scheme_c* ts = t8_element_array_get_scheme (elements);
  /* At first, we check whether any element has smaller id than the
   * given one. */
  const t8_element_t* query = t8_element_array_index_int (elements, 0);
  query_id = ts->t8_element_get_linear_id (query, maxlevel);
  if (query_id > element_id) {
    /* No element has id smaller than the given one */
    return -1;
  }

  /* We now perform the binary search */
  low = 0;
  high = t8_element_array_get_count (elements) - 1;
  while (low < high) {
    guess = (low + high + 1) / 2;
    query = t8_element_array_index_int (elements, guess);
    query_id = ts->t8_element_get_linear_id (query, maxlevel);
    if (query_id == element_id) {
      /* we are done */
      return guess;
    }
    else if (query_id > element_id) {
      /* look further left */
      high = guess - 1;
    }
    else {
      /* look further right, but keep guess in the search range */
      low = guess;
    }
  }
  T8_ASSERT (low == high);
  return low;
}

T8_EXTERN_C_END ();
