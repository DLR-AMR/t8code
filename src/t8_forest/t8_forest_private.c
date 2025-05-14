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

const t8_element_t*
t8_forest_get_tree_leaf_element (t8_tree_t tree, t8_locidx_t elem_in_tree)
{
  T8_ASSERT (tree != NULL);
  T8_ASSERT (0 <= elem_in_tree && elem_in_tree < t8_forest_get_tree_leaf_element_count (tree));
  return t8_element_array_index_locidx (&tree->elements, elem_in_tree);
}

t8_element_t*
t8_forest_get_tree_leaf_element_mutable (t8_tree_t tree, t8_locidx_t elem_in_tree)
{
  return (t8_element_t*) t8_forest_get_tree_leaf_element (tree, elem_in_tree);
}

const t8_element_array_t*
t8_forest_get_tree_leaf_element_array (const t8_forest_t forest, t8_locidx_t ltreeid)
{
  T8_ASSERT (t8_forest_is_committed (forest));
  T8_ASSERT (0 <= ltreeid && ltreeid < t8_forest_get_num_local_trees (forest));

  return &t8_forest_get_tree (forest, ltreeid)->elements;
}

t8_element_array_t*
t8_forest_get_tree_leaf_element_array_mutable (const t8_forest_t forest, t8_locidx_t ltreeid)
{
  return (t8_element_array_t*) t8_forest_get_tree_leaf_element_array (forest, ltreeid);
}
