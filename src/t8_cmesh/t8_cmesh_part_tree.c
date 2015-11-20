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

/** \file t8_cmesh_part_tree.c
 *
 * TODO: document this file
 */

#include <t8_cmesh/t8_cmesh_part_tree.h>

t8_ctree_t
t8_part_tree_get_tree (t8_part_tree_t P,t8_topidx_t tree)
{
  T8_ASSERT (0 <= tree && tree < P->num_trees);
  return ((t8_ctree_t) P->first_tree) + tree;
}

t8_cghost_t
t8_part_tree_get_ghost (t8_part_tree_t P,t8_topidx_t ghost)
{
  t8_cghost_t         first_ghost;

  T8_ASSERT (0 <= ghost && ghost < P->num_ghosts);
  first_ghost = (t8_cghost_t)
      (P->first_tree + (P->num_trees * sizeof (t8_ctree_struct_t)));
  return first_ghost + ghost;
}

void *
t8_part_tree_get_attribute (t8_part_tree_t P,size_t offset)
{
  return P->first_tree + (P->num_trees * sizeof (t8_ctree_struct_t)) +
      (P->num_ghosts * sizeof (t8_cghost_struct_t)) + offset;
}
