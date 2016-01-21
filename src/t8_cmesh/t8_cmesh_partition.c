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

/** \file t8_cmesh_partition.c
 *
 * TODO: document this file
 */

#include <t8_cmesh.h>
#include "t8_cmesh_types.h"
#include "t8_cmesh_trees.h"
#include "t8_cmesh_partition.h"

void
t8_cmesh_partition (t8_cmesh_t cmesh)
{
  t8_cmesh_t        cmesh_from;
  t8_gloidx_t       last_tree;

  T8_ASSERT (cmesh != NULL);
  T8_ASSERT (!cmesh->committed);
  T8_ASSERT (cmesh->set_partitioned);
  T8_ASSERT (cmesh->set_from != NULL);
  T8_ASSERT (cmesh->set_from->committed);

  cmesh_from = cmesh->set_from;
  if (cmesh->set_level >= 0) {
    /* Compute first and last tree index */
    t8_cmesh_uniform_bounds (cmesh_from, cmesh->set_level, &cmesh->first_tree,
                             NULL, &last_tree, NULL, &cmesh->last_tree_shared);
    cmesh->num_local_trees = last_tree - cmesh->first_tree + 1;
  }
  else {
    T8_ASSERT (cmesh->tree_per_proc != NULL);
  }
}
