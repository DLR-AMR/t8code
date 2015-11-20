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

/** \file t8_cmesh_trees.h
 *
 * TODO: document this file
 */

#ifndef T8_CMESH_PART_TREE_H
#define T8_CMESH_PART_TREE_H

#include <t8.h>
#include <t8_cmesh.h>
#include <t8_cmesh/t8_cmesh_types.h>

typedef t8_part_tree *t8_part_tree_t;
typedef t8_cmesh_trees *t8_cmesh_trees_t;

T8_EXTERN_C_BEGIN ();

/* allocate a t8_cmesh_tree struct and allocate memory for its entries.
 * No memory for ctrees or ghosts is allocated here */
/* TODO: document */
void              t8_cmesh_trees_init (t8_cmesh_trees_t * ptrees, int num_procs,
                                       t8_topidx_t num_trees,
                                       t8_topidx_t num_ghosts);

/* allocate the first_tree array of a given tree_part in a tree struct
 * with a given number of bytes */
void              t8_cmesh_trees_init_part (t8_cmesh_trees_t trees, int proc,
                                            size_t array_size);

t8_ctree_t        t8_cmesh_trees_get_tree (t8_cmesh_trees_t trees,
                                           t8_topidx_t tree);

t8_cghost_t       t8_cmesh_trees_get_ghost (t8_cmesh_trees_t trees,
                                            t8_topidx_t ghost);

void             *t8_cmesh_trees_get_attribute (t8_cmesh_trees_t trees,
                                                t8_topidx_t tree);

T8_EXTERN_C_END ();

#endif /* !T8_CMESH_PART_TREE_H */
