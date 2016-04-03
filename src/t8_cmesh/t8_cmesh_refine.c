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

/** \file t8_cmesh_refine.c
 *
 * TODO: document this file
 */

/* TODO: could this file be part of cmesh_commit.c? */

#include <t8_cmesh.h>
#include "t8_cmesh_types.h"
#include "t8_cmesh_trees.h"
#include "t8_cmesh_partition.h"


void            t8_cmesh_refine (t8_cmesh_t cmesh)
{
  t8_cmesh_t          cmesh_from;
  t8_locidx_t         itree;
  t8_ctree_t          tree, newtrees;
  t8_locidx_t        *tree_neighbors;
  int8_t             *ttf;
  int                 dim, factor, factor_ghosts, level;

  T8_ASSERT (cmesh != NULL);
  T8_ASSERT (cmesh->set_from != NULL);
  T8_ASSERT (cmesh->from_method == T8_CMESH_FROM_REFINE);
  T8_ASSERT (cmesh->set_from->committed);
  T8_ASSERT (cmesh->set_from->num_trees_per_eclass[T8_ECLASS_PYRAMID] == 0);

  cmesh_from = (t8_cmesh_t) cmesh->set_from;
  dim = cmesh_from->dimension;
  level = cmesh_from->set_level;
  /* The number of new trees per old tree
   * dim     factor (level = 1)
   *  0         1   (points)
   *  1         2   (lines)
   *  2         4   (quads and triangles)
   *  3         8   (Hexes, prisms, Tets)
   */
  factor = 1 << (dim * level); /
  /* Since we only consider face-ghosts, the numer of new ghosts per old ghosts is the number of
   * face-children of a face. */
  /* The number of new ghosts per old ghosts
   * dim     factor (level = 1)
   *  0         0   (points -> No ghosts)
   *  1         1   (lines -> boundaries are points)
   *  2         2   (quads and triangles -> boundaries are lines)
   *  3         4   (Hexes,prisms and Tets -> boundaries are quads/triangles)
   */

  factor_ghosts = 1 << ((dim - 1) * level); /* The number of new ghosts per old ghosts */

  cmesh->num_local_trees = cmesh_from->num_local_trees * factor;
  cmesh->num_trees = cmesh_from->num_trees * factor;
  cmesh->first_tree = cmesh_from->first_tree * factor;
  /* Check for locidx overflow */
  T8_ASSERT ((t8_gloidx_t) cmesh_from->num_local_trees * factor ==
             cmesh->num_local_trees);
  /************************/
  /* Create the new trees */
  /************************/
  /* Initialize trees struct with yet unknown number of ghosts */
  t8_cmesh_trees_init (cmesh->trees, 1, cmesh->num_local_trees, 0);
  for (itree = 0;itree < cmesh_from->num_local_trees;itree++) {
    tree = t8_cmesh_trees_get_tree_ext (cmesh_from->trees, itree,
                                        &tree_neighbors, &ttf);
  }
}
