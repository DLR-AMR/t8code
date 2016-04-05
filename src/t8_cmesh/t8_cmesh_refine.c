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


void
t8_cmesh_refine_tree (t8_cmesh_t cmesh, t8_cmesh_t cmesh_from,
                      t8_locidx_t treeid, t8_locidx_t firstnewtree)
{
  t8_ctree_t          tree, newtree;
  t8_locidx_t         itree;
  t8_locidx_t        *tree_neighbors, *ntree_neighbors;
  int8_t             *ttf, *nttf;
  int                 dim, factor, level;


  dim = cmesh_from->dimension;
  level = cmesh->set_level;
  T8_ASSERT (level == 1);  /* levels bigger than 1 are not yet implemented */
  /* The number of new trees per old tree
   * dim     factor (level = 1)
   *  0         1   (points)
   *  1         2   (lines)
   *  2         4   (quads and triangles)
   *  3         8   (Hexes, prisms, Tets)
   */
  factor = 1 << (dim * level);
  tree = t8_cmesh_trees_get_tree_ext (cmesh_from->trees, treeid,
                                      &tree_neighbors, &ttf);
  for (itree = 0;itree < factor;itree++) {
    newtree = t8_cmesh_trees_get_tree_ext (cmesh, firstnewtree + itree,
                                           &tree_neighbors, &nttf);
    newtree->eclass = tree->eclass;
    newtree->num_attributes = tree->num_attributes;
    newtree->treeid = firstnewtree + itree;
  }
}

void            t8_cmesh_refine (t8_cmesh_t cmesh)
{
  t8_cmesh_t          cmesh_from;
  t8_locidx_t         itree, firstnewtree;
  t8_ctree_t          tree;
  int                 dim, factor, factor_ghosts, level;

  T8_ASSERT (cmesh != NULL);
  T8_ASSERT (cmesh->set_from != NULL);
  T8_ASSERT (cmesh->from_method == T8_CMESH_FROM_REFINE);
  T8_ASSERT (cmesh->set_from->committed);
  T8_ASSERT (cmesh->set_from->num_trees_per_eclass[T8_ECLASS_PYRAMID] == 0);
  T8_ASSERT (cmesh->set_level == 1); /* levels bigger than 1 are not yet implemented */

  cmesh_from = (t8_cmesh_t) cmesh->set_from;
  dim = cmesh_from->dimension;
  level = cmesh->set_level;
  /* The number of new trees per old tree
   * dim     factor (level = 1)
   *  0         1   (points)
   *  1         2   (lines)
   *  2         4   (quads and triangles)
   *  3         8   (Hexes, prisms, Tets)
   */
  factor = 1 << (dim * level);
  cmesh->num_local_trees = cmesh_from->num_local_trees * factor;
  cmesh->num_trees = cmesh_from->num_trees * factor;
  /* Since we only consider face-ghosts, the numer of new ghosts per old ghosts is the number of
   * face-children of a face. */
  /* The number of new ghosts per old ghosts
   * dim     factor (level = 1)
   *  0         0   (points -> No ghosts)
   *  1         1   (lines -> boundaries are points)
   *  2         2   (quads and triangles -> boundaries are lines)
   *  3         4   (Hexes,prismsg and Tets -> boundaries are quads/triangles)
   */

  factor_ghosts = 1 << ((dim - 1) * level); /* The number of new ghosts per old ghosts */
  cmesh->num_ghosts = cmesh_from->num_ghosts * factor_ghosts;
  cmesh->first_tree = cmesh_from->first_tree * factor;
  /* Check for locidx overflow */
  T8_ASSERT ((t8_gloidx_t) cmesh_from->num_local_trees * factor ==
             cmesh->num_local_trees);
  /************************/
  /* Create the new trees */
  /************************/
  t8_cmesh_trees_init (cmesh->trees, 1, cmesh->num_local_trees,
                       cmesh->num_ghosts);
  /* Loop over all trees in cmesh_from and refine them. */
  for (itree = 0, firstnewtree = 0;itree < cmesh_from->num_local_trees;
       itree++, firstnewtree += factor) {
    t8_cmesh_refine_tree (cmesh, cmesh_from, itree, firstnewtree);
  }
}
