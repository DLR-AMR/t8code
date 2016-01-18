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
#include "t8_cmesh_types.h"

T8_EXTERN_C_BEGIN ();

/* allocate a t8_cmesh_tree struct and allocate memory for its entries.
 * No memory for ctrees or ghosts is allocated here */
/* TODO: document */
void                t8_cmesh_trees_init (t8_cmesh_trees_t * ptrees,
                                         int num_procs, t8_topidx_t num_trees,
                                         t8_topidx_t num_ghosts);

/* allocate the first_tree array of a given tree_part in a tree struct
 * with a given number of bytes */
void                t8_cmesh_trees_init_part (t8_cmesh_trees_t trees,
                                              int proc,
                                              t8_locidx_t first_tree,
                                              t8_topidx_t last_tree,
                                              t8_topidx_t num_ghosts,
                                              size_t attr_bytes);

/** Add a tree to a trees structure.
 * \param [in,out]  trees The trees structure to be updated.
 * \param [in]      tree_id The local id of the tree to be inserted.
 * \param [in]      proc  The mpirank of the process from which the tree was
 *                        received.
 * \param [in]      eclass The tree's element class.
 */
void                t8_cmesh_trees_add_tree (t8_cmesh_trees_t trees,
                                             t8_topidx_t tree_id, int proc,
                                             t8_eclass_t eclass);

/* TODO: document. Adds a face-connection */
void                t8_cmesh_tree_set_join (t8_cmesh_trees_t trees,
                                            t8_locidx_t id1, t8_locidx_t id2,
                                            int face1, int face2,
                                            int orientation);

/** Add a ghost to a trees structure.
 * \param [in,out]  trees The trees structure to be updated.
 * \param [in]      ghost_index The local id of the ghost to be inserted.
 * \param [in]      tree_id The global index of the ghost.
 * \param [in]      proc  The mpirank of the process from which the ghost was
 *                        received.
 * \param [in]      eclass The ghost's element class.
 */
void                t8_cmesh_trees_add_ghost (t8_cmesh_trees_t trees,
                                              t8_locidx_t ghost_index,
                                              t8_gloidx_t tree_id, int proc,
                                              t8_eclass_t eclass);

/* TODO: This function return NULL if the tree is not present.
 *       So far no error checking is done here. */
/** Return a pointer to a specific tree in a trees struct.
 * \param [in]      trees The tress structure where the tree is to be looked up.
 * \param [in]      tree  The local id of the tree.
 * \return                A pointer to the tree with local id \a tree.
 */
t8_ctree_t          t8_cmesh_trees_get_tree (t8_cmesh_trees_t trees,
                                             t8_topidx_t tree);

/* TODO: This function return NULL if the ghost is not present.
 *       So far no error checking is done here. */
/** Return a pointer to a specific ghost in a trees struct.
 * \param [in]      trees The tress structure where the tree is to be looked up.
 * \param [in]      ghost The local id of the ghost.
 * \return                A pointer to the ghost with local id \a ghost.
 */
t8_cghost_t         t8_cmesh_trees_get_ghost (t8_cmesh_trees_t trees,
                                              t8_locidx_t ghost);

/* TODO: document */
void                t8_cmesh_trees_init_attributes (t8_cmesh_trees_t trees,
                                                    t8_locidx_t tree_id,
                                                    size_t num_attributes);

/* TODO: These need to be rewritten with package_id and key */
void               *t8_cmesh_trees_get_attribute (t8_cmesh_trees_t trees,
                                                  t8_topidx_t tree_id,
                                                  int package_id, int key,
                                                  size_t * data_size);

/* TODO: These need to be rewritten with package_id and key */
/* TODO: this uses char * and cmesh_set_attribute uses void *. Unify! */
/* attr_tree_index is index of attr in tree's attribute array.
 * We assume that the attributes are already sorted! */
void                t8_cmesh_tree_add_attribute (t8_cmesh_trees_t trees,
                                                 int proc,
                                                 t8_topidx_t tree_id,
                                                 int package_id, int key,
                                                 char *attr, size_t size,
                                                 size_t offset,
                                                 int attr_tree_index);

int                 t8_cmesh_trees_is_equal (t8_cmesh_t cmesh,
                                             t8_cmesh_trees_t trees_a,
                                             t8_cmesh_trees_t trees_b);

void                t8_cmesh_trees_destroy (t8_cmesh_trees_t * trees);

T8_EXTERN_C_END ();

#endif /* !T8_CMESH_PART_TREE_H */
