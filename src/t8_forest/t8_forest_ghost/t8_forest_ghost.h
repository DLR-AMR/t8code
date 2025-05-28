/*
  This file is part of t8code.
  t8code is a C library to manage a collection (a forest) of multiple
  connected adaptive space-trees of general element classes in parallel.

  Copyright (C) 2024 the developers

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

/** \file t8_forest_ghost.h
 * Routines to manage and access the ghost layer structure of a forest.
 */

#ifndef T8_FOREST_GHOST_H
#define T8_FOREST_GHOST_H

#include <t8.h>
#include <t8_forest/t8_forest_types.h>
#include <t8_forest/t8_forest_general.h>

T8_EXTERN_C_BEGIN ();

/* The information stored for the remote trees.
 * Each remote process stores an array of these */
typedef struct
{
  t8_gloidx_t global_id;       /* global id of the tree */
  int mpirank;                 /* The mpirank of the remote process */
  t8_element_array_t elements; /* The remote elements of that tree */
  sc_array_t element_indices;  /* The (tree) local indices of the ghost elements. */
  t8_eclass_t eclass;          /* The trees element class */
} t8_ghost_remote_tree_t;

typedef struct
{
  int remote_rank;          /* The rank of the remote process */
  t8_locidx_t num_elements; /* The number of remote elements for this process */
  sc_array_t remote_trees;  /* Array of the remote trees of this process */
} t8_ghost_remote_t;

/* We enumerate the ghost trees by 0, 1, ..., num_ghost_trees - 1
 * In the context of a forest we add the number of local trees as offset,
 * so that we have a range of trees:
 *
 * | 0, 1, ..., num_trees - 1 | num_trees, ..., num_trees + num_ghosts - 1 |
 *
 *      local trees                           ghost trees
 *
 * For the functions in this header an argument lghost_tree always
 * means a number 0 <= lghost_tree < num_ghost_trees - 1
 */

/* TODO: document */
/* returns 0 if ghost structure doesnt exist */
t8_locidx_t
t8_forest_ghost_num_trees (const t8_forest_t forest);

/** Return the element offset of a ghost tree.
 * \param [in]      forest      The forest with constructed ghost layer.
 * \param [in]      lghost_tree A local ghost id of a ghost tree.
 * \return                      The element offset of this ghost tree.
 * \note forest must be committed before calling this function.
 */
t8_locidx_t
t8_forest_ghost_get_tree_element_offset (t8_forest_t forest, t8_locidx_t lghost_tree);

/** Given an index in the ghost_tree array, return this tree's number of elements 
 * \param [in]      forest      The forest with constructed ghost layer.
 * \param [in]      lghost_tree A local ghost id of a ghost tree.
 * \return                      The number of ghost elements in this tree.
*/
t8_locidx_t
t8_forest_ghost_tree_num_elements (t8_forest_t forest, t8_locidx_t lghost_tree);

/** Get a pointer to the ghost element array of a ghost tree.
 * \param [in]  forest    The forest. Ghost layer must exist.
 * \param [in]  lghost_tree The ghost tree id of a ghost tree.
 * \return                A pointer to the array of ghost elements of the tree.
 * \a forest must be committed before calling this function.
 */
t8_element_array_t *
t8_forest_ghost_get_tree_elements (const t8_forest_t forest, const t8_locidx_t lghost_tree);

/** Given a global tree compute the ghost local tree id of it.
 * \param [in]  forest    The forest. Ghost layer must exist.
 * \param [in]  gtreeid   A global tree in \a forest.
 * \return                If \a gtreeid is also a ghost tree, the index in
 *                        the ghost->ghost_trees array of the tree.
 *                        Otherwise a negative number.
 * \a forest must be committed before calling this function.
 * \see https://github.com/DLR-AMR/t8code/wiki/Tree-indexing for more details about tree indexing.
 */
t8_locidx_t
t8_forest_ghost_get_ghost_treeid (t8_forest_t forest, t8_gloidx_t gtreeid);

/* TODO: document */
t8_eclass_t
t8_forest_ghost_get_tree_class (const t8_forest_t forest, const t8_locidx_t lghost_tree);

/** Given a local ghost tree compute the global tree id of it.
 * \param [in]  forest    The forest. Ghost layer must exist.
 * \param [in]  lghost_tree The ghost tree id of a ghost tree.
 * \return                The global id of the local ghost tree \a lghost_tree.
 * \a forest must be committed before calling this function.
 * \see https://github.com/DLR-AMR/t8code/wiki/Tree-indexing for more details about tree indexing.
 */
t8_gloidx_t
t8_forest_ghost_get_global_treeid (const t8_forest_t forest, const t8_locidx_t lghost_tree);

/* TODO: document */
t8_element_t *
t8_forest_ghost_get_element (t8_forest_t forest, t8_locidx_t lghost_tree, t8_locidx_t lelement);

/** Return the first local ghost tree of a remote rank.
 * \param [in] forest   A forest with constructed ghost layer.
 * \param [in] remote   A remote rank of the ghost layer in \a forest.
 * \return              The ghost tree id of the first ghost tree that stores ghost
 *                      elements of \a remote.
 */
t8_locidx_t
t8_forest_ghost_remote_first_tree (t8_forest_t forest, int remote);

/** Return the local index of the first ghost element that belongs to a given remote rank.
 * \param [in] forest   A forest with constructed ghost layer.
 * \param [in] remote   A remote rank of the ghost layer in \a forest.
 * \return              The index i in the ghost elements of the first element of rank \a remote
 */
t8_locidx_t
t8_forest_ghost_remote_first_elem (t8_forest_t forest, int remote);

/** Increase the reference count of a ghost structure.
 * \param [in,out]  ghost     On input, this ghost structure must exist with
 *                            positive reference count.
 */
void
t8_forest_ghost_ref (t8_forest_ghost_t ghost);

/** Decrease the reference count of a ghost structure.
 * If the counter reaches zero, the ghost structure is destroyed.
 * See also \ref t8_forest_ghost_destroy, which is to be preferred when it is
 * known that the last reference to a cmesh is deleted.
 * \param [in,out]  pghost      On input, the ghost structure pointed to must
 *                              exist with positive reference count.
 *                              If the reference count reaches zero, the ghost
 *                              structure is destroyed and this pointer is set to NULL.
 *                              Otherwise, the pointer is not changed.
 */
void
t8_forest_ghost_unref (t8_forest_ghost_t *pghost);

/** Verify that a ghost structure has only one reference left and destroy it.
 * This function is preferred over \ref t8_ghost_unref when it is known that the last reference is to be deleted.
 * \param [in,out]  pghost     This ghost structure must have a reference count of one.
 *                             It can be in any state (committed or not).
 *                             Then it effectively calls \ref t8_forest_ghost_unref.
 */
void
t8_forest_ghost_destroy (t8_forest_ghost_t *pghost);

/** Create one layer of ghost elements for a forest.
 * \see t8_forest_set_ghost
 * \param [in,out]    forest     The forest.
 * \a forest must be committed before calling this function.
 */
void
t8_forest_ghost_create_ext (t8_forest_t forest);

/** Creating one layer of ghost elements for a forest. 
 * experimental version using the ghost_v3 algorithm
 * \param [in,out]    forest     The forest.
 * \a forest must be committed before calling this function.
 * \a forest->ghost_definition must have the \a type FACE and the \a version 1
 */
void
t8_forest_ghost_create_topdown (t8_forest_t forest);

/* Return the remote struct of a given remote rank */
t8_ghost_remote_t *
t8_forest_ghost_get_remote (t8_forest_t forest, int remote);

/** Return the array of remote ranks.
 * \param [in] forest   A forest with constructed ghost layer.
 * \param [in,out] num_remotes On output the number of remote ranks is stored here.
 * \return              The array of remote ranks in ascending order.
 */
int *
t8_forest_ghost_get_remotes (t8_forest_t forest, int *num_remotes);

T8_EXTERN_C_END ();

#endif /* !T8_FOREST_GHOST_H */
