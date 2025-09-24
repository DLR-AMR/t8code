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

/** \file t8_forest_ghost.h
 * We define the ghost routine to create a layer of halo elements
 * for a forest of trees in this file.
 */

/* TODO: begin documenting this file: make doxygen 2>&1 | grep t8_forest_ghost */

#ifndef T8_FOREST_GHOST_H
#define T8_FOREST_GHOST_H

#include <t8.h>
#include <t8_forest/t8_forest_types.h>

T8_EXTERN_C_BEGIN ();

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

/** 
 * Initialize a ghost type of a forest.
 * 
 * \param[out] pghost     Pointer to the forest's ghost.
 * \param[in]  ghost_type The type of the ghost elements, \see t8_ghost_type_t.
 */
void
t8_forest_ghost_init (t8_forest_ghost_t *pghost, t8_ghost_type_t ghost_type);

/** 
 * Return the number of trees in a ghost.
 *
 * \param[in] forest  The forest.
 * 
 * \return The number of trees in the forest's ghost (or 0 if ghost structure does not exist).
 */
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

/** Given an index in the ghost_tree array, return this tree's number of leaf elements
 * \param [in]  forest      The \a forest. Ghost layer must exist.
 * \param [in]  lghost_tree The ghost tree id of a ghost tree.
 * \return                  The number of ghost leaf elements of the tree.
 * \a forest must be committed before calling this function.
 */
t8_locidx_t
t8_forest_ghost_tree_num_leaf_elements (t8_forest_t forest, t8_locidx_t lghost_tree);

/** Get a pointer to the ghost leaf element array of a ghost tree.
 * \param [in]  forest    The forest. Ghost layer must exist.
 * \param [in]  lghost_tree The ghost tree id of a ghost tree.
 * \return                A pointer to the array of ghost leaf elements of the tree.
 * \a forest must be committed before calling this function.
 */
t8_element_array_t *
t8_forest_ghost_get_tree_leaf_elements (const t8_forest_t forest, const t8_locidx_t lghost_tree);

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

/**
  * Given an index in the ghost_tree array, return this tree's element class.
  * 
  * \param[in] forest       A committed forest.
  * \param[in] lghost_tree  The tree's local index in the ghost_tree array.
  * 
  * \return The element class of the given tree.
 */
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

/** Given an index into the ghost_trees array and for that tree an element index,
 *  return the corresponding element. 
 * \param [in]  forest      The \a forest. Ghost layer must exist.
 * \param [in]  lghost_tree The ghost tree id of a ghost tree.
 * \param [in]  lelement    The local id of the ghost leaf element considered.
 * \return                  A pointer to the ghost leaf element.
 * \a forest must be committed before calling this function.
 */
t8_element_t *
t8_forest_ghost_get_leaf_element (t8_forest_t forest, t8_locidx_t lghost_tree, t8_locidx_t lelement);

/** Return the array of remote ranks.
 * \param [in] forest   A forest with constructed ghost layer.
 * \param [in,out] num_remotes On output the number of remote ranks is stored here.
 * \return              The array of remote ranks in ascending order.
 */
int *
t8_forest_ghost_get_remotes (t8_forest_t forest, int *num_remotes);

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
 * This function is preferred over \ref t8_forest_ghost_unref when it is known that the last reference is to be deleted.
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
t8_forest_ghost_create (t8_forest_t forest);

/** Create one layer of ghost elements for a forest.
 * This version only works with balanced forests and is the original
 * algorithm from p4est: Scalable Algorithms For Parallel Adaptive
 *                Mesh Refinement On Forests of Octrees
 * \param [in,out]    forest     The balanced forest/
 * \a forest must be committed before calling this function.
 * \note The user should prefer \ref t8_forest_ghost_create even for balanced forests.
 */
void
t8_forest_ghost_create_balanced_only (t8_forest_t forest);

/**
 *  Experimental version of \ref t8_forest_ghost_create using the ghost_v3 algorithm 
 */
void
t8_forest_ghost_create_topdown (t8_forest_t forest);

T8_EXTERN_C_END ();

#endif /* !T8_FOREST_GHOST_H */
