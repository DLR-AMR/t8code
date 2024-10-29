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

/* TODO: comment */
void
t8_forest_ghost_init (t8_forest_ghost_t *pghost, t8_ghost_type_t ghost_type);

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

/* TODO: document */
t8_locidx_t
t8_forest_ghost_tree_num_elements (t8_forest_t forest, t8_locidx_t lghost_tree);

/** Retrieves a ghost element from its ghost tree given the element's linear id.
 *
 * \param [in] forest The forest object.
 * \param [in] lghost_tree The local index of the ghost tree.
 * \param [in] linear_id The linear id of the element.
 * \param [in] element_level The level of the element.
 * \param [out] loc_ghost_id The local id of the ghost. -1 if no ghost was found.
 * \return The ghost element. nullptr if no ghost was found.
 */
const t8_element_t *
t8_ghost_get_ghost_in_tree_from_linear_id (t8_forest_t forest, t8_locidx_t lghost_tree, t8_linearidx_t linear_id,
                                           int element_level, t8_locidx_t *loc_ghost_id);

/** Retrieves the local index of a ghost element in its specific ghost tree.
 *
 * \param [in] forest The forest object.
 * \param [in] lghost_tree The local index of the ghost tree.
 * \param [in] ghost_element The ghost element.
 * \return The local index of the element in the ghost element array of the ghost tree.
 *         -1 if no ghost element was found.
 */
t8_locidx_t
t8_ghost_get_ghost_id_in_tree (t8_forest_t forest, t8_locidx_t lghost_tree, t8_element_t *ghost_element);

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
 * This function is preferred over \ref t8_ghost_unref when it is known that the last reference is to be deleted.
 * \param [in,out]  pghost     This ghost structure must have a reference count of one.
 *                             It can be in any state (committed or not).
 *                             Then it effectively calls \ref t8_forest_ghost_unref.
 */
void
t8_forest_ghost_destroy (t8_forest_ghost_t *pghost);

/** Part of step 2 of the ghost_creat_ext 
 * for ghost_type face
 * Is declared, so that ghost_interface_face can use it
 * \see t8_forest_ghost_interface_faces::t8_ghost_step_2
 * \param [in,out]    forest     The forest.
 */
void
t8_forest_ghost_fill_remote_v3 (t8_forest_t forest);

/** Part of step 2 of the ghost_creat_ext 
 * for ghost_type face
 * Is declared, so that ghost_interface_face can use it
 * \see t8_forest_ghost_interface_faces::t8_ghost_step_2
 * \param [in,out]    forest     The forest.
 */
void
t8_forest_ghost_fill_remote (t8_forest_t forest, t8_forest_ghost_t ghost, int ghost_method);

/** Create one layer of ghost elements for a forest.
 * \see t8_forest_set_ghost
 * \param [in,out]    forest     The forest.
 * \a forest must be committed before calling this function.
 */
void
t8_forest_ghost_create_ext (t8_forest_t forest);

/** Create one layer of ghost elements for a forest.
 * \see t8_forest_set_ghost
 * \param [in,out]    forest     The forest.
 * \a forest must be committed before calling this function.
 * \a forest->ghost_interface must have the \a type FACE and the \a version 2
 */
void
t8_forest_ghost_create (t8_forest_t forest);

/** Create one layer of ghost elements for a forest.
 * This version only works with balanced forests and is the original
 * algorithm from p4est: Scalable Algorithms For Parallel Adaptive
 *                Mesh Refinement On Forests of Octrees
 * \param [in,out]    forest     The balanced forest/
 * \a forest must be committed before calling this function.
 * \a forest->ghost_interface must have the \a type FACE and the \a version 1
 * \note The user should prefer \ref t8_forest_ghost_create even for balanced forests.
 */
void
t8_forest_ghost_create_balanced_only (t8_forest_t forest);

/** Creating one layer of ghost elements for a forest. 
 * experimental version using the ghost_v3 algorithm
 * \param [in,out]    forest     The forest.
 * \a forest must be committed before calling this function.
 * \a forest->ghost_interface must have the \a type FACE and the \a version 1
 */
void
t8_forest_ghost_create_topdown (t8_forest_t forest);

T8_EXTERN_C_END ();

#endif /* !T8_FOREST_GHOST_H */
