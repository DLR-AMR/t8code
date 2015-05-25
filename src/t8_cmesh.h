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

/** \file t8_cmesh.h
 * We define the coarse mesh of trees in this file.
 */

#ifndef T8_CMESH_H
#define T8_CMESH_H

#include <t8.h>
#include <t8_eclass.h>

typedef struct t8_cmesh *t8_cmesh_t;

T8_EXTERN_C_BEGIN ();

/** Create a new cmesh with reference count one.
 * This cmesh needs to be specialized with the t8_cmesh_set_* calls.
 * Then it needs to be set up with \see t8_cmesh_construct.
 * \param [in,out] pcmesh       On input, this pointer must be non-NULL.
 *                              On return, this pointer set to the new cmesh.
 */
void                t8_cmesh_init (t8_cmesh_t * pcmesh);

/** Set the number of trees and number of trees per eclass for a cmesh.
 * It is not allowed to call this function after \see t8_cmesh_construct.
 * \param [in,out] cmesh        The cmesh to be updated.
 * \param [in]     num_trees    The number of trees to be set.
 * \param [in]     num_trees_per_eclass An array storing for each t8_eclass
 *                              the number of trees of this class.
 */
void                t8_cmesh_set_num_trees (t8_cmesh_t cmesh,
                                            t8_topidx_t num_trees,
                                            const t8_topidx_t
                                            num_trees_per_eclass
                                            [T8_ECLASS_LAST]);

/** Set the number of vertices for a cmesh.
 * TODO: remove vertex functions (we'll think of a safe parallel way).
 * It is not allowed to call this function after \see t8_cmesh_construct.
 * \param [in,out] cmesh        The cmesh to be updated.
 * \param [in]     num_vertices The number of vertices to be set.
 */
void                t8_cmesh_set_num_vertices (t8_cmesh_t cmesh,
                                               t8_topidx_t num_vertices);

/** Insert a new tree into the cmesh.
 * \param [in,out] cmesh        The cmesh to be updated.
 * \param [in]     tree         The global number of the tree.
 * \param [in]     tree_class   The element class of this tree.
 * \param [in]     vertices     An array storing for each corner of the
 *                              tree a vertex in the mesh.  TODO: remove
 */
void                t8_cmesh_insert_tree (t8_cmesh_t cmesh, t8_topidx_t tree,
                                          t8_eclass_t tree_class,
                                          t8_topidx_t * vertices);

/** After allocating and adding properties to a cmesh, finish its construction.
 * \param [in,out] cmesh        Must be created with \see t8_cmesh_init and
 *                              specialized with t8_cmesh_set_* calls first.
 */
void                t8_cmesh_construct (t8_cmesh_t cmesh);

/** Increase the reference counter of a cmesh.
 * \param [in,out] cmesh        On input, this cmesh must exist with positive
 *                              reference count.  It may be in any state.
 */
void                t8_cmesh_ref (t8_cmesh_t cmesh);

/** Decrease the reference counter of a cmesh.
 * If the counter reaches zero, this cmesh is destroyed.
 * \param [in,out] pcmesh       On input, the cmesh pointed to must exist
 *                              with positive reference count.  It may be in
 *                              any state.  If the reference count reaches
 *                              zero, the cmesh is destroyed and this pointer
 *                              set to NULL.
 *                              Otherwise, the pointer is not changed and
 *                              the cmesh is not modified in other ways.
 */
void                t8_cmesh_unref (t8_cmesh_t * pcmesh);

/** Create a coarse mesh that consists of a single triangle.
 * \return          A valid cmesh, as if _init and _construct had been called.
 */
t8_cmesh_t          t8_cmesh_new_tri (void);

/** Create a coarse mesh that consists of a single tetrahedron.
 * \return          A valid cmesh, as if _init and _construct had been called.
 */
t8_cmesh_t          t8_cmesh_new_tet (void);

/* TODO: add new_quad */

/** Create a coarse mesh that consists of a single hexahedron.
 * \return          A valid cmesh, as if _init and _construct had been called.
 */
t8_cmesh_t          t8_cmesh_new_hex (void);

T8_EXTERN_C_END ();

#endif /* !T8_CMESH_H */
