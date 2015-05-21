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

#include <sc_refcount.h>
#include <t8.h>
#include <t8_cmesh.h>

/** \file t8_cmesh.h
 *
 * TODO: document this file
 */

typedef struct t8_cmesh
{
  sc_refcount_t       rc; /**< The reference count of the cmesh. */
  t8_topidx_t         num_vertices; /**< the number of vertices that define the \a embedding of the forest (not the topology) */
  t8_topidx_t         num_trees;  /**< The number of trees */
  t8_topidx_t         num_trees_per_eclass[T8_ECLASS_LAST]; /**< Store for each elemet class the number of trees of this class. */
  t8_topidx_t        *tree_to_num_in_eclass; /**< Each tree gets a consecutive index inside the eclass it belongs to. */
  t8_eclass_t        *tree_to_eclass; /**< Store for each tree the element class it belongs to. */
  double             *vertices; /**< An array of size (3* num_vertices). */
  t8_topidx_t        *tree_to_vertex[T8_ECLASS_LAST]; /**< For each eclass a list of (num_trees_per_eclass[class] * num_vertices[class]) indices giving the vertices of each tree. */
}
t8_cmesh_struct_t;

void
t8_cmesh_init (t8_cmesh_t * pcmesh, t8_topidx_t num_vertices,
               t8_topidx_t num_trees,
               t8_topidx_t num_trees_per_eclass[T8_ECLASS_LAST],
               t8_eclass_t * tree_to_eclass)
{
  t8_cmesh_t          cmesh;
  t8_eclass_t         class_it;
  t8_topidx_t         tree_it;
  t8_topidx_t         counter[T8_ECLASS_LAST] = { };

  T8_ASSERT (pcmesh != NULL);
  T8_ASSERT (num_trees > 0);
  T8_ASSERT (tree_to_eclass != NULL);
  cmesh = *pcmesh = T8_ALLOC_ZERO (t8_cmesh_struct_t, 1);
  sc_refcount_init (&cmesh->rc);
  cmesh->num_vertices = num_vertices;
  cmesh->num_trees = num_trees;
  cmesh->tree_to_eclass = T8_ALLOC_ZERO (t8_eclass_t, num_trees);
  memcpy (cmesh->tree_to_eclass, tree_to_eclass,
          sizeof (t8_eclass_t) * num_trees);
  memcpy (cmesh->num_trees_per_eclass, num_trees_per_eclass,
          sizeof (t8_topidx_t) * T8_ECLASS_LAST);
  /* Allocate and fill vertices array */
  if (num_vertices > 0) {
    cmesh->vertices = T8_ALLOC_ZERO (double, 3 * cmesh->num_vertices);
    for (class_it = T8_ECLASS_FIRST; class_it < T8_ECLASS_LAST; class_it++) {
      if (num_trees_per_eclass[class_it] > 0) {
        cmesh->tree_to_vertex[class_it] =
          T8_ALLOC_ZERO (t8_topidx_t,
                         sizeof (t8_topidx_t) *
                         t8_eclass_num_vertices[class_it] *
                         num_trees_per_eclass[class_it]);
      }
    }
  }
  /* Fill tree_to_num_in_eclass */
  cmesh->tree_to_num_in_eclass = T8_ALLOC (t8_topidx_t, num_trees);
  for (tree_it = 0; tree_it < num_trees; tree_it++) {
    cmesh->tree_to_num_in_eclass[tree_it] =
      counter[tree_to_eclass[tree_it]]++;
  }
#ifdef T8_ENABLE_DEBUG
  for (class_it = T8_ECLASS_FIRST; class_it < T8_ECLASS_LAST; class_it++) {
    T8_ASSERT (counter[class_it] == num_trees_per_eclass[class_it]);
  }
#endif
}

static void
t8_cmesh_destroy (t8_cmesh_t * pcmesh)
{
  t8_cmesh_t          cmesh;
  t8_eclass_t         class_it;

  T8_ASSERT (pcmesh != NULL);
  cmesh = *pcmesh;
  T8_ASSERT (cmesh != NULL);
  T8_ASSERT (cmesh->rc.refcount == 0);

  for (class_it = T8_ECLASS_FIRST; class_it < T8_ECLASS_LAST; class_it++) {
    if (cmesh->num_trees_per_eclass[class_it] > 0) {
      T8_FREE (cmesh->tree_to_vertex[class_it]);
    }
  }
  T8_FREE (cmesh->tree_to_num_in_eclass);
  T8_FREE (cmesh->tree_to_eclass);
  T8_FREE (cmesh->vertices);
  T8_FREE (cmesh);

  *pcmesh = NULL;
}

void
t8_cmesh_ref (t8_cmesh_t cmesh)
{
  T8_ASSERT (cmesh != NULL);

  sc_refcount_ref (&cmesh->rc);
}

void
t8_cmesh_unref (t8_cmesh_t * pcmesh)
{
  t8_cmesh_t          cmesh;

  T8_ASSERT (pcmesh != NULL);
  cmesh = *pcmesh;
  T8_ASSERT (cmesh != NULL);

  if (sc_refcount_unref (&cmesh->rc)) {
    t8_cmesh_destroy (pcmesh);
  }
}

t8_cmesh_t
t8_cmesh_new_tet (void)
{
  t8_cmesh_t          cmesh;
  t8_topidx_t         num_trees_per_eclass[T8_ECLASS_LAST] = { };
  t8_eclass_t         tree_to_eclass;

  num_trees_per_eclass[T8_ECLASS_TET] = 1;
  tree_to_eclass = T8_ECLASS_TET;
  t8_cmesh_init (&cmesh, 4, 1, num_trees_per_eclass, &tree_to_eclass);

  return cmesh;
}
