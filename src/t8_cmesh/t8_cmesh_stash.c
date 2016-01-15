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

/** \file t8_cmesh_stash.c
 * We define the data structures and routines for temporary storage before commit
 */

#include <t8.h>
#include <t8_eclass.h>
#include "t8_cmesh_stash.h"

void
t8_stash_init (t8_stash_t * pstash)
{
  t8_stash_t          stash;
  T8_ASSERT (pstash != NULL);

  stash = *pstash = T8_ALLOC (t8_stash_struct_t, 1);
  sc_array_init (&stash->attributes, sizeof (t8_stash_attribute_struct_t));
  sc_array_init (&stash->classes, sizeof (t8_stash_class_struct_t));
  sc_array_init (&stash->joinfaces, sizeof (t8_stash_joinface_struct_t));
}

void
t8_stash_destroy (t8_stash_t * pstash)
{
  t8_stash_t          stash;
  t8_stash_attribute_struct_t *attr;
  size_t              attr_count;

  T8_ASSERT (pstash != NULL);
  stash = *pstash;
  sc_array_reset (&stash->classes);
  sc_array_reset (&stash->joinfaces);
  for (attr_count = 0; attr_count < stash->attributes.elem_count;
       attr_count++) {
    attr =
      (t8_stash_attribute_struct_t *) sc_array_index (&stash->attributes,
                                                      attr_count);
    if (attr->is_owned) {
      T8_FREE (attr->attr_data);
    }
  }
  sc_array_reset (&stash->attributes);
  T8_FREE (stash);
  pstash = NULL;
}

void
t8_stash_add_class (t8_stash_t stash, t8_gloidx_t id, t8_eclass_t eclass)
{
  t8_stash_class_struct_t *sclass;

  T8_ASSERT (stash != NULL);
  sclass = (t8_stash_class_struct_t *) sc_array_push (&stash->classes);
  sclass->eclass = eclass;
  sclass->id = id;
}

void
t8_stash_add_facejoin (t8_stash_t stash, t8_gloidx_t id1, t8_gloidx_t id2,
                       int face1, int face2, int orientation)
{
  t8_stash_joinface_struct_t *sjoin;

  T8_ASSERT (stash != NULL);
  sjoin = (t8_stash_joinface_struct_t *) sc_array_push (&stash->joinfaces);
  /* We inserte the face connection such that join->id1 is the smaller of
   * the two ids */
  sjoin->face1 = id1 <= id2 ? face1 : face2;
  sjoin->face2 = id1 <= id2 ? face2 : face1;;
  sjoin->id1 = id1 <= id2 ? id1 : id2;;
  sjoin->id2 = id1 <= id2 ? id2 : id1;;
  sjoin->orientation = orientation;
}

void
t8_stash_add_attribute (t8_stash_t stash, t8_gloidx_t id, int package_id,
                        int key, size_t size, void *attr, int data_persists)
{
  t8_stash_attribute_struct_t *sattr;

  T8_ASSERT (stash != NULL);
  sattr = (t8_stash_attribute_struct_t *) sc_array_push (&stash->attributes);
  sattr->attr_size = size;
  sattr->id = id;
  sattr->is_owned = data_persists ? 0 : 1;
  sattr->key = key;
  sattr->package_id = package_id;
  sattr->attr_data = data_persists ? attr : T8_ALLOC (char, size);
  if (!data_persists) {
    memcpy (sattr->attr_data, attr, size);
  }
}

size_t
t8_stash_get_attribute_size (t8_stash_t stash, size_t index)
{
  return ((t8_stash_attribute_struct_t *)
          sc_array_index (&stash->attributes, index))->attr_size;
}

void               *
t8_stash_get_attribute (t8_stash_t stash, size_t index)
{
  return ((t8_stash_attribute_struct_t *)
          sc_array_index (&stash->attributes, index))->attr_data;
}

t8_gloidx_t
t8_stash_get_attribute_tree_id (t8_stash_t stash, size_t index)
{
  return ((t8_stash_attribute_struct_t *)
          sc_array_index (&stash->attributes, index))->id;
}

int
t8_stash_get_attribute_key (t8_stash_t stash, size_t index)
{
  return ((t8_stash_attribute_struct_t *)
          sc_array_index (&stash->attributes, index))->key;
}

int
t8_stash_get_attribute_id (t8_stash_t stash, size_t index)
{
  return ((t8_stash_attribute_struct_t *)
          sc_array_index (&stash->attributes, index))->package_id;
}

int
t8_stash_attribute_is_owned (t8_stash_t stash, size_t index)
{
  return ((t8_stash_attribute_struct_t *)
          sc_array_index (&stash->attributes, index))->is_owned;
}

/* Compare two attribute entries A1 and A2.
 * A1 is smaller than A2 if and only if its treeid is smaller or (if equal)
 * its package id is smaller or (if also equal) its key is smaller.
 */
static int
t8_stash_attribute_compare (const void *v1, const void *v2)
{
  t8_stash_attribute_struct_t *A1 = (t8_stash_attribute_struct_t *) v1;
  t8_stash_attribute_struct_t *A2 = (t8_stash_attribute_struct_t *) v2;

  if (A1->id == A2->id) {
    if (A1->package_id == A2->package_id) {
      return A1->key < A2->key ? -1 : A1->key > A2->key;
    }
    return A1->package_id < A2->package_id ? -1 : 1;
  }
  return A1->id < A2->id ? -1 : A1->id > A2->id;
}

/* Sort the attribute entries in the order
 * (treeid, packageid, key)
 */
void
t8_stash_attribute_sort (t8_stash_t stash)
{
  sc_array_sort (&stash->attributes, t8_stash_attribute_compare);
}

/* bcast the data of stash on root to all procs.
 * On the other procs stash_init has to be called before */
t8_stash_t
t8_stash_bcast (t8_stash_t stash, int root, sc_MPI_Comm comm,
                size_t elem_counts[3])
{
  int                 mpirank, mpisize, mpiret;
  mpiret = sc_MPI_Comm_rank (comm, &mpirank);
  SC_CHECK_MPI (mpiret);
  mpiret = sc_MPI_Comm_size (comm, &mpisize);
  SC_CHECK_MPI (mpiret);

  if (mpirank != root) {
    sc_array_resize (&stash->attributes, elem_counts[0]);
    sc_array_resize (&stash->classes, elem_counts[1]);
    sc_array_resize (&stash->joinfaces, elem_counts[2]);
  }
  if (elem_counts[0] > 0) {
    mpiret = sc_MPI_Bcast (stash->attributes.array,
                           elem_counts[0] *
                           sizeof (t8_stash_attribute_struct_t), sc_MPI_BYTE,
                           0, comm);
    SC_CHECK_MPI (mpiret);
  }
  if (elem_counts[1] > 0) {
    mpiret = sc_MPI_Bcast (stash->classes.array,
                           elem_counts[1] * sizeof (t8_stash_class_struct_t),
                           sc_MPI_BYTE, 0, comm);
    SC_CHECK_MPI (mpiret);
  }
  if (elem_counts[2] > 0) {
    mpiret = sc_MPI_Bcast (stash->joinfaces.array,
                           elem_counts[2] *
                           sizeof (t8_stash_joinface_struct_t), sc_MPI_BYTE,
                           0, comm);
    SC_CHECK_MPI (mpiret);
  }
  return stash;
}

int
t8_stash_is_equal (t8_stash_t stash_a, t8_stash_t stash_b)
{
  if (stash_a == stash_b) {
    return 1;
  }
  if (stash_a == NULL || stash_b == NULL) {
    return 0;
  }
  return (sc_array_is_equal (&stash_a->attributes, &stash_b->attributes)
          && sc_array_is_equal (&stash_a->classes, &stash_b->classes)
          && sc_array_is_equal (&stash_a->joinfaces, &stash_b->joinfaces));
}
