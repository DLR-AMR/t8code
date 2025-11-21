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
#include <t8_cmesh/t8_cmesh_internal/t8_cmesh_stash.h>

void
t8_stash_init (t8_stash_t *pstash)
{
  t8_stash_t stash;
  T8_ASSERT (pstash != NULL);

  stash = *pstash = T8_ALLOC (t8_stash_struct_t, 1);
  sc_array_init (&stash->attributes, sizeof (t8_stash_attribute_struct_t));
  sc_array_init (&stash->classes, sizeof (t8_stash_class_struct_t));
  sc_array_init (&stash->joinfaces, sizeof (t8_stash_joinface_struct_t));
}

void
t8_stash_destroy (t8_stash_t *pstash)
{
  t8_stash_t stash;
  t8_stash_attribute_struct_t *attr;
  size_t attr_count;

  T8_ASSERT (pstash != NULL);
  stash = *pstash;
  sc_array_reset (&stash->classes);
  sc_array_reset (&stash->joinfaces);
  for (attr_count = 0; attr_count < stash->attributes.elem_count; attr_count++) {
    attr = (t8_stash_attribute_struct_t *) sc_array_index (&stash->attributes, attr_count);
    if (attr->is_owned) {
      T8_FREE (attr->attr_data);
    }
  }
  sc_array_reset (&stash->attributes);
  T8_FREE (stash);
  pstash = NULL;
}

void
t8_stash_add_class (t8_stash_t stash, const t8_gloidx_t id, const t8_eclass_t eclass)
{
  t8_stash_class_struct_t *sclass;

  T8_ASSERT (stash != NULL);
  sclass = (t8_stash_class_struct_t *) sc_array_push (&stash->classes);
  sclass->eclass = eclass;
  sclass->id = id;
}

/* returns -1 if treeid1 < treeid2
 *          0            =
 *         +1            >
 */
static int
t8_stash_class_compare (const void *c1, const void *c2)
{
  t8_stash_class_struct_t *class1, *class2;

  class1 = (t8_stash_class_struct_t *) c1;
  class2 = (t8_stash_class_struct_t *) c2;
  return class1->id < class2->id ? -1 : class1->id != class2->id;
}

void
t8_stash_class_sort (t8_stash_t stash)
{
  T8_ASSERT (stash != NULL);

  sc_array_sort (&stash->classes, t8_stash_class_compare);
}

static int
t8_stash_class_compare_index (const void *index, const void *c)
{
  t8_gloidx_t index1, index2;

  index1 = *((t8_gloidx_t *) index);
  index2 = ((t8_stash_class_struct_t *) c)->id;

  return index1 < index2 ? -1 : index1 != index2;
}

ssize_t
t8_stash_class_bsearch (const t8_stash_t stash, const t8_gloidx_t tree_id)
{
  return sc_array_bsearch (&stash->classes, &tree_id, t8_stash_class_compare_index);
}

void
t8_stash_add_facejoin (t8_stash_t stash, const t8_gloidx_t gid1, const t8_gloidx_t gid2, const int face1,
                       const int face2, const int orientation)
{
  t8_stash_joinface_struct_t *sjoin;

  T8_ASSERT (stash != NULL);
  sjoin = (t8_stash_joinface_struct_t *) sc_array_push (&stash->joinfaces);
  /* We insert the face connection such that join->id1 is the smaller of
   * the two ids */
  sjoin->face1 = gid1 <= gid2 ? face1 : face2;
  sjoin->face2 = gid1 <= gid2 ? face2 : face1;
  sjoin->id1 = gid1 <= gid2 ? gid1 : gid2;
  sjoin->id2 = gid1 <= gid2 ? gid2 : gid1;
  sjoin->orientation = orientation;
}

/* return -1 if face1.id1 < face2.id1
 *         0              =
 *        +1              >
 */
static int
t8_stash_facejoin_compare (const void *j1, const void *j2)
{
  t8_stash_joinface_struct_t *join1, *join2;

  join1 = (t8_stash_joinface_struct_t *) j1;
  join2 = (t8_stash_joinface_struct_t *) j2;

  return join1->id1 < join2->id1 ? -1 : join1->id1 != join2->id1;
}

void
t8_stash_joinface_sort (const t8_stash_t stash)
{
  T8_ASSERT (stash != NULL);

  sc_array_sort (&stash->joinfaces, t8_stash_facejoin_compare);
}

void
t8_stash_add_attribute (const t8_stash_t stash, const t8_gloidx_t id, const int package_id, const int key,
                        const size_t size, void *const attr, const int copy)
{
  t8_stash_attribute_struct_t *sattr;

  T8_ASSERT (stash != NULL);
  sattr = (t8_stash_attribute_struct_t *) sc_array_push (&stash->attributes);
  sattr->attr_size = size;
  sattr->id = id;
  sattr->is_owned = !copy ? 0 : 1;
  sattr->key = key;
  sattr->package_id = package_id;
  sattr->attr_data = !copy ? attr : T8_ALLOC (char, size);
  if (copy) {
    memcpy (sattr->attr_data, attr, size);
  }
}

size_t
t8_stash_get_attribute_size (const t8_stash_t stash, const size_t index)
{
  return ((t8_stash_attribute_struct_t *) sc_array_index (&stash->attributes, index))->attr_size;
}

void *
t8_stash_get_attribute (const t8_stash_t stash, const size_t index)
{
  return ((t8_stash_attribute_struct_t *) sc_array_index (&stash->attributes, index))->attr_data;
}

t8_gloidx_t
t8_stash_get_attribute_tree_id (const t8_stash_t stash, const size_t index)
{
  return ((t8_stash_attribute_struct_t *) sc_array_index (&stash->attributes, index))->id;
}

int
t8_stash_get_attribute_key (const t8_stash_t stash, const size_t index)
{
  return ((t8_stash_attribute_struct_t *) sc_array_index (&stash->attributes, index))->key;
}

int
t8_stash_get_attribute_id (const t8_stash_t stash, const size_t index)
{
  return ((t8_stash_attribute_struct_t *) sc_array_index (&stash->attributes, index))->package_id;
}

int
t8_stash_attribute_is_owned (const t8_stash_t stash, const size_t index)
{
  return ((t8_stash_attribute_struct_t *) sc_array_index (&stash->attributes, index))->is_owned;
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
t8_stash_attribute_sort (const t8_stash_t stash)
{
  sc_array_sort (&stash->attributes, t8_stash_attribute_compare);
}

static void
t8_stash_bcast_attributes (sc_array_t *attributes, const int root, sc_MPI_Comm comm)
{
  size_t num_atts, iatt, att_size, copied_bytes;
  t8_stash_attribute_struct_t *att;
  char *buffer;
  int mpirank, mpisize, mpiret;
  mpiret = sc_MPI_Comm_rank (comm, &mpirank);
  SC_CHECK_MPI (mpiret);
  mpiret = sc_MPI_Comm_size (comm, &mpisize);
  SC_CHECK_MPI (mpiret);

  T8_ASSERT (attributes != NULL);

  num_atts = attributes->elem_count;

  /* Count the number of byte of all attributes */
  att_size = 0;
  for (iatt = 0; iatt < num_atts; iatt++) {
    att = (t8_stash_attribute_struct_t *) sc_array_index (attributes, iatt);
    att_size += att->attr_size;
  }
  /* Allocate buffer */
  buffer = T8_ALLOC_ZERO (char, att_size);
  /* Copy all attributes to send buffer */
  if (mpirank == root) {
    copied_bytes = 0;
    for (iatt = 0; iatt < num_atts; iatt++) {
      att = (t8_stash_attribute_struct_t *) sc_array_index (attributes, iatt);
      memcpy (buffer + copied_bytes, att->attr_data, att->attr_size);
      copied_bytes += att->attr_size;
    }
  }
  /* broadcast buffer */
  sc_MPI_Bcast (buffer, att_size, sc_MPI_BYTE, root, comm);
  /* Copy attributes from buffer back to stash */
  if (mpirank != root) {
    copied_bytes = 0;
    for (iatt = 0; iatt < num_atts; iatt++) {
      att = (t8_stash_attribute_struct_t *) sc_array_index (attributes, iatt);
      att->attr_data = T8_ALLOC (char, att->attr_size);
      memcpy (att->attr_data, buffer + copied_bytes, att->attr_size);
      copied_bytes += att->attr_size;
    }
  }
  T8_FREE (buffer);
}

/* bcast the data of stash on root to all procs.
 * On the other procs stash_init has to be called before */
t8_stash_t
t8_stash_bcast (const t8_stash_t stash, const int root, sc_MPI_Comm comm, const size_t elem_counts[3])
{
  int mpirank, mpisize, mpiret;
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
    mpiret = sc_MPI_Bcast (stash->attributes.array, elem_counts[0] * sizeof (t8_stash_attribute_struct_t), sc_MPI_BYTE,
                           0, comm);
    SC_CHECK_MPI (mpiret);
    t8_stash_bcast_attributes (&stash->attributes, root, comm);
  }
  if (elem_counts[1] > 0) {
    mpiret
      = sc_MPI_Bcast (stash->classes.array, elem_counts[1] * sizeof (t8_stash_class_struct_t), sc_MPI_BYTE, 0, comm);
    SC_CHECK_MPI (mpiret);
  }
  if (elem_counts[2] > 0) {
    mpiret = sc_MPI_Bcast (stash->joinfaces.array, elem_counts[2] * sizeof (t8_stash_joinface_struct_t), sc_MPI_BYTE, 0,
                           comm);
    SC_CHECK_MPI (mpiret);
  }
  return stash;
}

int
t8_stash_is_equal (const t8_stash_t stash_a, const t8_stash_t stash_b)
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
