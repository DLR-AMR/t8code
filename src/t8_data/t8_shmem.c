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

/** \file t8_shmem.c
 *
 * TODO: document this file
 */

#include <t8_data/t8_shmem.h>

/* TODO: Think about whether we include a reference counter */

typedef struct t8_shmem_array
{
  void               *array;
  size_t              elem_size;
  size_t              elem_count;
  sc_MPI_Comm         comm;
#ifdef T8_ENABLE_DEBUG
  sc_shmem_type_t     shmem_type;
#endif
} t8_shmem_array_struct_t;

int
t8_shmem_set_type (sc_MPI_Comm comm, sc_shmem_type_t type)
{
  if (sc_shmem_get_type (comm) == SC_SHMEM_NOT_SET) {
    /* This communicator does not have a shmem type, so we set it */
    sc_shmem_set_type (comm, type);
    return 1;
  }
  else {
    /* There was already a type set, we do not change it */
    return 0;
  }
}

void
t8_shmem_array_init (t8_shmem_array_t * parray, size_t elem_size,
                     size_t elem_count, sc_MPI_Comm comm)
{
  t8_shmem_array_t    array;

  T8_ASSERT (parray != NULL);

  array = *parray = T8_ALLOC_ZERO (t8_shmem_array_struct_t, 1);
  array->array = sc_shmem_malloc (t8_get_package_id (), elem_size, elem_count,
                                  comm);
  array->comm = comm;
  array->elem_count = elem_count;
  array->elem_size = elem_size;
#ifdef T8_ENABLE_DEBUG
  array->shmem_type = T8_SHMEM_BEST_TYPE;
#endif
}

void
t8_shmem_array_copy (t8_shmem_array_t dest, t8_shmem_array_t source)
{
  size_t              bytes;
  SC_CHECK_ABORT (t8_shmem_array_get_elem_size (dest) ==
                  t8_shmem_array_get_elem_size (source),
                  "Try to copy shared memory arrays of different element size.\n");
  SC_CHECK_ABORT (t8_shmem_array_get_elem_count (dest) ==
                  t8_shmem_array_get_elem_count (source),
                  "Try to copy shared memory arrays of different element counts.\n");
  SC_CHECK_ABORT (t8_shmem_array_get_comm (dest) ==
                  t8_shmem_array_get_comm (source),
                  "Try to copy shared memory arrays with different communicators.\n");
  /* Get the number of bytes to copy */
  bytes =
    t8_shmem_array_get_elem_count (source) *
    t8_shmem_array_get_elem_size (source);
  sc_shmem_memcpy (dest->array, source->array, bytes, source->comm);
}

void
t8_shmem_array_allgather (void *sendbuf, int sendcount,
                          sc_MPI_Datatype sendtype,
                          t8_shmem_array_t recvarray, int recvcount,
                          sc_MPI_Datatype recvtype)
{
  T8_ASSERT (recvarray != NULL);
  T8_ASSERT (recvarray->array != NULL);

  sc_shmem_allgather (sendbuf, sendcount, sendtype, recvarray->array,
                      recvcount, recvtype, recvarray->comm);
}

sc_MPI_Comm
t8_shmem_array_get_comm (t8_shmem_array_t array)
{
  T8_ASSERT (array != NULL);
  return array->comm;
}

size_t
t8_shmem_array_get_elem_size (t8_shmem_array_t array)
{
  T8_ASSERT (array != NULL);
  return array->elem_size;
}

size_t
t8_shmem_array_get_elem_count (t8_shmem_array_t array)
{
  T8_ASSERT (array != NULL);
  return array->elem_count;
}

t8_gloidx_t        *
t8_shmem_array_get_gloidx_array (t8_shmem_array_t array)
{
  T8_ASSERT (array != NULL);
  T8_ASSERT (array->elem_size == sizeof (t8_gloidx_t));
  return (t8_gloidx_t *) array->array;
}

t8_gloidx_t
t8_shmem_array_get_gloidx (t8_shmem_array_t array, int index)
{
  T8_ASSERT (array != NULL);
  T8_ASSERT (array->array != NULL);
  T8_ASSERT (array->elem_size == sizeof (t8_gloidx_t));
  T8_ASSERT (0 <= index && (size_t) index < array->elem_count);

  return ((t8_gloidx_t *) array->array)[index];
}

void
t8_shmem_array_set_gloidx (t8_shmem_array_t array, int index,
                           t8_gloidx_t value)
{
  T8_ASSERT (array != NULL);
  T8_ASSERT (array->array != NULL);
  T8_ASSERT (array->elem_size == sizeof (t8_gloidx_t));
  T8_ASSERT (0 <= index && (size_t) index < array->elem_count);

  ((t8_gloidx_t *) array->array)[index] = value;
}

void *
t8_shmem_array_get_array (t8_shmem_array_t array)
{
  T8_ASSERT (array != NULL);
  return array->array;
}

void *
t8_shmem_array_index (t8_shmem_array_t array, size_t index)
{
  T8_ASSERT (array != NULL);
  T8_ASSERT (0 <= index && index < array->elem_count);

  return ((char *) array->array) + index * array->elem_size;
}

/* TODO: implement */
int
t8_shmem_array_is_equal (t8_shmem_array_t array_a, t8_shmem_array_t array_b)
{
  int                 retval;

  /* check if direct equality holds */
  if (array_a == array_b) {
    return 1;
  }
  /* check if one or both are NULL */
  if (array_a == NULL || array_b == NULL) {
    return array_a == array_b;
  }
  /* Compare metadata */
  retval = array_a->comm != array_b->comm ||
    array_a->elem_count != array_b->elem_count
    || array_a->elem_size != array_b->elem_size;
  if (retval != 0) {
    return 0;
  }
  /* Check the contents of the arrays */
  retval = memcmp (array_a->array, array_b->array,
                   array_a->elem_count * array_a->elem_size);
  if (retval == 0) {
    /* If memcmp returned 0, all checks were successful */
    return 1;
  }
  return 0;
}

void
t8_shmem_array_destroy (t8_shmem_array_t * parray)
{
  t8_shmem_array_t    array;
  T8_ASSERT (parray != NULL && *parray != NULL);
  array = *parray;
  sc_shmem_free (t8_get_package_id (), array->array, array->comm);
  T8_FREE (array);
  *parray = NULL;
}
