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

#include <t8_shmem.h>

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

void
t8_shmem_array_init (t8_shmem_array_t * parray, size_t elem_size,
                     size_t elem_count, sc_MPI_Comm comm)
{
  t8_shmem_array_t    array;

  T8_ASSERT (parray != NULL);

  array = *parray = T8_ALLOC_ZERO (t8_shmem_array_struct_t, 1);
  sc_shmem_set_type (comm, T8_SHMEM_BEST_TYPE);
  array->array = sc_shmem_malloc (t8_get_package_id (), elem_size, elem_count,
                                  comm);
  array->comm = comm;
  array->elem_count = elem_count;
  array->elem_size = elem_size;
#ifdef T8_ENABLE_DEBUG
  array->shmem_type = T8_SHMEM_BEST_TYPE;
#endif
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

void
t8_shmem_array_destroy (t8_shmem_array_t * parray)
{
  t8_shmem_array_t    array;
  T8_ASSERT (parray != NULL && *parray != NULL);
  array = *parray;
  sc_shmem_free (t8_get_package_id (), array->array, array->comm);
  T8_FREE (array);
  parray = NULL;
}
