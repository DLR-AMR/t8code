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

/** Shared memory array structure.
 * The array uses sc_shmem shared memory.*/
typedef struct t8_shmem_array
{
  void               *array;    /*< Pointer to the actual memory. */
  size_t              elem_size;        /*< Size of one entry in byte. */
  size_t              elem_count;       /*< Total count of entries. */
  sc_MPI_Comm         comm;     /*< MPI communicator. */
  int                 writing_possible; /*< True if we can currently write into this array. False if not. */
  int                 write_start_called;       /*< True if t8_shmem_array_start_writing was called and no call to t8_shmem_array_end_writing happened yet. */
#ifdef T8_ENABLE_DEBUG
  sc_shmem_type_t     shmem_type;       /*< Shared memory type of the communicator (at time of initializing the array). */
#endif
} t8_shmem_array_struct_t;

static int
t8_shmem_array_is_writing_possible (const t8_shmem_array_t array)
{
  return array->writing_possible;
}

#if T8_ENABLE_DEBUG
/* Check whether a shared memory array is initialized. */
int
t8_shmem_array_is_initialized (const t8_shmem_array_t array)
{
  return (array != NULL &&
          array->elem_size > 0 &&
          array->elem_count >= 0 &&
          array->array != NULL && array->comm != sc_MPI_COMM_NULL);
}
#endif

void
t8_shmem_init (sc_MPI_Comm comm)
{
  /* Check whether intranode and internode comms are set
   * for the current communicator. */
  sc_MPI_Comm         intranode;
  sc_MPI_Comm         internode;
  SC_CHECK_ABORT (comm != sc_MPI_COMM_NULL,
                  "Trying to initialize shared memory for NULL communicator.");

  sc_mpi_comm_get_node_comms (comm, &intranode, &internode);
  if (intranode == sc_MPI_COMM_NULL || internode == sc_MPI_COMM_NULL) {
    /* The inter/intra comms are not set. We need to set them to 
     * initialize shared memory usage. */
    sc_mpi_comm_get_and_attach (comm);
  }
}

void
t8_shmem_finalize (sc_MPI_Comm comm)
{
  /* Check whether intranode and internode comms are set
   * for the current communicator. */
  sc_MPI_Comm         intranode;
  sc_MPI_Comm         internode;

  sc_mpi_comm_get_node_comms (comm, &intranode, &internode);
  if (intranode != sc_MPI_COMM_NULL || internode != sc_MPI_COMM_NULL) {
    sc_mpi_comm_detach_node_comms (comm);
  }
}

void
t8_shmem_set_type (sc_MPI_Comm comm, sc_shmem_type_t type)
{
  /* Check whether intranode and internode comms are set
   * for the current communicator. */
  sc_MPI_Comm         intranode;
  sc_MPI_Comm         internode;

  sc_mpi_comm_get_node_comms (comm, &intranode, &internode);
#if T8_ENABLE_MPI
  if (intranode == sc_MPI_COMM_NULL || internode == sc_MPI_COMM_NULL) {
    t8_global_errorf
      ("WARNING: Trying to used shared memory but intranode and internode"
       " communicators are not set."
       " You should call t8_shmem_init before setting the shmem type.\n");
  }
#endif
  sc_shmem_set_type (comm, type);
}

void
t8_shmem_array_init (t8_shmem_array_t *parray, size_t elem_size,
                     size_t elem_count, sc_MPI_Comm comm)
{
  t8_shmem_array_t    array;
  /* Check whether intranode and internode comms are set
   * for the current communicator. */
  sc_MPI_Comm         intranode;
  sc_MPI_Comm         internode;
  SC_CHECK_ABORT (comm != sc_MPI_COMM_NULL,
                  "Trying to initialize shared memory array with NULL communicator.");

  sc_mpi_comm_get_node_comms (comm, &intranode, &internode);
#if T8_ENABLE_MPI
  if (intranode == sc_MPI_COMM_NULL || internode == sc_MPI_COMM_NULL) {
    t8_global_errorf
      ("WARNING: Trying to used shared memory but intranode and internode"
       " communicators are not set."
       " You should call t8_shmem_init before initializing a shared memory array.\n");
  }
#endif

  T8_ASSERT (parray != NULL);

  if (sc_shmem_get_type (comm) == SC_SHMEM_NOT_SET) {
    /* Set the shmem type to the best availble. */
    t8_shmem_set_type (comm, T8_SHMEM_BEST_TYPE);
  }
  array = *parray = T8_ALLOC_ZERO (t8_shmem_array_struct_t, 1);
  array->array = sc_shmem_malloc (t8_get_package_id (), elem_size, elem_count,
                                  comm);
  array->comm = comm;
  array->elem_count = elem_count;
  array->elem_size = elem_size;
  array->writing_possible = 0;
  array->write_start_called = 0;
#ifdef T8_ENABLE_DEBUG
  array->shmem_type = T8_SHMEM_BEST_TYPE;
#endif
}

int
t8_shmem_array_start_writing (t8_shmem_array_t array)
{
  T8_ASSERT (t8_shmem_array_is_initialized (array));

  if (sc_shmem_write_start (array->array, array->comm)) {
    array->writing_possible = 1;
  }
  else {
    array->writing_possible = 0;
  }
  array->write_start_called = 1;
  return array->writing_possible;
}

void
t8_shmem_array_end_writing (t8_shmem_array_t array)
{
  T8_ASSERT (t8_shmem_array_is_initialized (array));

  SC_CHECK_ABORT (array->write_start_called,
                  "End writing to shared array is only possible when t8_shmem_start_writing was called before.");
  sc_shmem_write_end (array->array, array->comm);
  array->write_start_called = 0;
  array->writing_possible = 0;
}

void
t8_shmem_array_copy (t8_shmem_array_t dest, t8_shmem_array_t source)
{
  size_t              bytes;
  T8_ASSERT (t8_shmem_array_is_initialized (dest));
  T8_ASSERT (t8_shmem_array_is_initialized (source));
  T8_ASSERT (!t8_shmem_array_is_writing_possible (dest));
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
t8_shmem_array_allgather (const void *sendbuf, int sendcount,
                          sc_MPI_Datatype sendtype,
                          t8_shmem_array_t recvarray, int recvcount,
                          sc_MPI_Datatype recvtype)
{
  T8_ASSERT (t8_shmem_array_is_initialized (recvarray));
  T8_ASSERT (!t8_shmem_array_is_writing_possible (recvarray));

  sc_shmem_allgather ((void *) sendbuf, sendcount, sendtype, recvarray->array,
                      recvcount, recvtype, recvarray->comm);
}

sc_MPI_Comm
t8_shmem_array_get_comm (t8_shmem_array_t array)
{
  T8_ASSERT (t8_shmem_array_is_initialized (array));
  return array->comm;
}

size_t
t8_shmem_array_get_elem_size (t8_shmem_array_t array)
{
  T8_ASSERT (t8_shmem_array_is_initialized (array));
  return array->elem_size;
}

size_t
t8_shmem_array_get_elem_count (t8_shmem_array_t array)
{
  T8_ASSERT (t8_shmem_array_is_initialized (array));
  return array->elem_count;
}

const t8_gloidx_t  *
t8_shmem_array_get_gloidx_array (t8_shmem_array_t array)
{
  T8_ASSERT (t8_shmem_array_is_initialized (array));
  T8_ASSERT (array->elem_size == sizeof (t8_gloidx_t));
  T8_ASSERT (!t8_shmem_array_is_writing_possible (array));
  return (const t8_gloidx_t *) array->array;
}

t8_gloidx_t        *
t8_shmem_array_get_gloidx_array_for_writing (t8_shmem_array_t array)
{
  T8_ASSERT (t8_shmem_array_is_initialized (array));
  T8_ASSERT (array->elem_size == sizeof (t8_gloidx_t));
  SC_CHECK_ABORT (t8_shmem_array_is_writing_possible (array),
                  "Writing not enabled for shmem array.");
  return (t8_gloidx_t *) array->array;
}

t8_gloidx_t
t8_shmem_array_get_gloidx (t8_shmem_array_t array, int index)
{
  T8_ASSERT (t8_shmem_array_is_initialized (array));
  T8_ASSERT (!t8_shmem_array_is_writing_possible (array));
  T8_ASSERT (array->elem_size == sizeof (t8_gloidx_t));
  T8_ASSERT (0 <= index && (size_t) index < array->elem_count);

  return ((t8_gloidx_t *) array->array)[index];
}

void
t8_shmem_array_set_gloidx (t8_shmem_array_t array, int index,
                           t8_gloidx_t value)
{
  T8_ASSERT (t8_shmem_array_is_initialized (array));
  T8_ASSERT (array->elem_size == sizeof (t8_gloidx_t));
  T8_ASSERT (0 <= index && (size_t) index < array->elem_count);
  SC_CHECK_ABORT (t8_shmem_array_is_writing_possible (array),
                  "Unauthorized write in shmem array. See t8_shmem_array_start_writing.");

  ((t8_gloidx_t *) array->array)[index] = value;
}

const void         *
t8_shmem_array_get_array (t8_shmem_array_t array)
{
  T8_ASSERT (t8_shmem_array_is_initialized (array));
  T8_ASSERT (!t8_shmem_array_is_writing_possible (array));
  return array->array;
}

const void         *
t8_shmem_array_index (t8_shmem_array_t array, size_t index)
{
  T8_ASSERT (t8_shmem_array_is_initialized (array));
  T8_ASSERT (!t8_shmem_array_is_writing_possible (array));
  T8_ASSERT (0 <= index && index < array->elem_count);

  return ((char *) array->array) + index * array->elem_size;
}

void               *
t8_shmem_array_index_for_writing (t8_shmem_array_t array, size_t index)
{
  T8_ASSERT (t8_shmem_array_is_initialized (array));
  T8_ASSERT (t8_shmem_array_is_writing_possible (array));
  T8_ASSERT (0 <= index && index < array->elem_count);

  return ((char *) array->array) + index * array->elem_size;
}

/* TODO: implement */
int
t8_shmem_array_is_equal (t8_shmem_array_t array_a, t8_shmem_array_t array_b)
{
  int                 retval;
  T8_ASSERT (t8_shmem_array_is_initialized (array_a));
  T8_ASSERT (t8_shmem_array_is_initialized (array_b));
  T8_ASSERT (!t8_shmem_array_is_writing_possible (array_a));
  T8_ASSERT (!t8_shmem_array_is_writing_possible (array_b));

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
t8_shmem_array_destroy (t8_shmem_array_t *parray)
{
  t8_shmem_array_t    array;
  T8_ASSERT (parray != NULL);
  T8_ASSERT (t8_shmem_array_is_initialized (*parray));
  T8_ASSERT (!t8_shmem_array_is_writing_possible (*parray));

  array = *parray;
  sc_shmem_free (t8_get_package_id (), array->array, array->comm);
  T8_FREE (array);
  *parray = NULL;
}
