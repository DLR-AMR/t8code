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

#include <t8_forest/t8_forest_partition.h>
#include <t8_forest/t8_forest_types.h>
#include <t8_cmesh/t8_cmesh_offset.h>
#include <t8_forest.h>

/* Given the element offset array and a rank, return the first
 * local element id of this rank */
static t8_gloidx_t
t8_forest_partition_first_element (t8_gloidx_t *offset, int rank)
{
  return offset[rank];
}

/* Given the element offset array and a rank, return the last
 * local element id of this rank */
static t8_gloidx_t
t8_forest_partition_last_element (t8_gloidx_t * offset, int rank)
{
  return offset[rank + 1] - 1;
}

/* For a committed forest create the array of element_offsets
 * and store it in forest->element_offsets
 */
static void
t8_forest_partition_create_offsets (t8_forest_t forest)
{
  sc_MPI_Comm         comm;
  t8_gloidx_t         first_local_element;

  T8_ASSERT (t8_forest_is_committed (forest));

  T8_ASSERT (forest->element_offsets == NULL);
  comm = forest->mpicomm;
  /* Set the shmem array type to comm */
  sc_shmem_set_type (comm, T8_SHMEM_BEST_TYPE);
  /* Initialize the offset array as a shmem array
   * holding mpisize+1 many t8_gloidx_t */
  t8_shmem_array_init (&forest->element_offsets, sizeof (t8_gloidx_t),
                       forest->mpisize + 1, comm);
  /* Calculate the global index of the first local element */
  t8_forest_get_first_local_element_id (forest);
  /* Collect all first global indices in the array */
  t8_shmem_array_allgather (&first_local_element, 1, T8_MPI_GLOIDX,
                            forest->element_offsets, 1, T8_MPI_GLOIDX);
}

/* Calculate the new element_offset for forest from
 * the element in forest->set_from assuming a partition without
 * element weights */
static void
t8_forest_partition_compute_new_offset (t8_forest_t forest)
{
  t8_forest_t         forest_from;
  sc_MPI_Comm         comm;
  t8_gloidx_t         new_first_element_id;

  T8_ASSERT (t8_forest_is_initialized (forest));
  T8_ASSERT (forest->set_from != NULL);

  forest_from = forest->set_from;
  comm = forest->mpicomm;

  T8_ASSERT (forest->element_offsets == NULL);
  /* Set the shmem array type to comm */
  sc_shmem_set_type (comm, T8_SHMEM_BEST_TYPE);
  /* Initialize the shmem array */
  t8_shmem_array_init (&forest->element_offsets, sizeof (t8_gloidx_t),
                       forest->mpisize + 1, comm);
  /* Calculate the first element index. We convert to doubles to
   * prevent overflow */
  new_first_element_id =
    (((double) forest->mpirank *
      (long double) forest_from->global_num_elements) /
     (double) forest->mpisize);
  T8_ASSERT (0 <= new_first_element_id &&
             new_first_element_id < forest_from->global_num_elements);
  /* Fill the offset array with the first element id of each process */
  t8_shmem_array_allgather (&new_first_element_id, 1, T8_MPI_GLOIDX,
                            forest->element_offsets, 1, T8_MPI_GLOIDX);
}

/* Find the owner of a given element.
 */
static int
t8_forest_partition_owner_of_element (int mpisize, t8_gloidx_t gelement,
                                      t8_gloidx_t * offset)
{
  /* Tree offsets are stored similar enough that we can exploit their function */
  /* In the element offset logic, an element cannot be owned by more than one
   * process, thus any owner must be the unique owner. */
  return t8_offset_any_owner_of_tree (mpisize, gelement, offset);
}

/* Compute the first and last rank that we need to receive elements from */
static void
t8_forest_partition_recvrange (t8_forest_t forest, int *recv_first,
                               int *recv_last)
{
  t8_gloidx_t         first_element, last_element;
  t8_gloidx_t        *offset_old, *offset_new;

  /* Get the old element offset array */
  offset_old =
    t8_shmem_array_get_gloidx_array (forest->set_from->element_offsets);
  /* Get the new element offset array */
  offset_new = t8_shmem_array_get_gloidx_array (forest->element_offsets);
  /* Compute new first and last element on this process from offset array */
  first_element =
    t8_forest_partition_first_element (offset_new, forest->mpisize);
  last_element =
    t8_forest_partition_last_element (offset_new, forest->mpisize);
  /* Calculate the first and last process we receive from */
  *recv_first = t8_forest_partition_owner_of_element (forest->mpisize,
                                                      first_element,
                                                      offset_old);
  *recv_last = t8_forest_partition_owner_of_element (forest->mpisize,
                                                     last_element,
                                                     offset_old);
}

/* Compute the first and last rank that we need to send elements to */
static void
t8_forest_partition_sendrange (t8_forest_t forest, int *send_first,
                               int *send_last)
{
  t8_gloidx_t         first_element, last_element;
  t8_gloidx_t        *offset_old, *offset_new;

  /* Get the old element offset array */
  offset_old =
    t8_shmem_array_get_gloidx_array (forest->set_from->element_offsets);
  /* Get the new element offset array */
  offset_new = t8_shmem_array_get_gloidx_array (forest->element_offsets);
  /* Compute old first and last element on this process from offset array */
  first_element =
    t8_forest_partition_first_element (offset_old, forest->mpisize);
  last_element =
    t8_forest_partition_last_element (offset_old, forest->mpisize);
  /* Calculate the first and last process we send to */
  *send_first = t8_forest_partition_owner_of_element (forest->mpisize,
                                                      first_element,
                                                      offset_new);
  *send_last = t8_forest_partition_owner_of_element (forest->mpisize,
                                                     last_element,
                                                     offset_new);
}

/* Carry out all sending of elements */
static void
t8_forest_partition_sendloop (t8_forest_t forest, int send_first,
                              int send_last, sc_MPI_Request ** requests,
                              int *num_request_alloc)
{
  int                 flag;

  /* Determine the number of requests for MPI communication.
   * If the current process is included in the sendrange we do not allocate
   * a request and buffer for it. */
  flag = forest->mpirank < send_first || forest->mpirank > send_last ? 1 : 0;
  *num_request_alloc = send_last - send_first + flag;
  *requests = T8_ALLOC (sc_MPI_Request, *num_request_alloc);
}

static void
t8_forest_partition_recvloop (t8_forest_t forest, int recv_first,
                              int recv_last)
{
}

/* Partition a forest from forest->set_from and the element offsets
 * set in forest->element_offsets
 */
static void
t8_forest_partition_given (t8_forest_t forest)
{
  t8_forest_t         forest_from;
  int                 send_first, send_last, recv_first, recv_last;
  sc_MPI_Request     *requests = NULL;
  int                 num_request_alloc;        /* The count of elemtns in the request array */

  T8_ASSERT (t8_forest_is_initialized (forest));
  T8_ASSERT (forest->set_from != NULL);
  forest_from = forest->set_from;
  T8_ASSERT (t8_forest_is_committed (forest_from));
  /* Compute the first and last rank that we send to */
  t8_forest_partition_sendrange (forest, &send_first, &send_last);
  /* Send all elements to other ranks */
  t8_forest_partition_sendloop (forest, send_first, send_last, &requests,
                                &num_request_alloc);
  /* Receive all element from other ranks */
  t8_forest_partition_recvrange (forest, &recv_first, &recv_last);
  t8_forest_partition_recvloop (forest, recv_first, recv_last);
}

/* Populate a forest with the partitioned elements of
 * forest->set_from.
 * Currently the elements are distributed evenly (each element
 * has the same weight).
 */
void
t8_forest_partition (t8_forest_t forest)
{
  t8_forest_t         forest_from;

  T8_ASSERT (!t8_forest_is_committed (forest));
  forest_from = forest->set_from;
  T8_ASSERT (t8_forest_is_committed (forest_from));

  if (forest_from->element_offsets == NULL) {
    /* We create the partition table of forest_from */
    t8_forest_partition_create_offsets (forest_from);
  }
  /* TODO: if offsets already exist on forest_from, check it for consistency */

  /* We now calculate the new element offsets */
  t8_forest_partition_compute_new_offset (forest);
  t8_forest_partition_given (forest);
}
