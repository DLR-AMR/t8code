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

/* For each tree that we send elements from to
 * other processes, we send this information to the
 * other process */
typedef struct
{
  t8_gloidx_t         gtree_id; /* The global id of that tree *//* TODO: we could optimize this out */
  t8_eclass_t         eclass;   /* The element class of that tree */
  t8_locidx_t         num_elements;     /* The number of elements from this tree that were sent */
} t8_forest_partition_tree_info_t;

/* Given the element offset array and a rank, return the first
 * local element id of this rank */
static              t8_gloidx_t
t8_forest_partition_first_element (t8_gloidx_t * offset, int rank)
{
  return offset[rank];
}

/* Query whether a given process is assigned no elements in
 * an offset array */
static int
t8_forest_partition_empty (t8_gloidx_t * offset, int rank)
{
  if (t8_forest_partition_first_element (offset, rank) >=
      t8_forest_partition_first_element (offset, rank + 1)) {
    return 1;
  }
  return 0;
}

/* Given the element offset array and a rank, return the last
 * local element id of this rank */
static              t8_gloidx_t
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
  int                 i, mpiret, mpisize;

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
  mpiret = sc_MPI_Comm_size (comm, &mpisize);
  SC_CHECK_MPI (mpiret);

  for (i = 0; i < mpisize; i++) {
    /* Calculate the first element index for each process. We convert to doubles to
     * prevent overflow */
    new_first_element_id =
      (((double) i *
        (long double) forest_from->global_num_elements) / (double) mpisize);
    T8_ASSERT (0 <= new_first_element_id &&
               new_first_element_id < forest_from->global_num_elements);
    t8_shmem_array_set_gloidx (forest->element_offsets, i,
                               new_first_element_id);
  }
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

/* Given a tree and its local id, the first and last element id that we need to send to a proc
 * and the first tree we send elements from,
 * calculate the first and last element of this tree that we need to send.
 * The returned element indices are local to the tree.
 */
static void
t8_forest_partition_tree_first_last_el (t8_tree_t tree,
                                        t8_locidx_t tree_id,
                                        t8_locidx_t first_element_send,
                                        t8_locidx_t last_element_send,
                                        t8_locidx_t current_tree,
                                        t8_locidx_t * first_tree_el,
                                        t8_locidx_t * last_tree_el)
{
  if (tree_id == current_tree) {
    /* For the first tree, the first element of that tree that we send is the
     * first element that we send */
    *first_tree_el = first_element_send - tree->elements_offset;
  }
  else {
    *first_tree_el = 0;
  }
  if (tree->elements_offset + tree->elements.elem_count > last_element_send) {
    /* For the last tree, the last element we send is the overall
     * last element that we send */
    *last_tree_el = last_element_send - tree->elements_offset;
  }
  else {
    /* Else it is the last element of this tree */
    *last_tree_el = tree->elements.elem_count - 1;
  }
}

/* Fill the send buffers for one send operation.
 * \param [in]  forest_from     The original forest
 * \param [in]  send_buffer     Unallocated send_buffer
 * \param [in]  current_tree    On input the id of the first tree that we need to send
 *                              elements from. On output the id of the next tree that
 *                              we would send elements from to the next process.
 * \param [in]  first_element_send The local id of the first element that we need to send.
 * \param [in]  last_element_send The local id of the last element that we need to send.
 */
/* The send buffer will look like this:
 *
 * | number of trees | padding | tree_1 info | ... | tree_n info | tree_1 elements | ... | tree_n elements |
 */
static void
t8_forest_partition_fill_buffer (t8_forest_t forest_from,
                                 char **send_buffer, int *buffer_alloc,
                                 t8_locidx_t * current_tree,
                                 t8_locidx_t first_element_send,
                                 t8_locidx_t last_element_send)
{
  t8_locidx_t         num_elements_send;
  t8_tree_t           tree;
  t8_locidx_t         current_element, tree_id, num_trees_send;
  t8_locidx_t         first_tree_element, last_tree_element;
  int                 temp_alloc, byte_alloc, tree_info_pos, element_pos;
  t8_forest_partition_tree_info_t *tree_info;
  t8_locidx_t        *pnum_trees_send;

  current_element = first_element_send;
  tree_id = *current_tree;
  temp_alloc = 0;
  num_trees_send = 0;
  /* At first we calculate the number of bytes that fit in the buffer */
  while (current_element < last_element_send) {
    /* Get the first tree that we send elements from */
    tree = t8_forest_get_tree (forest_from, tree_id);
    t8_forest_partition_tree_first_last_el (tree, tree_id, first_element_send,
                                            last_element_send, *current_tree,
                                            &first_tree_element,
                                            &last_tree_element);
    SC_CHECK_ABORT (tree->eclass == T8_ECLASS_PYRAMID, "Forest partition"
                    " is not implement for pyramidal elements.");
    /* We now know how many elements this tree will send */
    num_elements_send = last_tree_element - first_tree_element + 1;
    T8_ASSERT (num_elements_send > 0);
    temp_alloc += num_elements_send * sizeof (tree->elements.elem_size);
    current_element = last_tree_element + 1;
    num_trees_send++;
  }
  /* We calculate the total number of bytes that we need to allocate
   * and allocate the buffer */
  /* The buffer consists of the number of trees, ... */
  byte_alloc = sizeof (t8_locidx_t);
  /* padding, ... */
  byte_alloc += T8_ADD_PADDING (byte_alloc);
  /* Store the position of the first tree info struct in the buffer */
  tree_info_pos = byte_alloc;
  /* an info struct for each tree, ... */
  byte_alloc = num_trees_send * sizeof (t8_forest_partition_tree_info_t);
  /* Store the position of the first element in the buffer */
  element_pos = byte_alloc;
  /* and the bytes for each tree's elements */
  byte_alloc += temp_alloc;
  /* Note, that we do not add padding after the info structs and
   * each tree's elements, since these are multiples of structs and
   * structs are padded correctly */
  /* We allocate the buffer */
  *send_buffer = T8_ALLOC (char, byte_alloc);
  /* We store the number of trees at first in the send buffer */
  pnum_trees_send = (t8_locidx_t *) send_buffer;
  *pnum_trees_send = num_trees_send;
  for (tree_id = 0; tree_id < num_trees_send; tree_id++) {
    /* Get the first tree that we send elements from */
    tree = t8_forest_get_tree (forest_from, tree_id + *current_tree);
    t8_forest_partition_tree_first_last_el (tree, tree_id + *current_tree,
                                            first_element_send,
                                            last_element_send, *current_tree,
                                            &first_tree_element,
                                            &last_tree_element);
    /* We now know how many elements this tree will send */
    num_elements_send = last_tree_element - first_tree_element + 1;
    /* Get the tree info struct for this tree and fill it */
    tree_info = (t8_forest_partition_tree_info_t *)
      (send_buffer + tree_info_pos);
    tree_info->eclass = tree->eclass;
    tree_info->gtree_id = tree_id + *current_tree +
      forest_from->first_local_tree;
    tree_info->num_elements = num_elements_send;
    tree_info_pos += sizeof (t8_forest_partition_tree_info_t);
    /* We can now fill the send buffer with all elements of that tree */
    memcpy (*send_buffer + element_pos,
            t8_sc_array_index_locidx (&tree->elements, first_tree_element),
            num_elements_send * tree->elements.elem_size);
    element_pos += num_elements_send * tree->elements.elem_size;
  }
  *buffer_alloc = byte_alloc;
}

/* Carry out all sending of elements */
static void
t8_forest_partition_sendloop (t8_forest_t forest, int send_first,
                              int send_last, sc_MPI_Request ** requests,
                              int *num_request_alloc, char ***send_buffer)
{
  int                 iproc, mpiret;
  t8_gloidx_t         gfirst_element_send, glast_element_send;
  t8_gloidx_t         gfirst_local_element, glast_local_element;
  t8_locidx_t         first_element_send, last_element_send;
  t8_locidx_t         current_tree;
  t8_locidx_t         num_elements_send;
  t8_gloidx_t        *offset_to, *offset_from;
  t8_forest_t         forest_from;
  char               *buffer;
  int                 buffer_alloc;
  sc_MPI_Comm         comm;

  T8_ASSERT (t8_forest_is_initialized (forest));
  forest_from = forest->set_from;
  T8_ASSERT (t8_forest_is_committed (forest));

  comm = forest->mpicomm;
  /* Determine the number of requests for MPI communication. */
  /* TODO: Currently we also send to ourselves via MPI. We can optimize this. */
  *num_request_alloc = send_last - send_first;
  *requests = T8_ALLOC (sc_MPI_Request, *num_request_alloc);

  /* Allocate memory for pointers to the send buffers */
  *send_buffer = T8_ALLOC (char *, send_last - send_first + 1);

  /* Get the new and old offset array */
  offset_to = t8_shmem_array_get_gloidx_array (forest->element_offsets);
  offset_from =
    t8_shmem_array_get_gloidx_array (forest_from->element_offsets);

  /* Compute the global ids of the current first and last local element */
  gfirst_local_element = offset_from[forest->mpirank];
  glast_local_element = offset_from[forest->mpirank + 1];
  /* loop over all processes that we send to */
  for (iproc = send_first; iproc < send_last; iproc++) {
    /* At first, we compute the local index of the first and last element
     * that we send to proc */
    if (iproc == send_first) {
      /* If this is the first process we send to, the first element we send is
       * our very first element */
      first_element_send = 0;
      current_tree = 0;
    }
    else {
      /* Otherwise, the first element we send is the new first element on the
       * process */
      gfirst_element_send = offset_to[iproc];
      first_element_send = gfirst_element_send - gfirst_local_element;
      /* assert for overflow error */
      T8_ASSERT ((t8_gloidx_t) first_element_send ==
                 gfirst_element_send - gfirst_local_element);
    }
    if (iproc == send_last) {
      /* To the last process we send all our remaining elements */
      last_element_send = forest_from->local_num_elements - 1;
    }
    else {
      /* Otherwise, the last element we send to proc is the last
       * element on proc in the new partition. */
      glast_element_send = offset_to[iproc + 1] - 1;
      last_element_send = glast_element_send - gfirst_local_element;
    }
    num_elements_send = last_element_send - first_element_send + 1;
    if (num_elements_send < 0) {
      num_elements_send = 0;
    }
    /* We now know the local indices of the first and last element that
     * we send to proc. */
    buffer = *send_buffer[iproc - send_first];
    if (num_elements_send > 0) {
      /* Fill the buffer with the elements on the current tree */
      t8_forest_partition_fill_buffer (forest_from,
                                       &buffer, &buffer_alloc,
                                       &current_tree, first_element_send,
                                       last_element_send);
      /* Post the MPI Send.
       * TODO: This will also send to ourselves if proc==mpirank */
      mpiret = sc_MPI_Isend (buffer, buffer_alloc, sc_MPI_BYTE, iproc,
                             T8_MPI_PARTITION_FOREST, comm,
                             requests[iproc - send_first]);
      SC_CHECK_MPI (mpiret);

    }
    else {
      /* We do not send any elements to iproc (iproc is empty in new partition) */
      /* Set the request to NULL, such that it is ignored when we wait for
       * the requests to complete */
      *(*requests + iproc - send_first) = sc_MPI_REQUEST_NULL;
    }
  }
}

static void
t8_forest_partition_recvloop (t8_forest_t forest, int recv_first,
                              int recv_last)
{
  int                 iproc, num_receive, proc_recv;
  int                *recv_from_data;
  t8_forest_t         forest_from;
  t8_gloidx_t        *offset_from;
  char               *recv_buffer;
  sc_list_t          *receive_from;     /* Linked list storing the ranks from which
                                           we still expect a message */
  sc_link_t          *iterate;  /* Iterator through the linked list */
  sc_link_t          *prev;     /* We need to store iterates predecessor in order to
                                   be able to remove iterate when needed */
  int                 iprobe_flag, mpiret;
  sc_MPI_Comm         comm;
  sc_MPI_Status       status;

  /* Initial checks and inits */
  T8_ASSERT (t8_forest_is_initialized (forest));
  forest_from = forest->set_from;
  T8_ASSERT (t8_forest_is_committed (forest_from));
  offset_from =
    t8_shmem_array_get_gloidx_array (forest_from->element_offsets);
  comm = forest->mpicomm;

  num_receive = recv_last - recv_first + 1;
  /* Initialize the list of processes from which we receive */
  receive_from = sc_list_new (NULL);
  /* The array storing all receiver ranks */
  recv_from_data = T8_ALLOC (int, num_receive);
  for (iproc = recv_first; iproc < recv_last; iproc++) {
    recv_from_data[iproc - recv_first] = iproc;
    if (!t8_forest_partition_empty (offset_from, iproc)) {
      /* If we expect a message from this rank, we store it in our list */
      sc_list_append (receive_from, (void *) (recv_from_data + iproc));
    }
  }
  /****     Actual communication    ****/

  /* Until there is only one sender left we iprobe for a message for each
   * sender and if there is one we receive it and remove the sender from
   * the list.
   * The last message can be received via probe */
  while (receive_from->elem_count > 1) {
    iprobe_flag = 0;
    prev = NULL;
    t8_debugf ("Probing for %zd messages.\n", receive_from->elem_count);
    for (iterate = receive_from->first; iterate != NULL && iprobe_flag == 0;) {
      /* Iterate through all entries and probe for message */
      proc_recv = *(int *) iterate->data;
      mpiret = sc_MPI_Iprobe (proc_recv, T8_MPI_PARTITION_FOREST, comm,
                              &iprobe_flag, &status);
      SC_CHECK_MPI (mpiret);
      if (iprobe_flag == 0) {
        /* We do only continue if there is no message to receive */
        prev = iterate;
        iterate = iterate->next;
      }
    }
    if (iprobe_flag != 0) {
      /* There is a message to receive */
      T8_ASSERT (proc_recv == status.MPI_SOURCE);
      T8_ASSERT (status.MPI_TAG == T8_MPI_PARTITION_FOREST);
      T8_ASSERT (recv_first <= proc_recv && proc_recv <= recv_last);
      /* TODO: Receive elements here */
      sc_list_remove (receive_from, prev);
    }
  }
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
  char              **send_buffer;

  T8_ASSERT (t8_forest_is_initialized (forest));
  T8_ASSERT (forest->set_from != NULL);
  forest_from = forest->set_from;
  T8_ASSERT (t8_forest_is_committed (forest_from));
  /* Compute the first and last rank that we send to */
  t8_forest_partition_sendrange (forest, &send_first, &send_last);
  /* Send all elements to other ranks */
  t8_forest_partition_sendloop (forest, send_first, send_last, &requests,
                                &num_request_alloc, &send_buffer);
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
