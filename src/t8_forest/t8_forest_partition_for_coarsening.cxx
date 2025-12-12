/*
  This file is part of t8code.
  t8code is a C library to manage a collection (a forest) of multiple
  connected adaptive space-trees of general element classes in parallel.

  Copyright (C) 2025 the developers

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

/**
 * \file This file contains the main implementations of t8code's partition-for-coarsening feature.
*/

#include <t8_forest/t8_forest.h>
#include <t8_forest/t8_forest_types.h>
#include <t8_forest/t8_forest_private.h>
#include <t8_forest/t8_forest_pfc_message.hxx>
#include <t8_forest/t8_forest_pfc_helper.hxx>
#include <t8_data/t8_shmem.h>
#include <t8_element.h>
#include <t8_forest/t8_forest_partition_for_coarsening.h>
#include <t8_eclass.h>
#include <vector>
#include <algorithm>

/** Return the process owner of the given (global) element_id.
 *
 * \param[in]   partition   the current partitioning, given as array of element offsets.
 * \param[in]   mpisize     the number of MPI ranks
 * \param[in]   element_id  the global index of the element considered
 *
 * \return The ID of the process owning the given element.
 **/
constexpr t8_procidx_t
proc_owner (const t8_gloidx_t *partition, const t8_procidx_t mpisize, const t8_gloidx_t element_id)
{
  // Due to the pointer arithmetics, this may look more complicated than it is:
  // ----
  // (1.) The upper_bound returns a pointer to the first entry in partition that is bigger than element_id.
  // (2.) Subtracting 1 gives the pointer to the previous entry in partition, i.e., the entry corresponding
  //      to the process holding the considered element.
  // (3.) Subtracting the pointer partition gives the number of entries that would fit in between the two
  //      memory addresses, i.e., it translates the entry's memory address into its array index, i.e., the
  //      the process id.
  return (std::upper_bound (partition, partition + mpisize, element_id) - 1) - partition;
}

/**
 * Determine the range of processes that may potentially hold siblings of any process-local element.
 *
 * To do so, the range of relevant elements is computed first, essentially by extending the local
 * element range in both directions by max_num_siblings minus one.
 * With that, the corrdesponding processes are simple to find.
 *
 * \param[in]   partition_old     partition of the source forest (as C-style array)
 * \param[in]   mpirank           MPI rank of the current process
 * \param[in]   mpisize           MPI size of the current communicator
 * \param[out]  proc_range_begin  Id of the process holding the first relevant element.
 * \param[out]  proc_range_end    Id of the process after the last one holding a relevant element.
*/
void
get_relevant_process_range (const t8_gloidx_t *partition_old, const t8_procidx_t rank, const t8_procidx_t mpisize,
                            t8_procidx_t &proc_range_begin, t8_procidx_t &proc_range_end)
{

  // Get maximum number of siblings.
  const int max_num_siblings = T8_ECLASS_MAX_CHILDREN;

  // Get range of elements that are relevant to find all families of the process-local elements.
  // We have to consider all elements that may potentially be siblings of local elements.
  //    (Note: The SC_MAX and SC_MIN commands are only relevant for the first and last process.)
  const t8_gloidx_t relevant_eles_begin = SC_MAX (0, partition_old[rank] - (max_num_siblings - 1));
  const t8_gloidx_t relevant_eles_end
    = SC_MIN (partition_old[mpisize], partition_old[rank + 1] + (max_num_siblings - 1));

  // Get the process range holding these elements.
  proc_range_begin = proc_owner (partition_old, mpisize, relevant_eles_begin);
  proc_range_end = proc_owner (partition_old, mpisize, relevant_eles_end - 1) + 1;

  // Sanity checks
  T8_ASSERT (0 <= proc_range_begin);
  T8_ASSERT (proc_range_begin <= proc_range_end);
  T8_ASSERT (proc_range_end <= mpisize);
}

/**
 * Send PFC messages to all relevant processes and obtain the associated requests.
 *
 * It is important to understand that we compute the PFC corrections before the new partition
 * is applied. Consequently, all PFC procedures are based on the old partition of the source
 * forest \a forest_from, so the partition offsets we correct in general do not match
 * those of \a forest_from.
 * This in particular means that it does not suffice to only receive from one direction, e.g.,
 * from all processes of lower rank and then only check the lower process borders, see comment
 * in description of \ref t8_forest_pfc_family_range_around_border.
 *
 * \param[in]   forest_from    the old forest (forest->set_from of the new partitioned forest)
 * \param[out]  requests  the MPI requests as std::vector
 *                          on input:  empty
 *                          on output: contains the send requests
*/
template <typename MessageType>
static void
t8_forest_pfc_send_loop_range (const t8_forest_t forest_from, std::vector<sc_MPI_Request> &requests)
{
  // Assertions: The forest must be committed and the request vector empty.
  T8_ASSERT (t8_forest_is_committed (forest_from));
  T8_ASSERT (requests.size () == 0);

  // Get current offset vector of forest.
  const t8_gloidx_t *partition_old = t8_shmem_array_get_gloidx_array (forest_from->element_offsets);

  // Get range of elements that may potentially be a sibling of a process-local element.
  t8_procidx_t begin;
  t8_procidx_t end;
  get_relevant_process_range (partition_old, forest_from->mpirank, forest_from->mpisize, begin, end);

  // Loop over processes of relevant process range to send data if required.
  for (t8_procidx_t iproc = begin; iproc < end; iproc++) {

    // Skip empty processes and the own rank.
    if (partition_old[iproc] >= partition_old[iproc + 1] || iproc == forest_from->mpirank)
      continue;

    // Construct and fill message of type MessageType (see t8_forest_pfc_message).
    MessageType message (forest_from->scheme, iproc, forest_from->mpicomm);
    message.fill (forest_from);

    // Send the message (and obtain the associated requests).
    sc_MPI_Request request;
    message.mpi_Isend (forest_from, request);

    // Add to requests array
    requests.push_back (std::move (request));
  }
}

/**
 *  Determine the messages to be received.
 *
 * It is important to understand that we compute the PFC corrections before the new partition
 * is applied. Consequently, all PFC procedures are based on the old partition of the source
 * forest \a forest_from, so the partition offsets we correct in general do not match
 * those of \a forest_from.
 * This in particular means that it does not suffice to only receive from one direction, e.g.,
 * from all processes of lower rank and then only check the lower process borders, see comment
 * in description of \ref t8_forest_pfc_family_range_around_border.
 *
 * \param[in]   forest_from   the old forest (forest->set_from of the new partitioned forest)
 * \param[out]  requests      the MPI messages
*/
template <typename MessageType>
static void
t8_forest_pfc_recv_loop_range (const t8_forest_t forest_from, std::vector<MessageType> &messages)
{
  // Assertions: forest must be committed and messages empty.
  T8_ASSERT (t8_forest_is_committed (forest_from));
  T8_ASSERT (messages.size () == 0);

  // Get current offset vector of forest.
  const t8_gloidx_t *partition_old = t8_shmem_array_get_gloidx_array (forest_from->element_offsets);

  // Get range of elements that may potentially be a sibling of a process-local element.
  t8_procidx_t begin;
  t8_procidx_t end;
  get_relevant_process_range (partition_old, forest_from->mpirank, forest_from->mpisize, begin, end);

  // Loop over process range to receive messages.
  for (t8_procidx_t iproc = begin; iproc < end; iproc++) {

    // Skip empty partitions and the own rank.
    if (partition_old[iproc] >= partition_old[iproc + 1] || iproc == forest_from->mpirank)
      continue;

    // Receive message.
    MessageType message (forest_from->scheme, iproc, forest_from->mpicomm);
    t8_debugf ("receive message from %i\n", message.iproc);
    int buf_size;
    char *recv_buf;
    message.mpi_Recv (recv_buf, buf_size);

    // Unpack message
    int position = 0;
    message.unpack (recv_buf, buf_size, &position);
    T8_ASSERT (position == buf_size);
    T8_FREE (recv_buf);

    // Push to messages
    messages.push_back (std::move (message));
  }
}

/**
 * Determine whether a full family would be split by a process boundary in the new partitioning.
 *
 * It is important to note that we run this check for the new partition which is not applied yet.
 * Therefore, this function is based on the old partition of the source forest \a forest_from,
 * meaning the \a border_element_id we check will in general not be at a process boundary of
 * \a forest_from.
 * The \a messages received from other processes contain all elements that may potentially form
 * a family with those held by this process p in \a forest_from 's partitioning.
 * With this set of elements, we can identify all siblings of \a border_element_id and check whether
 * they form a full family. Note that even if they do not, we pass back the range of sibgling elements
 * via the output arguments \a family_begin and \a family_end.
 *
 * Note that we can see here why the two-directional communication, i.e., sending to / receiving
 * from processes with both higher and lower rank, is the simplest way to over cover all cases up to the
 * extreme ones, i.e., \a border_element_id being
 *    (a) the first local element of \a forest_from and the last sibling of a full family, or
 *    (b) the last local element of \a forest_from and the first sibling of a full family.
 * A one-directional communication pattern would be possible, but would require to double the message sizes
 * and make this function less convenient as we would have to check \a boder_element_id s not part
 * of our element range in \b forest_from.
 *
 * \param[in]   forest_from         the old forest (forest->set_from of the new partitioned forest)
 * \param[in]   border_element_id   the global ID of the new partitioning's border element we run the check for
 * \param[in]   messages            the PFC messages received from other processes
 * \param[out]  family_begin        the global element ID of the family's first member
 * \param[out]  family_end          the global element ID of the family's last member
 *
 * \return True (i.e., nonzero) if and only if a full family is found across the process borders.
*/
static int
t8_forest_pfc_family_range_around_border (const t8_forest_t forest_from, const t8_gloidx_t border_element_id,
                                          const std::vector<t8_forest_pfc_message> &messages, t8_gloidx_t &family_begin,
                                          t8_gloidx_t &family_end)
{
  // From the element with global ID border_element_id, determine
  // - the global tree ID
  // - the tree
  // - the element's index within the tree
  // - the element itself
  t8_gloidx_t gtree_id;
  t8_tree_t tree;
  t8_locidx_t index_in_tree;
  t8_element_t *element;
  t8_forest_pfc_helper_index_in_tree_from_globalid (forest_from, border_element_id, gtree_id, tree, index_in_tree,
                                                    element);

  // Get scheme and eclass from forest and tree.
  const t8_scheme_c *scheme = t8_forest_get_scheme (forest_from);
  const t8_eclass_t eclass = tree->eclass;

  // If the element is the root, return false because the root does not have any parent or siblings.
  if (scheme->element_get_level (eclass, element) == 0) {
    family_begin = border_element_id;
    family_end = border_element_id;
    return false;
  }

  // Allocate and determine parent element.
  t8_element_t *parent;
  t8_element_new (scheme, eclass, 1, &parent);
  scheme->element_get_parent (eclass, element, parent);

  // Get global ID of first (process-)local element
  const t8_gloidx_t first_tree_element
    = t8_forest_get_first_local_leaf_element_id (forest_from) + tree->elements_offset;

  // Determine range of global IDs forming the family of first_tree_element, by calling the helper function
  // t8_forest_pfc_extreme_local_sibling twice, i.e., searching in the direction of in- and decreasing indices.
  // Note: The end iterator is one behind the last family member.
  family_begin = first_tree_element + t8_forest_pfc_extreme_local_sibling (scheme, tree, index_in_tree, -1);
  family_end = first_tree_element + t8_forest_pfc_extreme_local_sibling (scheme, tree, index_in_tree, +1) + 1;

  // Check if other processes have the same parent as the current family, so we need to adjust our range.
  for (t8_procidx_t imessage = 0; imessage < (t8_procidx_t) messages.size (); imessage++) {
    t8_debugf ("process message from %i\n", messages[imessage].iproc);

    // On the same tree we can use our scheme to compare, because we know that the eclasses are equal.
    if (messages[imessage].itree == gtree_id
        && scheme->element_is_equal (eclass, parent, messages[imessage].get_parent ())) {

      // If parents are equal, extend lower or upper range border by num_siblings, depending on the send "direction",
      // i.e., towards lower- or higher-rank processes.
      if (messages[imessage].iproc < forest_from->mpirank) {
        family_begin -= messages[imessage].num_siblings;
      }
      else {
        family_end += messages[imessage].num_siblings;
      }
    }
  }

  // Determine the parent's number of children.
  const int num_children = scheme->element_get_num_children (eclass, parent);

  // Deallocate parent element
  t8_element_destroy (scheme, eclass, 1, &parent);

  // Return true if the considered family contains all children of the parent.
  return (family_end - family_begin == num_children);
}

/** Compute the process-local corrections of the given partition.
 *
 * \param[in]   forest_from             the old forest (forest->set_from of the new partitioned forest)
 * \param[in]   partition_new_shmem     the current partitioning (without PFC correction) as shared-memory array
 * \param[in]   messages                the PFC messages received from other processes
 * \param[out]  corrected_local_offsets a std::vector of t8_gloidx_t>
 *                                      on input:  empty
 *                                      on output: containing the corrections to be applied to the local offsets to obtain the PFC partitioning
*/
static void
t8_forest_pfc_correct_local_offsets (const t8_forest_t forest_from, const t8_shmem_array_t partition_new_shmem,
                                     const std::vector<t8_forest_pfc_message> &messages,
                                     std::vector<t8_gloidx_t> &corrected_local_offsets)
{
  T8_ASSERT (t8_forest_is_committed (forest_from));

  // Get current partitioning as array of t8_gloidx_t.
  const t8_gloidx_t *partition_new = t8_shmem_array_get_gloidx_array (partition_new_shmem);

  // Get offsets of current and next process in old partitioning.
  // (Note that here the next process is used solely to know the upper bounds.)
  const t8_gloidx_t old_offset_this_process
    = t8_shmem_array_get_gloidx (forest_from->element_offsets, forest_from->mpirank);
  const t8_gloidx_t old_offset_next_process
    = t8_shmem_array_get_gloidx (forest_from->element_offsets, forest_from->mpirank + 1);

  // Determine where these offsets would end up in the new partitioning.
  const t8_gloidx_t *min_local_element_pointer
    = std::lower_bound (partition_new, partition_new + forest_from->mpisize, old_offset_this_process);
  const t8_gloidx_t *next_min_local_element_pointer
    = std::lower_bound (partition_new, partition_new + forest_from->mpisize, old_offset_next_process);

  // Get the processes that will (without correction) get the elements the current process had in the old partitioning.
  const t8_gloidx_t min_local_proc = min_local_element_pointer - partition_new;
  const t8_gloidx_t next_min_local_proc = next_min_local_element_pointer - partition_new;

  // Loop over this range of processes to adjust local borders.
  for (t8_procidx_t border_rank = min_local_proc; border_rank < next_min_local_proc; border_rank++) {
    t8_gloidx_t family_begin, family_end;

    // Check if there is a full family split by the current border. If so, adjust border.
    if (t8_forest_pfc_family_range_around_border (forest_from, partition_new[border_rank], messages, family_begin,
                                                  family_end)) {
      // Find process owning first family member.
      const t8_procidx_t rank_family_begin = proc_owner (partition_new, forest_from->mpisize, family_begin);

      // Push corrected local offset to vector: Depending on the split rank, to beginning or end of family.
      const t8_gloidx_t new_offset = (border_rank <= rank_family_begin) ? family_begin : family_end;
      corrected_local_offsets.push_back (new_offset);
    }
    else {
      // No correction needed: Push current offset to vector.
      const t8_gloidx_t new_offset = partition_new[border_rank];
      corrected_local_offsets.push_back (new_offset);
    }
  }
}

T8_EXTERN_C_BEGIN ();
void
t8_forest_pfc_correction_offsets (t8_forest_t forest)
{

  // Initialization.
  const t8_forest_t forest_old = forest->set_from; /* committed */
  const t8_shmem_array_t partition_new = forest->element_offsets;
  std::vector<t8_gloidx_t> corrected_local_offsets;

  // Nothing to be done for empty processes.
  if (t8_forest_get_local_num_leaf_elements (forest_old) != 0) {

    // Send requests to other processes.
    std::vector<sc_MPI_Request> requests;
    t8_forest_pfc_send_loop_range<t8_forest_pfc_message> (forest_old, requests);

    // Receive messages from other processes.
    std::vector<t8_forest_pfc_message> messages;
    t8_forest_pfc_recv_loop_range<t8_forest_pfc_message> (forest_old, messages);

    // Wait for Isend requests.
    const int mpiret = sc_MPI_Waitall (requests.size (), requests.data (), sc_MPI_STATUSES_IGNORE);
    SC_CHECK_MPI (mpiret);

    // Compute process-local corrections to partitioning
    t8_forest_pfc_correct_local_offsets (forest_old, partition_new, messages, corrected_local_offsets);
  }

  // Allgatherv the corrected local offsets, resulting in the new partition.
  t8_shmem_array_allgatherv (corrected_local_offsets.data (), corrected_local_offsets.size (), T8_MPI_GLOIDX,
                             partition_new, T8_MPI_GLOIDX, forest->mpicomm);
}
T8_EXTERN_C_END ();
