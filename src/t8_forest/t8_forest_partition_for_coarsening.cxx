#include <t8_forest/t8_forest.h>
#include <t8_forest/t8_forest_types.h>
#include <t8_forest/t8_forest_private.h>
#include <t8_forest/t8_forest_pfc_message.hxx>
#include <t8_forest/t8_forest_pfc_helper.hxx>
#include <t8_data/t8_shmem.h>
#include <t8_element.h>
#include <t8_forest/t8_forest_partition_for_coarsening.h>
#include <vector>
#include <algorithm>

/** Return the process owner of the given (global) element_id
 * 
 * \param[in]   partition   the current partitioning, given as array of element offsets.
 * \param[in]   mpisize     the number of MPI ranks
 * \param[in]   element_id  the global index of the element considered
 * 
 * \return The ID of the process owning the given element.
 **/
t8_procidx_t
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

t8_procidx_t
proc_owner_end (const t8_gloidx_t *partition, t8_procidx_t mpisize, t8_gloidx_t element_end)
{
  /*std::upper_bound gives the first boundary that is bigger than element_id, so it is the process after the one owning element_id */
  t8_gloidx_t element_id = element_end - 1;
  return std::upper_bound (partition, partition + mpisize, element_id) - partition;
}

int
t8_forest_max_num_children ([[maybe_unused]] t8_forest_t forest)
{
  return 10;
}

/** 
 *  Send PFC messages to all relevant processes and obtain the associated requests.
 * 
 * \param[in]   forest    the current forest
 * \param[out]  requests  the MPI requests as std::vector
 *                          on input:  empty
 *                          on output: contains the send requests
*/
template <typename MessageType>
static void
t8_forest_pfc_send_loop_range (const t8_forest_t forest, std::vector<sc_MPI_Request> &requests)
{
  // Assertions: The forest must be committed and the request vector empty
  T8_ASSERT (t8_forest_is_committed (forest));
  T8_ASSERT (requests.size () == 0);

  // Initializations
  t8_procidx_t rank = forest->mpirank;
  t8_procidx_t mpisize = forest->mpisize;
  const int max_num_siblings = t8_forest_max_num_children (forest);

  // Get current offset vector of forest.
  const t8_gloidx_t *partition = t8_shmem_array_get_gloidx_array (forest->element_offsets);

  // Get range of elements that may end up on the current process due to the correction:
  //    For that, bounds of the current (equally-sized) partitioning have to be extended by
  //    the maximum number of siblings that may form a family plus one (because MYTODO)
  //    (Note: The SC_MAX and SC_MIN commands are only relevant for the first and last process.)
  t8_gloidx_t relevant_begin = SC_MAX (0, partition[rank] - (max_num_siblings - 1));
  t8_gloidx_t relevant_end = SC_MIN (partition[mpisize], partition[rank + 1] + max_num_siblings);

  // Determine range of processors holding the range of relevant elements.
  t8_procidx_t begin = proc_owner (partition, mpisize, relevant_begin);
  t8_procidx_t end = proc_owner_end (partition, mpisize, relevant_end);
  T8_ASSERT (0 <= begin);
  T8_ASSERT (begin <= end);
  T8_ASSERT (end <= forest->mpisize);

  // Loop over processes of relevant range.
  for (t8_procidx_t iproc = begin; iproc < end; iproc++) {

    // Skip empty processes and the own rank.
    if (partition[iproc] >= partition[iproc + 1] || iproc == forest->mpirank)
      continue;

    // Construct and fill message of type MessageType (see t8_forest_pfc_message_c).
    MessageType message (forest->scheme, iproc, forest->mpicomm);
    message.fill (forest);

    // Send the message (and obtain the associated requests).
    sc_MPI_Request request;
    message.mpi_Isend (forest, request);

    // Add to requests array
    requests.push_back (std::move (request));
  }
}

/** Determine the messages to be received.
 * 
 * \param[in]   forest    the current forest
 * \param[out]  requests  the MPI messages
*/
template <typename MessageType>
static void
t8_forest_pfc_recv_loop_range (const t8_forest_t forest, std::vector<MessageType> &messages)
{
  // Assertions: forest must be committed and messages empty
  T8_ASSERT (t8_forest_is_committed (forest));
  T8_ASSERT (messages.size () == 0);

  // Initialization
  t8_procidx_t rank = forest->mpirank;
  t8_procidx_t mpisize = forest->mpisize;
  const int max_num_siblings = t8_forest_max_num_children (forest);

  // Get current offset vector of forest.
  const t8_gloidx_t *partition = t8_shmem_array_get_gloidx_array (forest->element_offsets);

  // Get range of elements that may end up on the current process due to the correction:
  //    For that, bounds of the current (equally-sized) partitioning have to be extended by
  //    the maximum number of siblings that may form a family plus one (because MYTODO)
  //    (Note: The SC_MAX and SC_MIN commands are only relevant for the first and last process.)
  t8_gloidx_t relevant_begin = SC_MAX (0, partition[rank] - max_num_siblings);
  t8_gloidx_t relevant_end = SC_MIN (partition[mpisize], partition[rank + 1] + (max_num_siblings - 1));
  t8_procidx_t begin = proc_owner (partition, mpisize, relevant_begin);
  t8_procidx_t end = proc_owner_end (partition, mpisize, relevant_end);
  T8_ASSERT (0 <= begin);
  T8_ASSERT (begin <= end);
  T8_ASSERT (end <= forest->mpisize);

  // Loop over process range
  for (t8_procidx_t iproc = begin; iproc < end; iproc++) {
    if (partition[iproc] >= partition[iproc + 1] || iproc == forest->mpirank)
      continue;

    // Receive message.
    MessageType message (forest->scheme, iproc, forest->mpicomm);
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

/** Determine whether a full family is split by a process boundary.
 * 
 * \param[in]   forest              the forest
 * \param[in]   border_element_id   the global ID of the border element
 * \param[in]   messages            the PFC messages received from other processes
 * \param[out]  family_begin        the global element ID of the family's first member
 * \param[out]  family_end          the global element ID of the family's last member
 * 
 * \return True (i.e., nonzero) if and only if a full family is found across the process borders. 
*/
static int
t8_forest_pfc_family_range_around_border (const t8_forest_t forest, const t8_gloidx_t border_element_id,
                                          const std::vector<t8_forest_pfc_message_c> &messages,
                                          t8_gloidx_t &family_begin, t8_gloidx_t &family_end)
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
  t8_forest_pfc_helper_index_in_tree_from_globalid (forest, border_element_id, gtree_id, tree, index_in_tree, element);

  // Get scheme and eclass from forest and tree
  const t8_scheme_c *newscheme = t8_forest_get_scheme (forest);
  t8_eclass_t eclass = tree->eclass;

  // If the element is the root, return false because the root does not have any parent or siblings.
  if (newscheme->element_get_level (eclass, element) == 0) {
    family_begin = border_element_id;
    family_end = border_element_id;
    return false;
  }

  // Allocate and determine parent element.
  t8_element_t *parent;
  t8_element_new (newscheme, eclass, 1, &parent);
  newscheme->element_get_parent (eclass, element, parent);

  // Get global ID of first (process-)local element
  t8_gloidx_t first_tree_element = t8_forest_get_first_local_leaf_element_id (forest) + tree->elements_offset;

  // TODO: Isn't the tree offset always zero because it is the first tree?
  T8_ASSERT(tree->elements_offset==0);

  // Determine range of global IDs forming the family of first_tree_element, by calling the helper function
  // t8_forest_pfc_extreme_local_sibling twice, i.e., searching in the direction of in- and decreasing indices.
  // Note: The end iterator is one behind the last family member.
  family_begin = first_tree_element + t8_forest_pfc_extreme_local_sibling (newscheme, tree, index_in_tree, true);
  family_end = first_tree_element + t8_forest_pfc_extreme_local_sibling (newscheme, tree, index_in_tree, false) + 1;

  // Check if other processes have the same parent as the current family, so we need to adjust our range
  for (t8_procidx_t imessage = 0; imessage < (t8_procidx_t) messages.size (); imessage++) {
    t8_debugf ("process message from %i\n", messages[imessage].iproc);

    // On the same tree we can use our scheme to compare, because we know that the eclasses are equal.
    if (messages[imessage].itree == gtree_id
        && newscheme->element_is_equal (eclass, parent, messages[imessage].get_parent ())) {

      // If parents are equal, extend lower or upper range border by num_siblings, depending on the send "direction",
      // i.e., towards lower- or higher-rank processes.
      if (messages[imessage].iproc < forest->mpirank) {
        family_begin -= messages[imessage].num_siblings;
      }
      else {
        family_end += messages[imessage].num_siblings;
      }
    }
  }

  // Determine the parent's number of children.
  int num_children = newscheme->element_get_num_children (eclass, parent);

  // Deallocate parent element
  t8_element_destroy (newscheme, eclass, 1, &parent);

  // Return true if the considered family contains all children of the parent.
  return (family_end - family_begin == num_children);
}

/*  */

/** 
 * 
 * 
 * Possible Todo: Replace by all to rank with most elements
*/
static int
t8_forest_pfc_family_split_rank_all_to_first (const t8_shmem_array_t partition_new_shmem,
                                              [[maybe_unused]] const int rank, const t8_gloidx_t family_begin,
                                              [[maybe_unused]] const t8_gloidx_t family_end)
{
  // Get number of MPI ranks and the partition array
  int num_ranks = t8_shmem_array_get_elem_count (partition_new_shmem);
  const t8_gloidx_t *partition_new = t8_shmem_array_get_gloidx_array (partition_new_shmem);

  // Determine which processor holds the element with global ID family_begin according to new partition.
  const t8_gloidx_t *it = std::lower_bound (partition_new, partition_new + num_ranks, family_begin);

  // Return the processor index.
  return (it - partition_new);
}

/** Compute the process-local corrections of the given partition.
 * 
 * \param[in]   forest                  the forest
 * \param[in]   partition_new_shmem     the current partitioning (without PFC correcton) as shared-memory array 
 * \param[in]   messages                the PFC messages received from other processes
 * \param[out]  corrected_local_offsets a std::vector of t8_gloidx_t>
 *                                      on input:  empty
 *                                      on output: containing the corrections to be applied to the local offsets to obtain the PFC partitioning
*/
static void
t8_forest_pfc_correct_local_offsets (const t8_forest_t forest, const t8_shmem_array_t partition_new_shmem,
                                     const std::vector<t8_forest_pfc_message_c> &messages,
                                     std::vector<t8_gloidx_t> &corrected_local_offsets)
{
  T8_ASSERT (t8_forest_is_committed (forest));

  // Get current partitioning as array of t8_gloidx_t.
  const t8_gloidx_t *partition_new = t8_shmem_array_get_gloidx_array (partition_new_shmem);

  // Determine on which process the first IDs of the old partitioning would be according to the new one.

  const t8_gloidx_t *min_local_element_pointer
    = std::lower_bound (partition_new, partition_new + forest->mpisize,
                        t8_shmem_array_get_gloidx (forest->element_offsets, forest->mpirank));
  const t8_gloidx_t *next_min_local_element_pointer
    = std::lower_bound (partition_new, partition_new + forest->mpisize,
                        t8_shmem_array_get_gloidx (forest->element_offsets, forest->mpirank + 1));
  const t8_gloidx_t min_local_proc = min_local_element_pointer - partition_new;
  const t8_gloidx_t next_min_local_proc = next_min_local_element_pointer - partition_new;

  // TODO: Is the above maybe only needed for empty processes?
  // Otherwise simplify?
  T8_ASSERT(t8_shmem_array_get_gloidx (forest->element_offsets, forest->mpirank) != t8_shmem_array_get_gloidx (forest->element_offsets, forest->mpirank + 1));
  // T8_ASSERT(min_local_proc == forest->mpirank);
  T8_ASSERT(next_min_local_proc == forest->mpirank+1);

  /* adjust all local borders */
  for (t8_procidx_t border_rank = min_local_proc; border_rank < next_min_local_proc; border_rank++) {
    t8_gloidx_t family_begin, family_end;

    // Check if there is a full family split by the current border.
    if (t8_forest_pfc_family_range_around_border (forest, partition_new[border_rank], messages, family_begin,
                                                  family_end)) {
      // border needs to be adjusted
      t8_procidx_t rank = forest->mpirank;
      // determine rank that gets all elements
      t8_procidx_t split_rank
        = t8_forest_pfc_family_split_rank_all_to_first (forest->element_offsets, rank, family_begin, family_end);
      // correct local offset, possible TODO: update all local_offsets affected by this family
      t8_gloidx_t new_offset = (border_rank <= split_rank) ? family_begin : family_end;
      corrected_local_offsets.push_back (new_offset);
    }
    else {
      /* no correction needed*/
      t8_gloidx_t new_offset = partition_new[border_rank];
      corrected_local_offsets.push_back (new_offset);
    }
  }
}

T8_EXTERN_C_BEGIN ();
void
t8_forest_pfc_correction_offsets (t8_forest_t forest)
{

  // Initialization
  const t8_forest_t forest_old = forest->set_from; /* committed */
  const t8_shmem_array_t partition_new = forest->element_offsets;
  std::vector<t8_gloidx_t> corrected_local_offsets;

  // Nothing to be done for empty processes.
  if (t8_forest_get_local_num_leaf_elements (forest_old) != 0) {
    
    // Send requests to other processes.
    std::vector<sc_MPI_Request> requests;
    t8_forest_pfc_send_loop_range<t8_forest_pfc_message_c> (forest_old, requests);

    // Receive messages from other processes.
    std::vector<t8_forest_pfc_message_c> messages;
    t8_forest_pfc_recv_loop_range<t8_forest_pfc_message_c> (forest_old, messages);

    // Wait for Isend requests.
    int mpiret = sc_MPI_Waitall (requests.size (), requests.data (), sc_MPI_STATUSES_IGNORE);
    SC_CHECK_MPI (mpiret);

    // Compute process-local corrections to partitioning
    t8_forest_pfc_correct_local_offsets (forest_old, partition_new, messages, corrected_local_offsets);
  }

  // Allgatherv the corrected local offsets, resulting in the new partition.
  t8_shmem_array_allgatherv (corrected_local_offsets.data (), corrected_local_offsets.size (), T8_MPI_GLOIDX,
                             partition_new, T8_MPI_GLOIDX, forest->mpicomm);

}
T8_EXTERN_C_END ();
