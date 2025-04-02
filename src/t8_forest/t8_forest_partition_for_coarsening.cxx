#include <t8_forest/t8_forest.h>
#include <t8_forest/t8_forest_types.h>
#include <t8_forest/t8_forest_private.h>
#include <t8_forest/t8_forest_pfc_message.hxx>
#include <t8_forest/t8_forest_pfc_helper.hxx>
#include <t8_data/t8_shmem.h>
#include <t8_element_cxx.hxx>
#include <vector>
#include <algorithm>

t8_procidx_t
proc_owner (const t8_gloidx_t *partition, t8_procidx_t mpisize, t8_gloidx_t element_id)
{
  /*std::upper_bound gives the first boundary that is bigger than element_id, so subtracting one gives the last boundary, that is less than or equal to it , so its process owns element */
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
t8_forest_max_num_children (t8_forest_t forest)
{
  return 10;
}

template <typename MessageType>
static void
t8_forest_pfc_send_loop_range (t8_forest_t forest, std::vector<sc_MPI_Request> &requests)
{
  T8_ASSERT (t8_forest_is_committed (forest));
  T8_ASSERT (requests.size () == 0);
  t8_procidx_t rank = forest->mpirank;
  t8_procidx_t mpisize = forest->mpisize;

  const t8_gloidx_t *partition = t8_shmem_array_get_gloidx_array (forest->element_offsets);
  const int max_num_siblings = t8_forest_max_num_children (forest);
  /* why the difference between directions? */
  t8_gloidx_t relevant_begin = SC_MAX (0, partition[rank] - (max_num_siblings - 1));
  t8_gloidx_t relevant_end = SC_MIN (partition[mpisize], partition[rank + 1] + max_num_siblings);

  t8_procidx_t begin = proc_owner (partition, mpisize, relevant_begin);
  t8_procidx_t end = proc_owner_end (partition, mpisize, relevant_end);
  T8_ASSERT (0 <= begin);
  T8_ASSERT (begin <= end);
  T8_ASSERT (end <= forest->mpisize);

  for (t8_procidx_t iproc = begin; iproc < end; iproc++) {
    if (partition[iproc] >= partition[iproc + 1] || iproc == forest->mpirank)
      continue;

    /* Fill message */
    MessageType message (forest->scheme_cxx, iproc, forest->mpicomm);
    message.fill (forest);

    sc_MPI_Request request;
    message.mpi_Isend (forest, request);

    requests.push_back (std::move (request));
  }
}

template <typename MessageType>
static void
t8_forest_pfc_recv_loop_range (t8_forest_t forest, std::vector<MessageType> &messages)
{
  T8_ASSERT (t8_forest_is_committed (forest));
  T8_ASSERT (messages.size () == 0);
  t8_procidx_t rank = forest->mpirank;
  t8_procidx_t mpisize = forest->mpisize;

  const t8_gloidx_t *partition = t8_shmem_array_get_gloidx_array (forest->element_offsets);
  const int max_num_siblings = t8_forest_max_num_children (forest);
  /* why the difference between directions? */
  t8_gloidx_t relevant_begin = SC_MAX (0, partition[rank] - max_num_siblings);
  t8_gloidx_t relevant_end = SC_MIN (partition[mpisize], partition[rank + 1] + (max_num_siblings - 1));

  t8_procidx_t begin = proc_owner (partition, mpisize, relevant_begin);
  t8_procidx_t end = proc_owner_end (partition, mpisize, relevant_end);
  T8_ASSERT (0 <= begin);
  T8_ASSERT (begin <= end);
  T8_ASSERT (end <= forest->mpisize);

  for (t8_procidx_t iproc = begin; iproc < end; iproc++) {
    if (partition[iproc] >= partition[iproc + 1] || iproc == forest->mpirank)
      continue;

    MessageType message (forest->scheme_cxx, iproc, forest->mpicomm);
    t8_debugf ("receive message from %i\n", message.iproc);
    int buf_size;
    char *recv_buf;
    message.mpi_Recv (recv_buf, buf_size);
    int position = 0;
    message.unpack (recv_buf, buf_size, &position);
    T8_ASSERT (position == buf_size);
    T8_FREE (recv_buf);
    messages.push_back (std::move (message));
  }
}

static int
t8_forest_pfc_family_range_around_border (t8_forest_t forest, t8_gloidx_t border_element_id,
                                          std::vector<t8_forest_pfc_message_c> &messages, t8_gloidx_t &family_begin,
                                          t8_gloidx_t &family_end)
{
  t8_gloidx_t gtree_id;
  t8_tree_t tree;
  t8_eclass_scheme_c *scheme;
  t8_element_t *parent, *element;
  t8_locidx_t index_in_tree;
  t8_forest_pfc_helper_index_in_tree_from_globalid (forest, border_element_id, gtree_id, scheme, tree, index_in_tree,
                                                    element);

  if (scheme->t8_element_level (element) == 0) {
    family_begin = border_element_id;
    family_end = border_element_id;
    return false;
  }

  scheme->t8_element_new (1, &parent);
  scheme->t8_element_parent (element, parent);

  t8_gloidx_t first_tree_element = t8_forest_get_first_local_element_id (forest) + tree->elements_offset;
  family_begin = first_tree_element + t8_forest_pfc_extreme_local_sibling (scheme, tree, index_in_tree, true);
  /* end iterator is one behind last element */
  family_end = first_tree_element + t8_forest_pfc_extreme_local_sibling (scheme, tree, index_in_tree, false) + 1;

  /* check if other processes have the same parent as we do, so we need to adjust our range */
  for (t8_procidx_t imessage = 0; imessage < (t8_procidx_t) messages.size (); imessage++) {
    t8_debugf ("process message from %i\n", messages[imessage].iproc);
    /* on the same tree we can use our scheme to compare, because we know that the eclasses are equal */
    if (messages[imessage].itree == gtree_id && scheme->t8_element_equal (parent, messages[imessage].get_parent ())) {
      /* if our parents are equal, adjust left or right range border by num_siblings*/
      if (messages[imessage].iproc < forest->mpirank) {
        family_begin -= messages[imessage].num_siblings;
      }
      else {
        family_end += messages[imessage].num_siblings;
      }
    }
  }

  int num_children = scheme->t8_element_num_children (parent);
  scheme->t8_element_destroy (1, &parent);
  return (family_end - family_begin == num_children);
}

/* Maybe replace by all to rank with most elements */
static int
t8_forest_pfc_family_split_rank_all_to_first (t8_shmem_array_t partition_new_shmem, int rank, t8_gloidx_t family_begin,
                                              t8_gloidx_t family_end)
{
  int num_ranks = t8_shmem_array_get_elem_count (partition_new_shmem);
  const t8_gloidx_t *partition_new = t8_shmem_array_get_gloidx_array (partition_new_shmem);
  const t8_gloidx_t *it = std::lower_bound (partition_new, partition_new + num_ranks, family_begin);
  return (it - partition_new);
}

static void
t8_forest_pfc_correct_local_offsets (t8_forest_t forest, t8_shmem_array_t partition_new_shmem,
                                     std::vector<t8_forest_pfc_message_c> &messages,
                                     std::vector<t8_gloidx_t> &corrected_local_offsets)
{
  T8_ASSERT (t8_forest_is_committed (forest));
  const t8_gloidx_t *partition_new = t8_shmem_array_get_gloidx_array (partition_new_shmem);
  const t8_gloidx_t *min_local_element_pointer
    = std::lower_bound (partition_new, partition_new + forest->mpisize,
                        t8_shmem_array_get_gloidx (forest->element_offsets, forest->mpirank));
  const t8_gloidx_t *next_min_local_element_pointer
    = std::lower_bound (partition_new, partition_new + forest->mpisize,
                        t8_shmem_array_get_gloidx (forest->element_offsets, forest->mpirank + 1));

  const t8_gloidx_t min_local_proc = min_local_element_pointer - partition_new;
  const t8_gloidx_t next_min_local_proc = next_min_local_element_pointer - partition_new;

  /* adjust all local borders */
  for (t8_procidx_t border_rank = min_local_proc; border_rank < next_min_local_proc; border_rank++) {
    t8_gloidx_t family_begin, family_end;
    if (t8_forest_pfc_family_range_around_border (forest, partition_new[border_rank], messages, family_begin,
                                                  family_end)) {
      /* border needs to be adjusted */
      t8_procidx_t rank = forest->mpirank;
      /* determine rank that gets all elements */
      t8_procidx_t split_rank
        = t8_forest_pfc_family_split_rank_all_to_first (forest->element_offsets, rank, family_begin, family_end);
      /* correct local offset, possible TODO: update all local_offsets affected by this family */
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
  const t8_forest_t forest_old = forest->set_from; /* committed */
  const t8_shmem_array_t partition_new = forest->element_offsets;

  /* Send */
  std::vector<sc_MPI_Request> requests;
  t8_forest_pfc_send_loop_range<t8_forest_pfc_message_c> (forest_old, requests);

  std::vector<t8_forest_pfc_message_c> messages;
  t8_forest_pfc_recv_loop_range<t8_forest_pfc_message_c> (forest_old, messages);

  /* Wait for Isend requests */
  int mpiret = sc_MPI_Waitall (requests.size (), requests.data (), sc_MPI_STATUSES_IGNORE);
  SC_CHECK_MPI (mpiret);

  /* Compute local corrections to partition */
  std::vector<t8_gloidx_t> corrected_local_offsets;
  t8_forest_pfc_correct_local_offsets (forest_old, partition_new, messages, corrected_local_offsets);

  /** Allgather nums of local correnction */
  t8_shmem_array_t global_correction_counts;
  t8_shmem_array_init (&global_correction_counts, sizeof (int), forest->mpisize, forest->mpicomm);
  int num_corrections = corrected_local_offsets.size ();
  t8_shmem_array_allgather (&num_corrections, 1, sc_MPI_INT, global_correction_counts, forest->mpisize, sc_MPI_INT);

  /** Allgatherv the actual corrected offsets */
  t8_shmem_array_allgatherv (corrected_local_offsets.data (), corrected_local_offsets.size (), T8_MPI_GLOIDX,
                             partition_new, T8_MPI_GLOIDX, forest->mpicomm);

  /* Memory cleanup */
  t8_shmem_array_destroy (&global_correction_counts);
}
T8_EXTERN_C_END ();
