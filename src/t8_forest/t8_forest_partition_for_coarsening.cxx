#include <t8_forest/t8_forest.h>
#include <t8_forest/t8_forest_types.h>
#include <t8_forest/t8_forest_private.h>
#include <t8_data/t8_shmem.h>
#include <t8_element_cxx.hxx>
#include <vector>
#include <algorithm>

typedef int t8_procidx_t;
#define T8_MAX_CHILDREN 10

#define T8_PFC_MESSAGE 54323

/** Determine if we need to send information to a proc, so that they can determine if they need to adjust their lower bound in the elements_offset (partition)
 * Skip empty procs 
 * TODO: replace by binary search */
static void
t8_forest_pfc_determine_send_recv_range (t8_shmem_array_t partition_shmem, int rank, int *begin, int *end,
                                         int *num_interactions)
{

  int num_procs = t8_shmem_array_get_elem_count (partition_shmem) - 1;
  const t8_gloidx_t *partition = t8_shmem_array_get_gloidx_array (partition_shmem);

  *num_interactions = 0;
  if (partition[rank] >= partition[rank + 1])
    return; /* no elements, so we cannot send any */

  t8_gloidx_t min_relevant_element_id = SC_MAX (0, partition[rank] - T8_MAX_CHILDREN + 1);
  t8_debugf ("min_relevant_element_id %i\n", min_relevant_element_id);
  const t8_gloidx_t *min_relevant_element_id_pointer
    = std::upper_bound (partition, partition + num_procs, min_relevant_element_id) - 1;
  *begin = min_relevant_element_id_pointer - partition;

  t8_gloidx_t max_relevant_element_id = SC_MIN (partition[num_procs], partition[rank + 1] + T8_MAX_CHILDREN - 1);
  const t8_gloidx_t *max_relevant_element_id_pointer
    = std::lower_bound (partition, partition + num_procs, max_relevant_element_id);
  *end = max_relevant_element_id_pointer - partition;

  int num_nonzero = 0;
  for (const t8_gloidx_t *it = partition + *begin; it < partition + *end; it++) {
    if (*it < *(it + 1))
      num_nonzero++;
  }
  *num_interactions = num_nonzero - 1;
  t8_debugf ("begin: %i, end: %i, num_interactions: %i\n", *begin, *end, *num_interactions);

#if 0
  /* Unreachable proc_ids as initialisation for min and max search*/
  *send_lowest = num_procs;
  *send_highest = -1;
  for (t8_procidx_t iproc = 0; iproc < num_procs; iproc++) {
    /* Check if this process has relevant elements for process iproc */
    if (partition[iproc] < partition[iproc + 1] && /* This process has elements, so cannot decide about any border */
        /* The start element of iproc in the new partition is not too far to the left of the old range of our proc*/
        partition[rank] + 2 - T8_MAX_CHILDREN <= partition[iproc] && 
        /* The start element of iproc in the new partition is not too far to the right of the old range of our proc*/
        partition[iproc] < partition[rank + 1] + T8_MAX_CHILDREN - 1) {

      /* We need to send elements to process iproc, because their new start element is not too far away from our old elements */
      (*num_sends)++;
      *send_lowest = SC_MIN(*send_lowest, iproc);
      /* Do not need a max here because we traverse loop in order */
      *send_highest = iproc;
    }
  }
#endif
}

/** Determine the range of elements that we send/receive information from, i.e. they have elements close enough to our range, and are nonempty */

static void
t8_forest_pfc_determine_send_recv_range_old (t8_shmem_array_t partition_shmem, int rank, t8_procidx_t *receive_lowest,
                                             t8_procidx_t *receive_highest, t8_procidx_t *num_recvs)
{
  int num_procs = t8_shmem_array_get_elem_count (partition_shmem) - 1;
  const t8_gloidx_t *partition = t8_shmem_array_get_gloidx_array (partition_shmem);

  /* We do not receive any elements if we do not have any elements in the new partition*/
  *num_recvs = -1;
  if (partition[rank] >= partition[rank + 1])
    return;

  t8_gloidx_t min_relevant_element_id = SC_MAX (0, partition[rank] - T8_MAX_CHILDREN + 1);
  t8_gloidx_t max_relevant_element_id = SC_MIN (partition[num_procs], partition[rank + 1] + T8_MAX_CHILDREN - 1);
  t8_debugf ("min: %i, max:%i\n", min_relevant_element_id, max_relevant_element_id);

  /* Unreachable proc_ids as initialisation for min and max search*/
  *receive_lowest = num_procs;
  *receive_highest = -1;
  for (t8_procidx_t iproc = 0; iproc < num_procs; iproc++) {

    if (partition[iproc] < partition[iproc + 1] && /*This process has elements*/
        /* The start element of our proc in the partition is not too far to the left of the range of iproc*/
        min_relevant_element_id < partition[iproc + 1] &&
        /* The start element of our proc in the partition is not too far to the right of the range of iproc*/
        partition[iproc] < max_relevant_element_id) {
      t8_debugf ("partition[%i] = %i\n", rank, partition[rank]);
      /* process iproc has elements, that need to be send to our proc, because their old elements are not too far away from our new start element*/
      (*num_recvs)++;
      *receive_lowest = SC_MIN (*receive_lowest, iproc);
      /* Do not need a max here because we traverse loop in order */
      *receive_highest = iproc;
    }
  }
}

static t8_locidx_t
t8_forest_pfc_extreme_sibling (t8_eclass_scheme_c *scheme, t8_tree_t tree, t8_locidx_t start_element_id_in_tree,
                               int reverse)
{
  t8_debugf ("tree: %p, elid_intree: %i, reverse: %i\n", tree, start_element_id_in_tree, reverse);
  t8_element_t *parent_possible_sibling, *parent_start, *start_element;
  scheme->t8_element_new (1, &parent_possible_sibling);
  scheme->t8_element_new (1, &parent_start);

  start_element = t8_forest_get_tree_element (tree, start_element_id_in_tree);
  scheme->t8_element_debug_print (start_element);

  scheme->t8_element_parent (start_element, parent_start);
  int num_children = scheme->t8_element_num_children (parent_start);

  t8_locidx_t extreme_sibling_id_in_tree = start_element_id_in_tree;
  t8_locidx_t extreme_check_id_in_tree
    = reverse ? SC_MAX (0, start_element_id_in_tree - num_children)
              : SC_MIN (start_element_id_in_tree + num_children, t8_forest_get_tree_element_count (tree)) - 1;
  int increment = reverse ? -1 : 1;

  for (t8_locidx_t ielem = start_element_id_in_tree; (ielem - extreme_check_id_in_tree) * increment <= 0;
       ielem += increment) {
    const t8_element_t *possible_sibling = t8_forest_get_tree_element (tree, ielem);
    if (scheme->t8_element_level (possible_sibling)) { /*We are not allowed to compute parent of root*/
      scheme->t8_element_parent (possible_sibling, parent_possible_sibling);
      if (scheme->t8_element_equal (parent_start, parent_possible_sibling)) {
        extreme_sibling_id_in_tree = ielem;
      }
      else {
        break;
      }
    }
  }
  scheme->t8_element_destroy (1, &parent_possible_sibling);
  scheme->t8_element_destroy (1, &parent_start);
  return extreme_sibling_id_in_tree;
}

static void
t8_forest_pfc_helper_global_elementid_get_stuff (t8_forest_t forest, t8_gloidx_t gelement_id, t8_gloidx_t *gtree_id,
                                                 t8_eclass_scheme_c **scheme, t8_tree_t *tree,
                                                 t8_locidx_t *index_in_tree, t8_element_t **element)
{
  t8_gloidx_t global_id_of_first_local_element = t8_forest_get_first_local_element_id (forest);

  t8_locidx_t lelement_id = (t8_locidx_t) (gelement_id - global_id_of_first_local_element);

  t8_locidx_t ltree_id;
  *element = t8_forest_get_element (forest, lelement_id, &ltree_id);

  *tree = t8_forest_get_tree (forest, ltree_id);
  *index_in_tree = lelement_id - (*tree)->elements_offset;
  T8_ASSERT (*element == t8_forest_get_tree_element (*tree, *index_in_tree));

  t8_locidx_t first_local_tree_id = t8_forest_get_first_local_tree_id (forest);
  *gtree_id = first_local_tree_id + ltree_id;
  t8_eclass_t eclass = (*tree)->eclass;
  *scheme = t8_forest_get_eclass_scheme (forest, eclass);
}

/* get min element id of the  with same parent as start*/
static t8_gloidx_t
t8_forest_pfc_helper_min_local_sibling_closest_to_start (t8_eclass_scheme_c *scheme, t8_tree_t tree,
                                                         t8_locidx_t index_in_tree)
{
  const int reverse = 1;
  return t8_forest_pfc_extreme_sibling (scheme, tree, index_in_tree, reverse);
}

/* get min element id of the  with same parent as start*/
static t8_gloidx_t
t8_forest_pfc_helper_max_local_sibling_closest_to_start (t8_eclass_scheme_c *scheme, t8_tree_t tree,
                                                         t8_locidx_t index_in_tree)
{
  const int reverse = 0;
  return t8_forest_pfc_extreme_sibling (scheme, tree, index_in_tree, reverse);
}

class t8_forest_pfc_message_c {
 public:
  t8_gloidx_t itree;  /*actually sent*/
  t8_eclass_t eclass; /*actually sent*/
  int num_siblings;   /*actually sent*/
  t8_eclass_scheme_c *scheme;
  void
  pack (void *buf, int buf_size, int *position, sc_MPI_Comm comm);
  void
  unpack (void *buf, int buf_size, int *position, sc_MPI_Comm comm, t8_scheme_cxx_t *scheme_cxx);
  int
  pack_size (sc_MPI_Comm comm);
  void
  pack_and_send (t8_procidx_t iproc, sc_MPI_Comm comm, sc_MPI_Request *request);
  void
  receive_and_unpack (t8_procidx_t iproc, sc_MPI_Comm comm);
  void
  fill (t8_forest_t forest, t8_procidx_t iproc);
  const t8_element_t *
  get_parent ()
  {
    T8_ASSERT (parent);
    return parent;
  }
  t8_forest_pfc_message_c ()
    : itree (0), eclass (T8_ECLASS_ZERO), num_siblings (0), scheme (NULL), parent (NULL), allocated_parent (0)
  {
  }
  t8_forest_pfc_message_c (const t8_forest_pfc_message_c &other) = delete;
  t8_forest_pfc_message_c (t8_forest_pfc_message_c &&other)
    : itree { other.itree }, eclass (other.eclass), num_siblings (other.num_siblings), scheme (other.scheme),
      parent (other.parent), allocated_parent (other.allocated_parent)
  {
    if (allocated_parent) {
      other.parent = NULL;
      other.allocated_parent = false;
    }
  }
  ~t8_forest_pfc_message_c ()
  {
    if (allocated_parent) {
      scheme->t8_element_destroy (1, &parent);
      parent = NULL; /*necessary?*/
    }
  }

 private:
  t8_element_t *parent; /*is unpacked and then sent*/
  bool allocated_parent;
};

void
t8_forest_pfc_message_c::pack (void *buf, int buf_size, int *position, sc_MPI_Comm comm)
{
  sc_MPI_Pack (&itree, 1, T8_MPI_GLOIDX, buf, buf_size, position, comm);
  t8_debugf ("packed itree %i\n", itree);

  int eclass_int = (int) eclass;
  sc_MPI_Pack (&eclass_int, 1, sc_MPI_INT, buf, buf_size, position, comm);
  t8_debugf ("packed eclass %i\n", eclass);

  scheme->t8_element_pack (parent, 1, buf, buf_size, position, comm);
  t8_debugf ("packed parent\n");
  scheme->t8_element_debug_print (parent);

  sc_MPI_Pack (&num_siblings, 1, sc_MPI_INT, buf, buf_size, position, comm);
  t8_debugf ("packed num_siblings %i\n", num_siblings);
}

/**Allocates memory for the element, if it is from the same eclass if our scheme */
void
t8_forest_pfc_message_c::unpack (void *buf, int buf_size, int *position, sc_MPI_Comm comm, t8_scheme_cxx_t *scheme_cxx)
{
  sc_MPI_Unpack (buf, buf_size, position, &itree, 1, T8_MPI_GLOIDX, comm);
  t8_debugf ("unpacked itree %i\n", itree);
  int eclass_int;
  sc_MPI_Unpack (buf, buf_size, position, &eclass_int, 1, sc_MPI_INT, comm);
  eclass = (t8_eclass_t) eclass_int;
  t8_debugf ("unpacked eclass %i\n", eclass);
  scheme = scheme_cxx->eclass_schemes[eclass_int];
  scheme->t8_element_new (1, &parent);
  allocated_parent = true;
  scheme->t8_element_unpack (buf, buf_size, position, parent, 1, comm);
  t8_debugf ("unpacked parent\n");
  scheme->t8_element_debug_print (parent);
  sc_MPI_Unpack (buf, buf_size, position, &num_siblings, 1, sc_MPI_INT, comm);
  t8_debugf ("unpacked num_siblings %i\n", num_siblings);
}

/* gives the pack size for this element, i.e. only if there actually is an element, its space needs are added. */
int
t8_forest_pfc_message_c::pack_size (sc_MPI_Comm comm)
{
  int message_size = 0;
  int datasize;
  sc_MPI_Pack_size (1, T8_MPI_GLOIDX, comm, &datasize);
  message_size += datasize;
  sc_MPI_Pack_size (1, sc_MPI_INT, comm, &datasize);
  message_size += datasize;
  if (parent != NULL) {
    scheme->t8_element_pack_size (1, comm, &datasize);
    message_size += datasize;
  }
  sc_MPI_Pack_size (1, sc_MPI_INT, comm, &datasize);
  message_size += datasize;
  return message_size;
}

void
t8_forest_pfc_message_c::pack_and_send (t8_procidx_t iproc, sc_MPI_Comm comm, sc_MPI_Request *request)
{
  int position = 0;
  int buf_size = pack_size (comm);
  char *buf = T8_ALLOC (char, buf_size);

  pack (buf, buf_size, &position, comm);

  int mpiret = sc_MPI_Isend (buf, position, sc_MPI_PACKED, iproc, T8_PFC_MESSAGE, comm, request);
  SC_CHECK_MPI (mpiret);
  T8_FREE (buf);
}

void
t8_forest_pfc_message_c::fill (t8_forest_t forest, t8_procidx_t iproc)
{
  t8_procidx_t rank = forest->mpirank;
  const t8_gloidx_t *partition = t8_shmem_array_get_gloidx_array (forest->element_offsets);
  t8_gloidx_t closest_to_rank_gid = (iproc <= rank) ? partition[rank] : partition[rank + 1] - 1;
  t8_tree_t tree;
  t8_locidx_t index_in_tree;
  t8_element_t *element_closest_to_receiver;
  t8_forest_pfc_helper_global_elementid_get_stuff (forest, closest_to_rank_gid, &itree, &scheme, &tree, &index_in_tree,
                                                   &element_closest_to_receiver);
  /* if we are already the root element, we cannot be part of a split family, so we send any(the root) element and no num_siblings */

  eclass = scheme->eclass;
  if (scheme->t8_element_level (element_closest_to_receiver) == 0) {
    parent = element_closest_to_receiver;
    num_siblings = 0;
  }
  else {
    scheme->t8_element_new (1, &parent);
    allocated_parent = true;
    scheme->t8_element_parent (element_closest_to_receiver, parent);
    if (iproc > rank) {
      t8_locidx_t min_id = t8_forest_pfc_helper_min_local_sibling_closest_to_start (scheme, tree, index_in_tree);
      num_siblings = index_in_tree - min_id + 1;
      t8_debugf ("left: index_in_tree %i, min_id %i, num_siblings %i\n", index_in_tree, min_id, num_siblings);
    }
    else {
      T8_ASSERT (iproc != rank);
      t8_locidx_t max_id = t8_forest_pfc_helper_max_local_sibling_closest_to_start (scheme, tree, index_in_tree);
      num_siblings = max_id - index_in_tree + 1;
      t8_debugf ("right: index_in_tree %i, max_id %i, num_siblings %i\n", index_in_tree, max_id, num_siblings);
    }
  }
}

/**Also possible to give num_sends and array_of_procs_to_send_to, which would then be the result of determine send_range */
static void
t8_forest_pfc_send_loop (t8_forest_t forest, t8_procidx_t begin, t8_procidx_t end,
                         std::vector<sc_MPI_Request> &requests)
{
  t8_procidx_t rank = forest->mpirank;
  t8_procidx_t num_sends = requests.size ();
  const t8_gloidx_t *partition = t8_shmem_array_get_gloidx_array (forest->element_offsets);

  /* array index of send messages, starts with -1 because it is increased at the beginning of the loop whenever we actually have to send elements */
  t8_procidx_t parent_message_index = -1;
  for (t8_procidx_t iproc = begin; iproc < end; iproc++) {
    /* This process has no elements in the old partition and therefore does not need to be involved in the correction process */
    if (partition[iproc] >= partition[iproc + 1] || iproc == rank)
      continue;
    /* increase parent_index because we found a new actual process*/
    parent_message_index++;

    t8_forest_pfc_message_c message;
    message.fill (forest, iproc);
    message.pack_and_send (iproc, forest->mpicomm, &(requests[parent_message_index]));
  }
  T8_ASSERT (parent_message_index + 1 == num_sends);
}

static void
t8_forest_pfc_probe_to_get_count (t8_procidx_t iproc, sc_MPI_Comm comm, t8_procidx_t *count)
{
  sc_MPI_Status status;
  sc_MPI_Probe (iproc, T8_PFC_MESSAGE, comm, &status);
  sc_MPI_Get_count (&status, sc_MPI_PACKED, count);
}

static void
t8_forest_pfc_recv_loop (t8_forest_t forest, t8_procidx_t begin, t8_procidx_t end,
                         std::vector<t8_forest_pfc_message_c> &messages_left,
                         std::vector<t8_forest_pfc_message_c> &messages_right)
{
  t8_procidx_t rank = forest->mpirank;
  const t8_gloidx_t *partition = t8_shmem_array_get_gloidx_array (forest->element_offsets);

  /* array index of send messages, starts with -1 because it is increased at the beginning of the loop whenever we actually have to send elements */
  t8_procidx_t parent_message_index = -1;
  for (t8_procidx_t iproc = begin; iproc < end; iproc++) {
    /* This process has no elements that it could send */
    if (partition[iproc] >= partition[iproc + 1] || iproc == rank)
      continue;
    /* increase parent_index because we found a new actual process*/
    parent_message_index++;

    int message_size;
    t8_forest_pfc_probe_to_get_count (iproc, forest->mpicomm, &message_size);
    t8_debugf ("found %i bytes in message from %i\n", message_size, iproc);
    char *recv_buf = T8_ALLOC (char, message_size);

    int mpiret = sc_MPI_Recv (recv_buf, message_size, sc_MPI_PACKED, iproc, T8_PFC_MESSAGE, forest->mpicomm,
                              sc_MPI_STATUS_IGNORE);
    SC_CHECK_MPI (mpiret);

    t8_forest_pfc_message_c message;
    int position = 0;
    message.unpack (recv_buf, message_size, &position, forest->mpicomm, forest->scheme_cxx);
    T8_FREE (recv_buf);
    if (iproc < rank) {
      messages_left.push_back (std::move (message));
    }
    else {
      T8_ASSERT (iproc > rank);
      messages_right.push_back (std::move (message));
    }
  }
}

int
t8_forest_pfc_split_rank_all_to_first (t8_shmem_array_t partition_new_shmem, int rank,
                                       t8_gloidx_t min_family_element_id, t8_gloidx_t max_family_element_id)
{
  int num_ranks = t8_shmem_array_get_elem_count (partition_new_shmem);
  const t8_gloidx_t *partition_new = t8_shmem_array_get_gloidx_array (partition_new_shmem);
  const t8_gloidx_t *it = std::lower_bound (partition_new, partition_new + num_ranks, min_family_element_id);
  return (it - partition_new);
}

static void
t8_forest_pfc_corrected_local_offsets (t8_forest_t forest, t8_shmem_array_t partition_new_shmem,
                                       std::vector<t8_forest_pfc_message_c> &messages_left,
                                       std::vector<t8_forest_pfc_message_c> &messages_right, int *num_corrections,
                                       t8_gloidx_t **corrected_local_offsets)
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

  *num_corrections = next_min_local_proc - min_local_proc;
  t8_debugf ("num_corrections %i\n", *num_corrections);
  *corrected_local_offsets = T8_ALLOC (t8_gloidx_t, *num_corrections);

  for (t8_procidx_t border_rank = min_local_proc; border_rank < next_min_local_proc; border_rank++) {
    if (border_rank == 0) {
      (*corrected_local_offsets)[border_rank] = 0;
      continue;
    }
    t8_gloidx_t gtree_id;
    t8_tree_t tree;
    t8_eclass_scheme_c *scheme;
    t8_element_t *parent, *element;
    t8_locidx_t index_in_tree;
    t8_forest_pfc_helper_global_elementid_get_stuff (forest, partition_new[border_rank], &gtree_id, &scheme, &tree,
                                                     &index_in_tree, &element);
    scheme->t8_element_new (1, &parent);
    scheme->t8_element_parent (element, parent);
    t8_debugf ("borde_rank:%i\n", border_rank);
    t8_debugf ("parent\n");
    scheme->t8_element_debug_print (parent);

    t8_gloidx_t min_family_element_id = t8_forest_get_first_local_element_id (forest)
                                        + t8_forest_pfc_helper_min_local_sibling_closest_to_start (
                                          scheme, tree, index_in_tree); /* min quadrant id with same parent */
    t8_debugf ("min_family_element_id:%i\n", min_family_element_id);
    t8_procidx_t num_messages = (t8_procidx_t) messages_left.size ();
    t8_debugf ("num_messages: %i\n", num_messages);
    for (int iproc = 0; iproc < num_messages; iproc++) {
      t8_debugf ("message from iproc %i\n", iproc);
      if (messages_left[iproc].itree == gtree_id) {
        const t8_element_t *iproc_parent = messages_left[iproc].get_parent ();
        t8_debugf ("iproc left parent\n");
        scheme->t8_element_debug_print (iproc_parent);
        if (scheme->t8_element_equal (parent, messages_left[iproc].get_parent ())) {
          min_family_element_id -= messages_left[iproc].num_siblings;
        }
      }
    }

    t8_gloidx_t max_family_element_id = t8_forest_get_first_local_element_id (forest)
                                        + t8_forest_pfc_helper_max_local_sibling_closest_to_start (
                                          scheme, tree, index_in_tree); /* min quadrant id with same parent */
    t8_debugf ("max_family_element_id:%i\n", max_family_element_id);
    num_messages = (t8_procidx_t) messages_right.size ();
    t8_debugf ("num_messages: %i\n", num_messages);
    for (int iproc = 0; iproc < num_messages; iproc++) {
      t8_debugf ("message from iproc %i\n", iproc);
      if (messages_right[iproc].itree == gtree_id) {
        const t8_element_t *iproc_parent = messages_right[iproc].get_parent ();
        t8_debugf ("iproc right parent\n");
        scheme->t8_element_debug_print (iproc_parent);
        if (scheme->t8_element_equal (parent, iproc_parent)) {
          max_family_element_id += messages_right[iproc].num_siblings;
          t8_debugf ("adapted max_family_element_id %i\n", max_family_element_id);
        }
      }
    }
    int num_children = scheme->t8_element_num_children (parent);
    scheme->t8_element_destroy (1, &parent);

    t8_debugf ("found max:%i, found min:%i, numchildren: %i\n", max_family_element_id, min_family_element_id,
               num_children);
    if (max_family_element_id - min_family_element_id + 1 == num_children) {
      /* compute correction */
      t8_procidx_t rank = forest->mpirank;
      t8_procidx_t split_rank = t8_forest_pfc_split_rank_all_to_first (forest->element_offsets, rank,
                                                                       min_family_element_id, max_family_element_id);
      t8_debugf ("splitrank: %i\n", split_rank);
      (*corrected_local_offsets)[border_rank - min_local_proc]
        = (rank <= split_rank) ? min_family_element_id : max_family_element_id + 1;
      t8_debugf ("corrected_local_offsets[%i]: %i\n", border_rank - min_local_proc,
                 (*corrected_local_offsets)[border_rank - min_local_proc]);
      /* possible TODO: update multiple corrected_local_offsets*/
    }
    else {
      (*corrected_local_offsets)[border_rank - min_local_proc] = partition_new[border_rank];
    }
  }
}

T8_EXTERN_C_BEGIN ();
void
t8_forest_pfc_correction_offsets (t8_forest_t forest)
{
  /** Determine procs that we need to send to (i.e. they have overlapping elements or close to)*/
  t8_procidx_t num_procs = forest->mpisize;
  t8_procidx_t rank = forest->mpirank;
  t8_procidx_t num_nonzero, begin, end;
  const t8_shmem_array_t partition_old = forest->set_from->element_offsets;
  t8_forest_pfc_determine_send_recv_range (partition_old, rank, &begin, &end, &num_nonzero);

  std::vector<t8_forest_pfc_message_c> messages_left;
  std::vector<t8_forest_pfc_message_c> messages_right;
  if (num_nonzero) {
    std::vector<sc_MPI_Request> requests (num_nonzero);
    t8_forest_pfc_send_loop (forest->set_from, begin, end, requests);

    /** Determine procs that we need to receive from */
    t8_procidx_t num_recvs, lowest_recv, highest_recv;
    t8_forest_pfc_determine_send_recv_range_old (partition_old, rank, &lowest_recv, &highest_recv, &num_recvs);

    T8_ASSERT (begin == lowest_recv);
    T8_ASSERT (end == highest_recv + 1);
    T8_ASSERT (num_recvs == num_nonzero);
    /** receive and unpack parent, tree and numsiblings(left/right)*/

    t8_forest_pfc_recv_loop (forest->set_from, begin, end, messages_left, messages_right);

    int mpiret = sc_MPI_Waitall (num_nonzero, requests.data (), sc_MPI_STATUSES_IGNORE);
    SC_CHECK_MPI (mpiret);
  }

  int num_local_corrections;
  t8_gloidx_t *corrected_local_offsets;
  t8_forest_pfc_corrected_local_offsets (forest->set_from, forest->element_offsets, messages_left, messages_right,
                                         &num_local_corrections, &corrected_local_offsets);

  /** Allgather corrected partition */
  t8_shmem_array_t global_correction_counts;
  t8_shmem_array_init (&global_correction_counts, sizeof (t8_gloidx_t), num_procs + 1, forest->mpicomm);

  t8_shmem_array_allgather (&num_local_corrections, 1, sc_MPI_INT, global_correction_counts, num_procs, sc_MPI_INT);
  t8_shmem_array_allgatherv (corrected_local_offsets, num_local_corrections, T8_MPI_GLOIDX, forest->element_offsets,
                             T8_MPI_GLOIDX, forest->mpicomm);
  t8_shmem_array_destroy (&global_correction_counts);
  T8_FREE (corrected_local_offsets);
}
T8_EXTERN_C_END ();
