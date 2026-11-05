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

/** \file t8_shmem.cxx
 * Implementation of \ref t8_shmem_array and of the C interface in t8_shmem.h.
 *
 * Every communicator gets a \ref t8_shmem_topology that is computed once and cached as an MPI attribute.
 * In node-shared mode, intranode rank 0 of every node owns an MPI-3 shared memory window that all processes of
 * the node read. A fill gathers the contributions of a node on its root, the roots exchange the node blocks
 * directly into their windows, and a final intranode barrier publishes the data.
 */

#include <t8_data/t8_shmem.hxx>

#include <climits>
#include <cstdlib>
#include <cstring>
#include <string>
#include <vector>

#if T8_ENABLE_MPI && defined(SC_ENABLE_MPICOMMSHARED) && defined(SC_ENABLE_MPIWINSHARED)
/** Defined if MPI-3 shared memory windows can be used. */
#define T8_SHMEM_WITH_WINDOWS 1
#endif

/** The processes of a communicator grouped into shared memory nodes. All values are identical on all processes
 * of a node, except for the ranks. */
struct t8_shmem_topology
{
  int mpirank = 0;               /**< The rank in the communicator. */
  int mpisize = 1;               /**< The size of the communicator. */
  int intrarank = 0;             /**< The rank in the node. */
  int intrasize = 1;             /**< The number of processes of the node. */
  bool can_share = false;        /**< True if the processes of a node can share memory, false otherwise. */
  bool use_node_shared = false;  /**< True if new arrays are node-shared, false if they are replicated. */
  bool equal_node_sizes = true;  /**< True if all nodes have the same number of processes, false otherwise. */
  std::vector<int> node_sizes;   /**< On node roots: the number of processes of each node. Empty otherwise. */
  std::vector<int> node_offsets; /**< On node roots: the rank of the first process of each node. Empty otherwise. */
#if T8_SHMEM_WITH_WINDOWS
  MPI_Comm intranode = MPI_COMM_NULL; /**< The processes of this node. MPI_COMM_NULL if \a can_share is false. */
  MPI_Comm roots = MPI_COMM_NULL;     /**< The roots of all nodes. MPI_COMM_NULL on non-roots. */

  /** Free the communicators, unless MPI is already finalized. */
  ~t8_shmem_topology ()
  {
    int finalized = 0;
    MPI_Finalized (&finalized);
    if (finalized) {
      return;
    }
    if (roots != MPI_COMM_NULL) {
      MPI_Comm_free (&roots);
    }
    if (intranode != MPI_COMM_NULL) {
      MPI_Comm_free (&intranode);
    }
  }
#endif
};

/** Check that a number of bytes fits into an int as required by MPI.
 * \param [in] bytes  The number of bytes.
 * \return \a bytes as int. */
static int
t8_shmem_mpi_count (const size_t bytes)
{
  SC_CHECK_ABORT (bytes <= (size_t) INT_MAX, "Shared memory array too large for MPI communication.");
  return static_cast<int> (bytes);
}

#if T8_SHMEM_WITH_WINDOWS

/** The keyval of the topology attribute. */
static int t8_shmem_topology_keyval = MPI_KEYVAL_INVALID;

/** Delete callback of the topology attribute. Only drops the reference of the attribute, since arrays may still
 * hold the topology. */
static int
t8_shmem_topology_delete ([[maybe_unused]] MPI_Comm comm, [[maybe_unused]] int keyval, void *attribute,
                          [[maybe_unused]] void *extra_state)
{
  delete static_cast<std::shared_ptr<t8_shmem_topology> *> (attribute);
  return MPI_SUCCESS;
}

/** Return the color of the intranode communicator that simulates multiple nodes in tests.
 * Reads the environment variable T8_SHMEM_SIMULATE_NODES in debugging mode. Accepted values are
 * "block:<k>" (k consecutive ranks per node), "cyclic:<k>" (k nodes, round robin) and "replicated".
 * \param [in] rank              The rank in the communicator.
 * \param [out] force_replicated True if the value is "replicated", false otherwise.
 * \return The color, or 0 if nothing is simulated. */
static int
t8_shmem_simulated_node_color ([[maybe_unused]] const int rank, bool *force_replicated)
{
  *force_replicated = false;
#if T8_ENABLE_DEBUG
  const char *value = std::getenv ("T8_SHMEM_SIMULATE_NODES");
  if (value == nullptr || *value == '\0') {
    return 0;
  }
  const std::string layout (value);
  if (layout == "replicated") {
    *force_replicated = true;
    return 0;
  }
  const size_t colon = layout.find (':');
  SC_CHECK_ABORTF (colon != std::string::npos, "Invalid T8_SHMEM_SIMULATE_NODES value %s.", value);
  const int k = std::atoi (layout.c_str () + colon + 1);
  SC_CHECK_ABORTF (k > 0, "Invalid T8_SHMEM_SIMULATE_NODES value %s.", value);
  if (layout.compare (0, colon, "block") == 0) {
    return rank / k;
  }
  SC_CHECK_ABORTF (layout.compare (0, colon, "cyclic") == 0, "Invalid T8_SHMEM_SIMULATE_NODES value %s.", value);
  return rank % k;
#else
  return 0;
#endif
}

/** Collectively compute the topology of a communicator.
 * \param [in] comm  The communicator.
 * \return The topology. */
static std::shared_ptr<t8_shmem_topology>
t8_shmem_topology_create (sc_MPI_Comm comm)
{
  auto topology = std::make_shared<t8_shmem_topology> ();
  int mpiret = MPI_Comm_rank (comm, &topology->mpirank);
  SC_CHECK_MPI (mpiret);
  mpiret = MPI_Comm_size (comm, &topology->mpisize);
  SC_CHECK_MPI (mpiret);
  if (topology->mpisize == 1) {
    return topology;
  }
  const int rank = topology->mpirank;

  /* Group the processes into nodes. */
  bool force_replicated;
  const int color = t8_shmem_simulated_node_color (rank, &force_replicated);
  MPI_Comm shared;
  mpiret = MPI_Comm_split_type (comm, MPI_COMM_TYPE_SHARED, rank, MPI_INFO_NULL, &shared);
  SC_CHECK_MPI (mpiret);
  mpiret = MPI_Comm_split (shared, color, rank, &topology->intranode);
  SC_CHECK_MPI (mpiret);
  mpiret = MPI_Comm_free (&shared);
  SC_CHECK_MPI (mpiret);
  mpiret = MPI_Comm_rank (topology->intranode, &topology->intrarank);
  SC_CHECK_MPI (mpiret);
  mpiret = MPI_Comm_size (topology->intranode, &topology->intrasize);
  SC_CHECK_MPI (mpiret);

  /* The ranks of a node must be contiguous, since a node block is exchanged as one piece. */
  int rank_bounds[2] = { -rank, rank };
  mpiret = MPI_Allreduce (MPI_IN_PLACE, rank_bounds, 2, MPI_INT, MPI_MAX, topology->intranode);
  SC_CHECK_MPI (mpiret);
  const bool contiguous = rank_bounds[1] + rank_bounds[0] + 1 == topology->intrasize;

  /* Loads and stores to the window are only defined in the unified memory model. */
  void *probe_base;
  MPI_Win probe;
  mpiret = MPI_Win_allocate_shared (0, 1, MPI_INFO_NULL, topology->intranode, &probe_base, &probe);
  SC_CHECK_MPI (mpiret);
  int *memory_model;
  int flag;
  mpiret = MPI_Win_get_attr (probe, MPI_WIN_MODEL, &memory_model, &flag);
  SC_CHECK_MPI (mpiret);
  const bool unified = flag && *memory_model == MPI_WIN_UNIFIED;
  mpiret = MPI_Win_free (&probe);
  SC_CHECK_MPI (mpiret);

  /* Agree on the decision. */
  int agreed[4] = { topology->intrasize, -topology->intrasize, !contiguous, !unified || force_replicated };
  mpiret = MPI_Allreduce (MPI_IN_PLACE, agreed, 4, MPI_INT, MPI_MAX, comm);
  SC_CHECK_MPI (mpiret);
  topology->equal_node_sizes = agreed[0] == -agreed[1];
  topology->can_share = agreed[0] > 1 && !agreed[2] && !agreed[3];
  if (agreed[2]) {
    t8_global_productionf ("The processes of a shared memory node are not contiguous. Shared memory arrays are "
                           "replicated on every process.\n");
  }
  if (!topology->can_share) {
    mpiret = MPI_Comm_free (&topology->intranode);
    SC_CHECK_MPI (mpiret);
    topology->intrarank = 0;
    topology->intrasize = 1;
    return topology;
  }
  topology->use_node_shared = true;

  /* The node roots exchange the node blocks. */
  mpiret = MPI_Comm_split (comm, topology->intrarank == 0 ? 0 : MPI_UNDEFINED, rank, &topology->roots);
  SC_CHECK_MPI (mpiret);
  if (topology->roots != MPI_COMM_NULL) {
    int num_nodes;
    mpiret = MPI_Comm_size (topology->roots, &num_nodes);
    SC_CHECK_MPI (mpiret);
    topology->node_sizes.resize (num_nodes);
    topology->node_offsets.resize (num_nodes);
    mpiret
      = MPI_Allgather (&topology->intrasize, 1, MPI_INT, topology->node_sizes.data (), 1, MPI_INT, topology->roots);
    SC_CHECK_MPI (mpiret);
    mpiret
      = MPI_Allgather (&topology->mpirank, 1, MPI_INT, topology->node_offsets.data (), 1, MPI_INT, topology->roots);
    SC_CHECK_MPI (mpiret);
  }
  return topology;
}
#endif /* T8_SHMEM_WITH_WINDOWS */

/** Return the topology of a communicator and compute it if necessary. Collective if it is not computed yet.
 * \param [in] comm  The communicator.
 * \return The topology. */
static std::shared_ptr<t8_shmem_topology>
t8_shmem_get_topology (sc_MPI_Comm comm)
{
  SC_CHECK_ABORT (comm != sc_MPI_COMM_NULL, "Shared memory needs a valid communicator.");
#if T8_SHMEM_WITH_WINDOWS
  int mpiret;
  if (t8_shmem_topology_keyval == MPI_KEYVAL_INVALID) {
    mpiret
      = MPI_Comm_create_keyval (MPI_COMM_NULL_COPY_FN, t8_shmem_topology_delete, &t8_shmem_topology_keyval, nullptr);
    SC_CHECK_MPI (mpiret);
  }
  std::shared_ptr<t8_shmem_topology> *attribute;
  int flag;
  mpiret = MPI_Comm_get_attr (comm, t8_shmem_topology_keyval, &attribute, &flag);
  SC_CHECK_MPI (mpiret);
  if (flag) {
    return *attribute;
  }
  std::shared_ptr<t8_shmem_topology> topology = t8_shmem_topology_create (comm);
  mpiret = MPI_Comm_set_attr (comm, t8_shmem_topology_keyval, new std::shared_ptr<t8_shmem_topology> (topology));
  SC_CHECK_MPI (mpiret);
  return topology;
#else
  /* Without shared memory windows every process has its own copy. */
  auto topology = std::make_shared<t8_shmem_topology> ();
  int mpiret = sc_MPI_Comm_rank (comm, &topology->mpirank);
  SC_CHECK_MPI (mpiret);
  mpiret = sc_MPI_Comm_size (comm, &topology->mpisize);
  SC_CHECK_MPI (mpiret);
  return topology;
#endif
}

t8_shmem_storage::t8_shmem_storage (sc_MPI_Comm comm_in, const size_t elem_size_in, const size_t count_in)
  : elem_size (elem_size_in), count (count_in), comm (comm_in), topology (t8_shmem_get_topology (comm_in))
{
  T8_ASSERT (elem_size > 0);
#if T8_ENABLE_DEBUG && T8_ENABLE_MPI
  /* All processes must allocate the same number of elements. */
  long long count_bounds[2] = { -(long long) count, (long long) count };
  const int mpiret = sc_MPI_Allreduce (MPI_IN_PLACE, count_bounds, 2, sc_MPI_LONG_LONG_INT, sc_MPI_MAX, comm);
  SC_CHECK_MPI (mpiret);
  T8_ASSERT (count_bounds[0] == -count_bounds[1]);
#endif
  const size_t bytes = elem_size * count;
#if T8_SHMEM_WITH_WINDOWS
  if (topology->use_node_shared && count > 0) {
    node_shared = true;
    writer = topology->intrarank == 0;
    void *local_base;
    int mpiret = MPI_Win_allocate_shared (writer ? (MPI_Aint) bytes : 0, 1, MPI_INFO_NULL, topology->intranode,
                                          &local_base, &window);
    SC_CHECK_MPI (mpiret);
    MPI_Aint window_size;
    int disp_unit;
    mpiret = MPI_Win_shared_query (window, 0, &window_size, &disp_unit, &data_ptr);
    SC_CHECK_MPI (mpiret);
    T8_ASSERT ((size_t) window_size >= bytes);
    mpiret = MPI_Win_lock_all (MPI_MODE_NOCHECK, window);
    SC_CHECK_MPI (mpiret);
    if (writer) {
      /* Value-initialize the content. The next fill publishes it together with the filled values. */
      std::memset (data_ptr, 0, bytes);
    }
    return;
  }
#endif
  writer = true;
  if (bytes > 0) {
    data_ptr = reinterpret_cast<std::byte *> (T8_ALLOC_ZERO (char, bytes));
  }
}

void
t8_shmem_storage::free_storage () noexcept
{
  T8_ASSERT (!writing);
#if T8_SHMEM_WITH_WINDOWS
  if (node_shared) {
    MPI_Win_unlock_all (window);
    MPI_Win_free (&window);
  }
#endif
  if (!node_shared && data_ptr != nullptr) {
    T8_FREE (data_ptr);
  }
  data_ptr = nullptr;
  count = 0;
  node_shared = false;
  published = false;
  topology.reset ();
}

t8_shmem_storage::~t8_shmem_storage ()
{
  free_storage ();
}

t8_shmem_storage::t8_shmem_storage (t8_shmem_storage &&other) noexcept
  : data_ptr (other.data_ptr), elem_size (other.elem_size), count (other.count), comm (other.comm),
    topology (std::move (other.topology)), node_shared (other.node_shared), writer (other.writer),
    writing (other.writing), published (other.published)
#if T8_ENABLE_MPI
    ,
    window (other.window)
#endif
{
  other.data_ptr = nullptr;
  other.count = 0;
  other.node_shared = false;
  other.writing = false;
  other.published = false;
#if T8_ENABLE_MPI
  other.window = MPI_WIN_NULL;
#endif
}

t8_shmem_storage &
t8_shmem_storage::operator= (t8_shmem_storage &&other) noexcept
{
  if (this != &other) {
    free_storage ();
    std::swap (data_ptr, other.data_ptr);
    std::swap (elem_size, other.elem_size);
    std::swap (count, other.count);
    std::swap (comm, other.comm);
    std::swap (topology, other.topology);
    std::swap (node_shared, other.node_shared);
    std::swap (writer, other.writer);
    std::swap (writing, other.writing);
    std::swap (published, other.published);
#if T8_ENABLE_MPI
    std::swap (window, other.window);
#endif
  }
  return *this;
}

bool
t8_shmem_storage::start_writing ()
{
  SC_CHECK_ABORT (!writing, "A write epoch of this shared memory array is already active.");
#if T8_SHMEM_WITH_WINDOWS
  if (node_shared && published) {
    /* Wait until all processes of the node have finished reading the old content. */
    int mpiret = MPI_Win_sync (window);
    SC_CHECK_MPI (mpiret);
    mpiret = MPI_Barrier (topology->intranode);
    SC_CHECK_MPI (mpiret);
  }
#endif
  writing = true;
  return writer;
}

void
t8_shmem_storage::end_writing ()
{
  SC_CHECK_ABORT (writing, "End writing to a shared memory array is only possible after start_writing.");
#if T8_SHMEM_WITH_WINDOWS
  if (node_shared) {
    /* Make the written data visible to all processes of the node. */
    int mpiret = MPI_Win_sync (window);
    SC_CHECK_MPI (mpiret);
    mpiret = MPI_Barrier (topology->intranode);
    SC_CHECK_MPI (mpiret);
    mpiret = MPI_Win_sync (window);
    SC_CHECK_MPI (mpiret);
  }
#endif
  writing = false;
  published = true;
}

bool
t8_shmem_storage::gather_fill (const void *local_value, const size_t offset)
{
  SC_CHECK_ABORT (!writing, "A write epoch of this shared memory array is already active.");
  const size_t mpisize = topology->mpisize;
  SC_CHECK_ABORT (offset + mpisize <= count, "Shared memory array too small for an allgather.");
  const int bytes = t8_shmem_mpi_count (elem_size);
  t8_shmem_mpi_count (elem_size * mpisize);
  writing = true;
  std::byte *const dest = data_ptr + offset * elem_size;
#if T8_SHMEM_WITH_WINDOWS
  if (node_shared) {
    /* The root receives only after all processes of the node have entered the gather, i.e. finished reading. */
    const int intrasize = topology->intrasize;
    std::vector<std::byte> node_block (writer ? elem_size * intrasize : 0);
    int mpiret = MPI_Gather (local_value, bytes, MPI_BYTE, node_block.data (), bytes, MPI_BYTE, 0, topology->intranode);
    SC_CHECK_MPI (mpiret);
    if (writer) {
      if (topology->equal_node_sizes) {
        mpiret = MPI_Allgather (node_block.data (), bytes * intrasize, MPI_BYTE, dest, bytes * intrasize, MPI_BYTE,
                                topology->roots);
      }
      else {
        std::vector<int> counts (topology->node_sizes);
        std::vector<int> displs (topology->node_offsets);
        for (size_t inode = 0; inode < counts.size (); ++inode) {
          counts[inode] *= bytes;
          displs[inode] *= bytes;
        }
        mpiret = MPI_Allgatherv (node_block.data (), bytes * intrasize, MPI_BYTE, dest, counts.data (), displs.data (),
                                 MPI_BYTE, topology->roots);
      }
      SC_CHECK_MPI (mpiret);
    }
    return writer;
  }
#endif
  const int mpiret
    = sc_MPI_Allgather (const_cast<void *> (local_value), bytes, sc_MPI_BYTE, dest, bytes, sc_MPI_BYTE, comm);
  SC_CHECK_MPI (mpiret);
  return writer;
}

/** Compute the displacements of a gatherv from the counts.
 * \param [in] counts  The number of bytes of every process.
 * \param [in] total_limit  The maximum allowed total number of bytes.
 * \return The displacements and the total number of bytes. */
static std::pair<std::vector<int>, int>
t8_shmem_displacements (const std::vector<int> &counts, const size_t total_limit)
{
  std::vector<int> displs (counts.size ());
  size_t total = 0;
  for (size_t iproc = 0; iproc < counts.size (); ++iproc) {
    displs[iproc] = t8_shmem_mpi_count (total);
    total += counts[iproc];
  }
  SC_CHECK_ABORT (total <= total_limit, "Shared memory array too small for an allgatherv.");
  return { displs, t8_shmem_mpi_count (total) };
}

bool
t8_shmem_storage::gatherv_fill (const void *local_values, const size_t local_count)
{
  SC_CHECK_ABORT (!writing, "A write epoch of this shared memory array is already active.");
  int bytes = t8_shmem_mpi_count (elem_size * local_count);
  writing = true;
  int mpiret;
#if T8_SHMEM_WITH_WINDOWS
  if (node_shared) {
    /* Gather the counts of the node. The root receives only after all processes of the node have entered, i.e.
     * finished reading, even if some processes send no data. */
    std::vector<int> intra_counts (writer ? topology->intrasize : 0);
    mpiret = MPI_Gather (&bytes, 1, MPI_INT, intra_counts.data (), 1, MPI_INT, 0, topology->intranode);
    SC_CHECK_MPI (mpiret);
    const auto [intra_displs, node_bytes] = t8_shmem_displacements (intra_counts, elem_size * count);
    std::vector<std::byte> node_block (node_bytes);
    mpiret = MPI_Gatherv (local_values, bytes, MPI_BYTE, node_block.data (), intra_counts.data (), intra_displs.data (),
                          MPI_BYTE, 0, topology->intranode);
    SC_CHECK_MPI (mpiret);
    if (writer) {
      std::vector<int> node_counts (topology->node_sizes.size ());
      mpiret = MPI_Allgather (&node_bytes, 1, MPI_INT, node_counts.data (), 1, MPI_INT, topology->roots);
      SC_CHECK_MPI (mpiret);
      const auto [node_displs, total_bytes] = t8_shmem_displacements (node_counts, elem_size * count);
      mpiret = MPI_Allgatherv (node_block.data (), node_bytes, MPI_BYTE, data_ptr, node_counts.data (),
                               node_displs.data (), MPI_BYTE, topology->roots);
      SC_CHECK_MPI (mpiret);
    }
    return writer;
  }
#endif
  std::vector<int> counts (topology->mpisize);
  mpiret = sc_MPI_Allgather (&bytes, 1, sc_MPI_INT, counts.data (), 1, sc_MPI_INT, comm);
  SC_CHECK_MPI (mpiret);
  const auto [displs, total_bytes] = t8_shmem_displacements (counts, elem_size * count);
  mpiret = sc_MPI_Allgatherv (const_cast<void *> (local_values), bytes, sc_MPI_BYTE, data_ptr, counts.data (),
                              const_cast<int *> (displs.data ()), sc_MPI_BYTE, comm);
  SC_CHECK_MPI (mpiret);
  return writer;
}

/* The C interface. The C handle is a t8_shmem_array<t8_gloidx_t>. */

int
t8_shmem_init (sc_MPI_Comm comm)
{
  const std::shared_ptr<t8_shmem_topology> topology = t8_shmem_get_topology (comm);
  return topology->can_share ? topology->intrasize : 1;
}

void
t8_shmem_finalize ([[maybe_unused]] sc_MPI_Comm comm)
{
#if T8_SHMEM_WITH_WINDOWS
  if (t8_shmem_topology_keyval != MPI_KEYVAL_INVALID && comm != sc_MPI_COMM_NULL) {
    const int mpiret = MPI_Comm_delete_attr (comm, t8_shmem_topology_keyval);
    SC_CHECK_MPI (mpiret);
  }
#endif
}

void
t8_shmem_set_type (sc_MPI_Comm comm, sc_shmem_type_t type)
{
  const std::shared_ptr<t8_shmem_topology> topology = t8_shmem_get_topology (comm);
  const bool replicated = type == SC_SHMEM_BASIC || type == SC_SHMEM_PRESCAN;
  topology->use_node_shared = topology->can_share && !replicated;
}

#if T8_ENABLE_DEBUG
int
t8_shmem_array_is_initialized (const t8_shmem_array_t array)
{
  return array != nullptr && array->get_comm () != sc_MPI_COMM_NULL;
}
#endif

void
t8_shmem_array_init (t8_shmem_array_t *parray, size_t elem_size, size_t elem_count, sc_MPI_Comm comm)
{
  T8_ASSERT (parray != nullptr);
  SC_CHECK_ABORT (elem_size == sizeof (t8_gloidx_t),
                  "The C interface of t8_shmem_array only supports elements of type t8_gloidx_t.");
  *parray = new t8_shmem_array<t8_gloidx_t> (comm, elem_count);
}

int
t8_shmem_array_start_writing (t8_shmem_array_t array)
{
  T8_ASSERT (t8_shmem_array_is_initialized (array));
  return array->start_writing ();
}

void
t8_shmem_array_end_writing (t8_shmem_array_t array)
{
  T8_ASSERT (t8_shmem_array_is_initialized (array));
  array->end_writing ();
}

void
t8_shmem_array_set_gloidx (t8_shmem_array_t array, int index, t8_gloidx_t value)
{
  T8_ASSERT (t8_shmem_array_is_initialized (array));
  T8_ASSERT (0 <= index && (size_t) index < array->size ());
  array->data_for_writing ()[index] = value;
}

void
t8_shmem_array_copy (t8_shmem_array_t dest, t8_shmem_array_t source)
{
  T8_ASSERT (t8_shmem_array_is_initialized (dest));
  T8_ASSERT (t8_shmem_array_is_initialized (source));
  SC_CHECK_ABORT (dest->size () == source->size (), "Try to copy shared memory arrays of different element counts.");
  SC_CHECK_ABORT (dest->get_comm () == source->get_comm (),
                  "Try to copy shared memory arrays with different communicators.");
  dest->write ([source] (std::span<t8_gloidx_t> data) { std::ranges::copy (*source, data.begin ()); });
}

void
t8_shmem_array_allgather (const void *sendbuf, int sendcount, sc_MPI_Datatype sendtype, t8_shmem_array_t recvarray,
                          int recvcount, sc_MPI_Datatype recvtype)
{
  T8_ASSERT (t8_shmem_array_is_initialized (recvarray));
  SC_CHECK_ABORT (sendcount == 1 && recvcount == 1 && sc_mpi_sizeof (sendtype) == sizeof (t8_gloidx_t)
                    && sc_mpi_sizeof (recvtype) == sizeof (t8_gloidx_t),
                  "The C interface of t8_shmem_array only supports allgather of one t8_gloidx_t.");
  recvarray->allgather (*static_cast<const t8_gloidx_t *> (sendbuf));
}

void
t8_shmem_array_allgatherv (void *sendbuf, const int sendcount, sc_MPI_Datatype sendtype, t8_shmem_array_t recvarray,
                           sc_MPI_Datatype recvtype, [[maybe_unused]] sc_MPI_Comm comm)
{
  T8_ASSERT (t8_shmem_array_is_initialized (recvarray));
  T8_ASSERT (comm == recvarray->get_comm ());
  SC_CHECK_ABORT (sendcount >= 0 && sc_mpi_sizeof (sendtype) == sizeof (t8_gloidx_t)
                    && sc_mpi_sizeof (recvtype) == sizeof (t8_gloidx_t),
                  "The C interface of t8_shmem_array only supports allgatherv of t8_gloidx_t.");
  recvarray->allgatherv (std::span<const t8_gloidx_t> (static_cast<const t8_gloidx_t *> (sendbuf), sendcount));
}

void
t8_shmem_array_prefix (const void *sendbuf, t8_shmem_array_t recvarray, const int count, sc_MPI_Datatype type,
                       sc_MPI_Op op, [[maybe_unused]] sc_MPI_Comm comm)
{
  T8_ASSERT (t8_shmem_array_is_initialized (recvarray));
  T8_ASSERT (comm == recvarray->get_comm ());
  SC_CHECK_ABORT (count == 1 && sc_mpi_sizeof (type) == sizeof (t8_gloidx_t) && op == sc_MPI_SUM,
                  "The C interface of t8_shmem_array only supports the prefix sum of one t8_gloidx_t.");
  recvarray->prefix_sum (*static_cast<const t8_gloidx_t *> (sendbuf));
}

sc_MPI_Comm
t8_shmem_array_get_comm (t8_shmem_array_t array)
{
  T8_ASSERT (t8_shmem_array_is_initialized (array));
  return array->get_comm ();
}

size_t
t8_shmem_array_get_elem_size ([[maybe_unused]] t8_shmem_array_t array)
{
  T8_ASSERT (t8_shmem_array_is_initialized (array));
  return sizeof (t8_gloidx_t);
}

size_t
t8_shmem_array_get_elem_count (t8_shmem_array_t array)
{
  T8_ASSERT (t8_shmem_array_is_initialized (array));
  return array->size ();
}

const t8_gloidx_t *
t8_shmem_array_get_gloidx_array (t8_shmem_array_t array)
{
  T8_ASSERT (t8_shmem_array_is_initialized (array));
  return array->data ();
}

t8_gloidx_t *
t8_shmem_array_get_gloidx_array_for_writing (t8_shmem_array_t array)
{
  T8_ASSERT (t8_shmem_array_is_initialized (array));
  return array->data_for_writing ().data ();
}

t8_gloidx_t
t8_shmem_array_get_gloidx (t8_shmem_array_t array, int index)
{
  T8_ASSERT (t8_shmem_array_is_initialized (array));
  T8_ASSERT (0 <= index);
  return (*array)[index];
}

const void *
t8_shmem_array_get_array (t8_shmem_array_t array)
{
  T8_ASSERT (t8_shmem_array_is_initialized (array));
  return array->data ();
}

const void *
t8_shmem_array_index (t8_shmem_array_t array, size_t index)
{
  T8_ASSERT (t8_shmem_array_is_initialized (array));
  return &(*array)[index];
}

void *
t8_shmem_array_index_for_writing (t8_shmem_array_t array, size_t index)
{
  T8_ASSERT (t8_shmem_array_is_initialized (array));
  T8_ASSERT (index < array->size ());
  return array->data_for_writing ().data () + index;
}

int
t8_shmem_array_is_equal (t8_shmem_array_t array_a, t8_shmem_array_t array_b)
{
  if (array_a == array_b) {
    return 1;
  }
  if (array_a == nullptr || array_b == nullptr) {
    return 0;
  }
  return array_a->get_comm () == array_b->get_comm () && *array_a == *array_b;
}

void
t8_shmem_array_destroy (t8_shmem_array_t *parray)
{
  T8_ASSERT (parray != nullptr);
  T8_ASSERT (t8_shmem_array_is_initialized (*parray));
  delete *parray;
  *parray = nullptr;
}

int
t8_shmem_array_binary_search (t8_shmem_array_t array, const t8_gloidx_t value, const int size,
                              int (*compare) (t8_shmem_array_t, const int, const t8_gloidx_t))
{
  int low = 0;
  int high = size;

  while (low <= high) {
    const int mid = low + (high - low) / 2;
    const int eval = compare (array, mid, value);
    if (eval == 0) {
      /* mid points to a value that equals \a value */
      return mid;
    }
    else if (eval < 0) {
      /* \a value is smaller than the value at mid */
      high = mid - 1;
    }
    else {
      /* \a value is greater than the value at mid */
      low = mid + 1;
    }
  }
  return -1; /* Not found */
}
