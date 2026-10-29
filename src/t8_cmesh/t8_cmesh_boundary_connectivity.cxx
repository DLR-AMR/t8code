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

#include <t8_cmesh/t8_cmesh.h>
#include <t8_cmesh/t8_cmesh_boundary_connectivity.hxx>
#include <t8_cmesh/t8_cmesh_internal/t8_cmesh_types.h>
#include <t8_cmesh/t8_cmesh_internal/t8_cmesh_offset.h>
#include <t8_cmesh/t8_cmesh_vertex_connectivity/t8_cmesh_vertex_connectivity.hxx>
#include <t8_data/t8_shmem.h>
#include <algorithm>
#include <limits>
#include <map>
#include <set>
#include <unordered_map>

/** The label key of a vertex or edge is 12 * global tree id + vertex or edge, 12 is the maximum number of both. */
constexpr t8_gloidx_t T8_BOUNDARY_KEY_STRIDE = 12;

/** Union-find with the parity of each element relative to the representative of its set. */
struct t8_parity_union_find
{
  explicit t8_parity_union_find (const size_t size): parent (size), parity (size, 0)
  {
    for (size_t i = 0; i < size; i++) {
      parent[i] = i;
    }
  }

  /** The representative of the set of \a i, and in \a i_parity the parity of \a i relative to it. */
  size_t
  find (const size_t i, int *i_parity)
  {
    if (parent[i] == i) {
      *i_parity = 0;
      return i;
    }
    int parent_parity;
    const size_t root = find (parent[i], &parent_parity);
    parity[i] ^= parent_parity;
    parent[i] = root;
    *i_parity = parity[i];
    return root;
  }

  /** Join the sets of \a i and \a j, such that the parity of \a j relative to \a i is \a relative_parity.
   * If both are already in the same set, the existing parity is kept. */
  void
  join (const size_t i, const size_t j, const int relative_parity)
  {
    int i_parity, j_parity;
    const size_t i_root = find (i, &i_parity);
    const size_t j_root = find (j, &j_parity);
    if (i_root != j_root) {
      parent[j_root] = i_root;
      parity[j_root] = i_parity ^ j_parity ^ relative_parity;
    }
  }

  std::vector<size_t> parent;
  std::vector<int> parity;
};

/** The eclass of a local tree or ghost. */
static t8_eclass_t
t8_cmesh_local_or_ghost_class (const t8_cmesh_t cmesh, const t8_locidx_t ltreeid)
{
  const t8_locidx_t num_local_trees = t8_cmesh_get_num_local_trees (cmesh);
  return ltreeid < num_local_trees ? t8_cmesh_get_tree_class (cmesh, ltreeid)
                                   : t8_cmesh_get_ghost_class (cmesh, ltreeid - num_local_trees);
}

/** The number of vertices (bdy_dim 0) or edges (bdy_dim 1, only for 3D trees) of an eclass. */
static int
t8_cmesh_num_tree_boundaries (const t8_eclass_t eclass, const int bdy_dim)
{
  if (bdy_dim == 0) {
    return t8_eclass_num_vertices[eclass];
  }
  return t8_eclass_to_dimension[eclass] == 3 ? t8_eclass_num_edges[eclass] : 0;
}

/** The corner permutation of a face connection, seen from the tree with \a eclass. */
static void
t8_cmesh_face_corner_permutation (const t8_eclass_t eclass, const int face, const t8_eclass_t neigh_eclass,
                                  const int neigh_face, const int orientation, int *perm)
{
  const int sign = t8_eclass_face_orientation[eclass][face] == t8_eclass_face_orientation[neigh_eclass][neigh_face];
  /* The orientation is defined relative to the smaller face, see t8_forest_element_face_neighbor. */
  const int eclass_compare = t8_eclass_compare (eclass, neigh_eclass);
  const int is_smaller = eclass_compare == 0 ? face <= neigh_face : eclass_compare == -1;
  t8_eclass_face_orientation_to_corner_permutation ((t8_eclass_t) t8_eclass_face_types[eclass][face], orientation,
                                                    sign, is_smaller, perm);
}

/** Send an array of global indices to each process and receive one from each process.
 * \param [in] send   The data for each process, of length mpisize.
 * \param [in] comm   The communicator.
 * \return            The data received from each process.
 */
static std::vector<std::vector<t8_gloidx_t>>
t8_cmesh_exchange (const std::vector<std::vector<t8_gloidx_t>> &send, sc_MPI_Comm comm)
{
  const int mpisize = send.size ();
  std::vector<int> send_counts (mpisize), recv_counts (mpisize);
  for (int rank = 0; rank < mpisize; rank++) {
    send_counts[rank] = send[rank].size ();
  }
  int mpiret = sc_MPI_Alltoall (send_counts.data (), 1, sc_MPI_INT, recv_counts.data (), 1, sc_MPI_INT, comm);
  SC_CHECK_MPI (mpiret);
  std::vector<std::vector<t8_gloidx_t>> recv (mpisize);
  /* Reserved, since MPI keeps pointers to the requests */
  std::vector<sc_MPI_Request> requests;
  requests.reserve (2 * mpisize);
  for (int rank = 0; rank < mpisize; rank++) {
    if (recv_counts[rank] > 0) {
      recv[rank].resize (recv_counts[rank]);
      requests.emplace_back ();
      mpiret = sc_MPI_Irecv (recv[rank].data (), recv_counts[rank], T8_MPI_GLOIDX, rank,
                             T8_MPI_CMESH_BOUNDARY_CONNECTIVITY, comm, &requests.back ());
      SC_CHECK_MPI (mpiret);
    }
  }
  for (int rank = 0; rank < mpisize; rank++) {
    if (send_counts[rank] > 0) {
      requests.emplace_back ();
      mpiret = sc_MPI_Isend ((void *) send[rank].data (), send_counts[rank], T8_MPI_GLOIDX, rank,
                             T8_MPI_CMESH_BOUNDARY_CONNECTIVITY, comm, &requests.back ());
      SC_CHECK_MPI (mpiret);
    }
  }
  mpiret = sc_MPI_Waitall (requests.size (), requests.data (), sc_MPI_STATUSES_IGNORE);
  SC_CHECK_MPI (mpiret);
  return recv;
}

void
t8_cmesh_boundary_connectivity::build_classes (const int bdy_dim, sc_MPI_Comm comm)
{
  classes &result = bdy_dim == 0 ? vertices : edges;
  const bool partitioned = t8_cmesh_is_partitioned (cmesh);
  const int num_procs = partitioned ? cmesh->mpisize : 1;
  const int my_rank = partitioned ? cmesh->mpirank : 0;
  const t8_gloidx_t *tree_offsets = partitioned ? t8_shmem_array_get_gloidx_array (cmesh->tree_offsets) : NULL;
  /* Exchange data between the processes, or keep it on this process if the cmesh is replicated */
  auto exchange = [&] (const std::vector<std::vector<t8_gloidx_t>> &send) {
    return partitioned ? t8_cmesh_exchange (send, comm) : send;
  };
  /* The process that collects the data of a tree, the first one that has it as local tree */
  auto tree_owner = [&] (const t8_gloidx_t gtreeid) {
    int some_owner = -1;
    return partitioned ? t8_offset_first_owner_of_tree (num_procs, gtreeid, tree_offsets, &some_owner) : 0;
  };

  /* The vertices or edges of all local trees and ghosts */
  const t8_locidx_t num_local_trees = t8_cmesh_get_num_local_trees (cmesh);
  const t8_locidx_t num_trees = num_local_trees + (partitioned ? t8_cmesh_get_num_ghosts (cmesh) : 0);
  std::vector<t8_locidx_t> offset (num_trees + 1, 0);
  for (t8_locidx_t itree = 0; itree < num_trees; itree++) {
    offset[itree + 1]
      = offset[itree] + t8_cmesh_num_tree_boundaries (t8_cmesh_local_or_ghost_class (cmesh, itree), bdy_dim);
  }
  t8_parity_union_find sets (offset[num_trees]);

  /* Identify the vertices or edges of all known face connections */
  for (t8_locidx_t itree = 0; itree < num_trees; itree++) {
    const t8_eclass_t eclass = t8_cmesh_local_or_ghost_class (cmesh, itree);
    if (offset[itree + 1] == offset[itree]) {
      continue;
    }
    for (int iface = 0; iface < t8_eclass_num_faces[eclass]; iface++) {
      int neigh_face, orientation;
      const t8_locidx_t neigh_tree = t8_cmesh_get_face_neighbor (cmesh, itree, iface, &neigh_face, &orientation);
      if (neigh_tree < 0 || neigh_tree >= num_trees || (neigh_tree == itree && neigh_face == iface)) {
        continue;
      }
      const t8_eclass_t neigh_eclass = t8_cmesh_local_or_ghost_class (cmesh, neigh_tree);
      int perm[T8_ECLASS_MAX_CORNERS_2D];
      t8_cmesh_face_corner_permutation (eclass, iface, neigh_eclass, neigh_face, orientation, perm);
      const int num_face_corners = t8_eclass_num_vertices[t8_eclass_face_types[eclass][iface]];
      /* The vertex of the neighbor tree at each tree vertex of this face, or -1. */
      int neigh_vertex[T8_ECLASS_MAX_CORNERS];
      std::fill (neigh_vertex, neigh_vertex + T8_ECLASS_MAX_CORNERS, -1);
      for (int icorner = 0; icorner < num_face_corners; icorner++) {
        const int vertex = t8_face_vertex_to_tree_vertex[eclass][iface][icorner];
        neigh_vertex[vertex] = t8_face_vertex_to_tree_vertex[neigh_eclass][neigh_face][perm[icorner]];
        if (bdy_dim == 0) {
          sets.join (offset[itree] + vertex, offset[neigh_tree] + neigh_vertex[vertex], 0);
        }
      }
      if (bdy_dim == 0) {
        continue;
      }
      /* The face edges are mapped onto the face edges of the neighbor face with the same vertices. */
      for (int iedge = 0; iedge < num_face_corners; iedge++) {
        const int edge = t8_face_edge_to_tree_edge[eclass][iface][iedge];
        const int image_first = neigh_vertex[t8_edge_vertex_to_tree_vertex[eclass][edge][0]];
        const int image_second = neigh_vertex[t8_edge_vertex_to_tree_vertex[eclass][edge][1]];
        [[maybe_unused]] bool found = false;
        for (int ineigh_edge = 0; ineigh_edge < num_face_corners; ineigh_edge++) {
          const int neigh_edge = t8_face_edge_to_tree_edge[neigh_eclass][neigh_face][ineigh_edge];
          const int *neigh_edge_vertices = t8_edge_vertex_to_tree_vertex[neigh_eclass][neigh_edge];
          if ((neigh_edge_vertices[0] == image_first && neigh_edge_vertices[1] == image_second)
              || (neigh_edge_vertices[0] == image_second && neigh_edge_vertices[1] == image_first)) {
            sets.join (offset[itree] + edge, offset[neigh_tree] + neigh_edge, neigh_edge_vertices[0] != image_first);
            found = true;
          }
        }
        T8_ASSERT (found);
      }
    }
  }

  /* In a replicated cmesh, identify the vertices with the same global vertex id, which also covers trees that only
   * share vertices. */
  if (bdy_dim == 0 && !partitioned && t8_cmesh_uses_vertex_connectivity (cmesh)) {
    std::unordered_map<t8_gloidx_t, size_t> first_vertex;
    for (t8_locidx_t itree = 0; itree < num_local_trees; itree++) {
      for (t8_locidx_t ivertex = offset[itree]; ivertex < offset[itree + 1]; ivertex++) {
        const t8_gloidx_t global_vertex = t8_cmesh_get_global_vertex_of_tree (cmesh, itree, ivertex - offset[itree]);
        const auto [entry, is_new] = first_vertex.try_emplace (global_vertex, ivertex);
        if (!is_new) {
          sets.join (entry->second, ivertex, 0);
        }
      }
    }
  }

  /* The label of each set is the smallest key of a vertex or edge in it, with the parity of the representative
   * relative to the vertex or edge of the label. In a partitioned cmesh, sets of different processes with a common
   * tree get the same label below. */
  const t8_gloidx_t no_label = std::numeric_limits<t8_gloidx_t>::max ();
  std::vector<t8_gloidx_t> label (offset[num_trees], no_label);
  std::vector<int> label_parity (offset[num_trees], 0);
  std::vector<t8_gloidx_t> gtreeids (num_trees);
  for (t8_locidx_t itree = 0; itree < num_trees; itree++) {
    gtreeids[itree] = t8_cmesh_get_global_id (cmesh, itree);
  }
  /* Lower the label of the set of a vertex or edge, given its label and its parity relative to it. */
  auto lower_label = [&] (const size_t index, const t8_gloidx_t new_label, const int parity) {
    int index_parity;
    const size_t root = sets.find (index, &index_parity);
    if (new_label < label[root]) {
      label[root] = new_label;
      label_parity[root] = parity ^ index_parity;
      return true;
    }
    return false;
  };
  /* The label of a vertex or edge and its parity relative to the one of the label */
  auto get_label = [&] (const size_t index, t8_gloidx_t *index_label, int *parity) {
    int index_parity;
    const size_t root = sets.find (index, &index_parity);
    *index_label = label[root];
    *parity = index_parity ^ label_parity[root];
  };
  for (t8_locidx_t itree = 0; itree < num_trees; itree++) {
    for (t8_locidx_t index = offset[itree]; index < offset[itree + 1]; index++) {
      lower_label (index, T8_BOUNDARY_KEY_STRIDE * gtreeids[itree] + index - offset[itree], 0);
    }
  }

  if (partitioned) {
    /* Each process sends the labels of the trees of which it is not the first owner to the first owner, which
     * lowers its labels and replies with them, until no label changes. */
    int changed_global = 1;
    while (changed_global) {
      int changed = 0;
      std::vector<std::vector<t8_gloidx_t>> requests (num_procs);
      std::vector<std::vector<size_t>> request_indices (num_procs);
      for (t8_locidx_t itree = 0; itree < num_trees; itree++) {
        const int owner = tree_owner (gtreeids[itree]);
        if (owner == my_rank) {
          continue;
        }
        for (t8_locidx_t index = offset[itree]; index < offset[itree + 1]; index++) {
          t8_gloidx_t index_label;
          int parity;
          get_label (index, &index_label, &parity);
          requests[owner].insert (requests[owner].end (),
                                  { gtreeids[itree], index - offset[itree], index_label, parity });
          request_indices[owner].push_back (index);
        }
      }
      const auto received = exchange (requests);
      /* First lower all labels, then reply with the lowered labels */
      std::vector<std::vector<size_t>> received_indices (num_procs);
      for (int rank = 0; rank < num_procs; rank++) {
        for (size_t ientry = 0; ientry < received[rank].size (); ientry += 4) {
          const t8_locidx_t ltreeid = t8_cmesh_get_local_id (cmesh, received[rank][ientry]);
          T8_ASSERT (0 <= ltreeid && ltreeid < num_local_trees);
          const size_t index = offset[ltreeid] + received[rank][ientry + 1];
          changed |= lower_label (index, received[rank][ientry + 2], received[rank][ientry + 3]);
          received_indices[rank].push_back (index);
        }
      }
      std::vector<std::vector<t8_gloidx_t>> replies (num_procs);
      for (int rank = 0; rank < num_procs; rank++) {
        for (const size_t index : received_indices[rank]) {
          t8_gloidx_t index_label;
          int parity;
          get_label (index, &index_label, &parity);
          replies[rank].insert (replies[rank].end (), { index_label, parity });
        }
      }
      const auto replied = exchange (replies);
      for (int rank = 0; rank < num_procs; rank++) {
        for (size_t ientry = 0; ientry < request_indices[rank].size (); ientry++) {
          changed |= lower_label (request_indices[rank][ientry], replied[rank][2 * ientry],
                                  replied[rank][2 * ientry + 1]);
        }
      }
      const int mpiret = sc_MPI_Allreduce (&changed, &changed_global, 1, sc_MPI_INT, sc_MPI_LOR, comm);
      SC_CHECK_MPI (mpiret);
    }
  }

  /* Collect the members of each class on the owner of the tree of the label. The first owner of each local tree
   * sends its vertices or edges as members, and each process requests the classes of its local trees. */
  std::vector<std::vector<t8_gloidx_t>> collect (num_procs);
  for (t8_locidx_t itree = 0; itree < num_local_trees; itree++) {
    const bool is_first_owner = tree_owner (gtreeids[itree]) == my_rank;
    const t8_eclass_t eclass = t8_cmesh_get_tree_class (cmesh, itree);
    for (t8_locidx_t index = offset[itree]; index < offset[itree + 1]; index++) {
      t8_gloidx_t index_label;
      int parity;
      get_label (index, &index_label, &parity);
      std::vector<t8_gloidx_t> &message = collect[tree_owner (index_label / T8_BOUNDARY_KEY_STRIDE)];
      if (is_first_owner) {
        message.insert (message.end (), { 0, index_label, gtreeids[itree], eclass, index - offset[itree], parity });
      }
      message.insert (message.end (), { 1, index_label });
    }
  }
  const auto collected = exchange (collect);
  std::map<t8_gloidx_t, std::vector<member>> class_members;
  std::vector<std::set<t8_gloidx_t>> requested (num_procs);
  for (int rank = 0; rank < num_procs; rank++) {
    for (size_t ientry = 0; ientry < collected[rank].size ();) {
      if (collected[rank][ientry] == 0) {
        const auto *entry = collected[rank].data () + ientry;
        class_members[entry[1]].push_back (
          { entry[2], (t8_eclass_t) entry[3], (int) entry[4], (int) entry[5] });
        ientry += 6;
      }
      else {
        requested[rank].insert (collected[rank][ientry + 1]);
        ientry += 2;
      }
    }
  }
  std::vector<std::vector<t8_gloidx_t>> send_classes (num_procs);
  for (int rank = 0; rank < num_procs; rank++) {
    for (const t8_gloidx_t class_label : requested[rank]) {
      const std::vector<member> &members = class_members[class_label];
      send_classes[rank].insert (send_classes[rank].end (), { class_label, (t8_gloidx_t) members.size () });
      for (const member &m : members) {
        send_classes[rank].insert (send_classes[rank].end (), { m.gtreeid, m.eclass, m.bdy_id, m.parity });
      }
    }
  }
  const auto received_classes = exchange (send_classes);
  std::unordered_map<t8_gloidx_t, int> class_index;
  for (int rank = 0; rank < num_procs; rank++) {
    for (size_t ientry = 0; ientry < received_classes[rank].size ();) {
      const t8_gloidx_t class_label = received_classes[rank][ientry];
      const t8_gloidx_t num_members = received_classes[rank][ientry + 1];
      ientry += 2;
      class_index[class_label] = result.members.size ();
      result.members.emplace_back ();
      for (t8_gloidx_t imember = 0; imember < num_members; imember++, ientry += 4) {
        const auto *entry = received_classes[rank].data () + ientry;
        result.members.back ().push_back ({ entry[0], (t8_eclass_t) entry[1], (int) entry[2], (int) entry[3] });
      }
    }
  }

  /* The class and parity of each vertex or edge of the local trees */
  result.offset.assign (offset.begin (), offset.begin () + num_local_trees + 1);
  result.class_of.resize (offset[num_local_trees]);
  result.parity.resize (offset[num_local_trees]);
  for (t8_locidx_t index = 0; index < offset[num_local_trees]; index++) {
    t8_gloidx_t index_label;
    get_label (index, &index_label, &result.parity[index]);
    result.class_of[index] = class_index.at (index_label);
  }
}

t8_cmesh_boundary_connectivity::t8_cmesh_boundary_connectivity (const t8_cmesh_t cmesh, sc_MPI_Comm comm)
  : cmesh (cmesh)
{
  T8_ASSERT (t8_cmesh_is_committed (cmesh));
  build_classes (0, comm);
  build_classes (1, comm);
}

std::vector<t8_cmesh_boundary_neighbor>
t8_cmesh_boundary_connectivity::get_neighbors (const t8_locidx_t ltreeid, const int bdy_dim, const int bdy_id) const
{
  std::vector<t8_cmesh_boundary_neighbor> neighbors;
  const t8_eclass_t eclass = t8_cmesh_get_tree_class (cmesh, ltreeid);
  const t8_gloidx_t gtreeid = t8_cmesh_get_global_id (cmesh, ltreeid);
  const int dim = t8_eclass_to_dimension[eclass];
  T8_ASSERT (0 <= bdy_dim && bdy_dim < dim);
  if (bdy_dim == dim - 1) {
    T8_ASSERT (0 <= bdy_id && bdy_id < t8_eclass_num_faces[eclass]);
    int neigh_face, orientation;
    const t8_locidx_t neigh_tree = t8_cmesh_get_face_neighbor (cmesh, ltreeid, bdy_id, &neigh_face, &orientation);
    if (neigh_tree >= 0 && !(neigh_tree == ltreeid && neigh_face == bdy_id)) {
      const t8_eclass_t neigh_eclass = t8_cmesh_local_or_ghost_class (cmesh, neigh_tree);
      t8_cmesh_boundary_neighbor neighbor { t8_cmesh_get_global_id (cmesh, neigh_tree), neigh_eclass, neigh_face, {} };
      t8_cmesh_face_corner_permutation (eclass, bdy_id, neigh_eclass, neigh_face, orientation, neighbor.perm);
      neighbors.push_back (neighbor);
    }
    return neighbors;
  }
  T8_ASSERT (bdy_dim == 0 || (bdy_dim == 1 && dim == 3));
  const classes &bdy_classes = bdy_dim == 0 ? vertices : edges;
  const t8_locidx_t index = bdy_classes.offset[ltreeid] + bdy_id;
  T8_ASSERT (index < bdy_classes.offset[ltreeid + 1]);
  for (const member &m : bdy_classes.members[bdy_classes.class_of[index]]) {
    if (m.gtreeid != gtreeid || m.bdy_id != bdy_id) {
      const int flip = m.parity != bdy_classes.parity[index];
      neighbors.push_back ({ m.gtreeid, m.eclass, m.bdy_id, { flip, !flip } });
    }
  }
  return neighbors;
}

void
t8_cmesh_build_boundary_connectivity (const t8_cmesh_t cmesh, sc_MPI_Comm comm)
{
  T8_ASSERT (t8_cmesh_is_committed (cmesh));
  if (cmesh->boundary_connectivity == NULL) {
    cmesh->boundary_connectivity = new t8_cmesh_boundary_connectivity (cmesh, comm);
  }
}

std::vector<t8_cmesh_boundary_neighbor>
t8_cmesh_get_boundary_neighbors (const t8_cmesh_t cmesh, const t8_locidx_t ltreeid, const int bdy_dim,
                                 const int bdy_id)
{
  T8_ASSERT (t8_cmesh_is_committed (cmesh));
  T8_ASSERT (t8_cmesh_treeid_is_local_tree (cmesh, ltreeid));
  if (cmesh->boundary_connectivity == NULL) {
    SC_CHECK_ABORT (!t8_cmesh_is_partitioned (cmesh),
                    "The boundary connectivity of a partitioned cmesh must be built with "
                    "t8_cmesh_build_boundary_connectivity.\n");
    cmesh->boundary_connectivity = new t8_cmesh_boundary_connectivity (cmesh, sc_MPI_COMM_SELF);
  }
  return cmesh->boundary_connectivity->get_neighbors (ltreeid, bdy_dim, bdy_id);
}
