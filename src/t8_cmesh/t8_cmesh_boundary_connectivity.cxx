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
#include <t8_cmesh/t8_cmesh_vertex_connectivity/t8_cmesh_vertex_connectivity.hxx>
#include <algorithm>
#include <unordered_map>

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

/** The corner permutation of a face connection, seen from \a ltreeid. */
static void
t8_cmesh_face_corner_permutation (const t8_cmesh_t cmesh, const t8_locidx_t ltreeid, const int face,
                                  const t8_locidx_t neigh_ltreeid, const int neigh_face, const int orientation,
                                  int *perm)
{
  const t8_eclass_t eclass = t8_cmesh_get_tree_class (cmesh, ltreeid);
  const t8_eclass_t neigh_eclass = t8_cmesh_get_tree_class (cmesh, neigh_ltreeid);
  const int sign = t8_eclass_face_orientation[eclass][face] == t8_eclass_face_orientation[neigh_eclass][neigh_face];
  /* The orientation is defined relative to the smaller face, see t8_forest_element_face_neighbor. */
  const int eclass_compare = t8_eclass_compare (eclass, neigh_eclass);
  const int is_smaller = eclass_compare == 0 ? face <= neigh_face : eclass_compare == -1;
  t8_eclass_face_orientation_to_corner_permutation ((t8_eclass_t) t8_eclass_face_types[eclass][face], orientation,
                                                    sign, is_smaller, perm);
}

t8_cmesh_boundary_connectivity::t8_cmesh_boundary_connectivity (const t8_cmesh_t cmesh): cmesh (cmesh)
{
  T8_ASSERT (t8_cmesh_is_committed (cmesh));
  SC_CHECK_ABORT (!t8_cmesh_is_partitioned (cmesh), "Boundary connectivity requires a replicated cmesh.\n");
  const t8_locidx_t num_trees = t8_cmesh_get_num_local_trees (cmesh);
  vertex_offset.resize (num_trees + 1, 0);
  edge_offset.resize (num_trees + 1, 0);
  for (t8_locidx_t itree = 0; itree < num_trees; itree++) {
    const t8_eclass_t eclass = t8_cmesh_get_tree_class (cmesh, itree);
    vertex_offset[itree + 1] = vertex_offset[itree] + t8_eclass_num_vertices[eclass];
    edge_offset[itree + 1] = edge_offset[itree] + (t8_eclass_to_dimension[eclass] == 3 ? t8_eclass_num_edges[eclass] : 0);
  }
  t8_parity_union_find vertices (vertex_offset[num_trees]);
  t8_parity_union_find edges (edge_offset[num_trees]);

  /* Identify the vertices and edges of all connected faces. */
  for (t8_locidx_t itree = 0; itree < num_trees; itree++) {
    const t8_eclass_t eclass = t8_cmesh_get_tree_class (cmesh, itree);
    for (int iface = 0; iface < t8_eclass_num_faces[eclass]; iface++) {
      int neigh_face, orientation;
      const t8_locidx_t neigh_tree = t8_cmesh_get_face_neighbor (cmesh, itree, iface, &neigh_face, &orientation);
      if (neigh_tree < 0 || (neigh_tree == itree && neigh_face == iface)) {
        continue;
      }
      const t8_eclass_t neigh_eclass = t8_cmesh_get_tree_class (cmesh, neigh_tree);
      int perm[T8_ECLASS_MAX_CORNERS_2D];
      t8_cmesh_face_corner_permutation (cmesh, itree, iface, neigh_tree, neigh_face, orientation, perm);
      const int num_face_corners = t8_eclass_num_vertices[t8_eclass_face_types[eclass][iface]];
      /* The vertex of the neighbor tree at each tree vertex of this face, or -1. */
      int neigh_vertex[T8_ECLASS_MAX_CORNERS];
      std::fill (neigh_vertex, neigh_vertex + T8_ECLASS_MAX_CORNERS, -1);
      for (int icorner = 0; icorner < num_face_corners; icorner++) {
        const int vertex = t8_face_vertex_to_tree_vertex[eclass][iface][icorner];
        neigh_vertex[vertex] = t8_face_vertex_to_tree_vertex[neigh_eclass][neigh_face][perm[icorner]];
        vertices.join (vertex_offset[itree] + vertex, vertex_offset[neigh_tree] + neigh_vertex[vertex], 0);
      }
      if (t8_eclass_to_dimension[eclass] < 3) {
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
            edges.join (edge_offset[itree] + edge, edge_offset[neigh_tree] + neigh_edge,
                        neigh_edge_vertices[0] != image_first);
            found = true;
          }
        }
        T8_ASSERT (found);
      }
    }
  }

  /* Identify the vertices with the same global vertex id, which also covers trees that only share vertices. */
  if (t8_cmesh_uses_vertex_connectivity (cmesh)) {
    std::unordered_map<t8_gloidx_t, size_t> first_vertex;
    for (t8_locidx_t itree = 0; itree < num_trees; itree++) {
      const t8_eclass_t eclass = t8_cmesh_get_tree_class (cmesh, itree);
      for (int ivertex = 0; ivertex < t8_eclass_num_vertices[eclass]; ivertex++) {
        const t8_gloidx_t global_vertex = t8_cmesh_get_global_vertex_of_tree (cmesh, itree, ivertex);
        const auto [entry, is_new] = first_vertex.try_emplace (global_vertex, vertex_offset[itree] + ivertex);
        if (!is_new) {
          vertices.join (entry->second, vertex_offset[itree] + ivertex, 0);
        }
      }
    }
  }

  /* Number the classes and collect their members. */
  vertex_class.resize (vertex_offset[num_trees]);
  edge_class.resize (edge_offset[num_trees]);
  edge_parity.resize (edge_offset[num_trees]);
  std::unordered_map<size_t, int> vertex_root_to_class, edge_root_to_class;
  for (t8_locidx_t itree = 0; itree < num_trees; itree++) {
    for (t8_locidx_t ivertex = vertex_offset[itree]; ivertex < vertex_offset[itree + 1]; ivertex++) {
      int parity;
      const size_t root = vertices.find (ivertex, &parity);
      const auto [entry, is_new] = vertex_root_to_class.try_emplace (root, vertex_class_members.size ());
      if (is_new) {
        vertex_class_members.emplace_back ();
      }
      vertex_class[ivertex] = entry->second;
      vertex_class_members[entry->second].emplace_back (itree, ivertex - vertex_offset[itree]);
    }
    for (t8_locidx_t iedge = edge_offset[itree]; iedge < edge_offset[itree + 1]; iedge++) {
      const size_t root = edges.find (iedge, &edge_parity[iedge]);
      const auto [entry, is_new] = edge_root_to_class.try_emplace (root, edge_class_members.size ());
      if (is_new) {
        edge_class_members.emplace_back ();
      }
      edge_class[iedge] = entry->second;
      edge_class_members[entry->second].emplace_back (itree, iedge - edge_offset[itree]);
    }
  }
}

std::vector<t8_cmesh_boundary_neighbor>
t8_cmesh_boundary_connectivity::get_neighbors (const t8_locidx_t ltreeid, const int bdy_dim, const int bdy_id) const
{
  std::vector<t8_cmesh_boundary_neighbor> neighbors;
  const t8_eclass_t eclass = t8_cmesh_get_tree_class (cmesh, ltreeid);
  const int dim = t8_eclass_to_dimension[eclass];
  T8_ASSERT (0 <= bdy_dim && bdy_dim < dim);
  if (bdy_dim == dim - 1) {
    T8_ASSERT (0 <= bdy_id && bdy_id < t8_eclass_num_faces[eclass]);
    int neigh_face, orientation;
    const t8_locidx_t neigh_tree = t8_cmesh_get_face_neighbor (cmesh, ltreeid, bdy_id, &neigh_face, &orientation);
    if (neigh_tree >= 0 && !(neigh_tree == ltreeid && neigh_face == bdy_id)) {
      t8_cmesh_boundary_neighbor neighbor { neigh_tree, neigh_face, {} };
      t8_cmesh_face_corner_permutation (cmesh, ltreeid, bdy_id, neigh_tree, neigh_face, orientation, neighbor.perm);
      neighbors.push_back (neighbor);
    }
  }
  else if (bdy_dim == 0) {
    T8_ASSERT (0 <= bdy_id && bdy_id < t8_eclass_num_vertices[eclass]);
    for (const auto &[neigh_tree, neigh_vertex] : vertex_class_members[vertex_class[vertex_offset[ltreeid] + bdy_id]]) {
      if (neigh_tree != ltreeid || neigh_vertex != bdy_id) {
        neighbors.push_back ({ neigh_tree, neigh_vertex, { 0 } });
      }
    }
  }
  else {
    T8_ASSERT (bdy_dim == 1 && dim == 3);
    T8_ASSERT (0 <= bdy_id && bdy_id < t8_eclass_num_edges[eclass]);
    const t8_locidx_t edge = edge_offset[ltreeid] + bdy_id;
    for (const auto &[neigh_tree, neigh_edge] : edge_class_members[edge_class[edge]]) {
      if (neigh_tree != ltreeid || neigh_edge != bdy_id) {
        const int flip = edge_parity[edge] != edge_parity[edge_offset[neigh_tree] + neigh_edge];
        neighbors.push_back ({ neigh_tree, neigh_edge, { flip, !flip } });
      }
    }
  }
  return neighbors;
}

std::vector<t8_cmesh_boundary_neighbor>
t8_cmesh_get_boundary_neighbors (const t8_cmesh_t cmesh, const t8_locidx_t ltreeid, const int bdy_dim,
                                 const int bdy_id)
{
  T8_ASSERT (t8_cmesh_is_committed (cmesh));
  T8_ASSERT (t8_cmesh_treeid_is_local_tree (cmesh, ltreeid));
  if (cmesh->boundary_connectivity == NULL) {
    cmesh->boundary_connectivity = new t8_cmesh_boundary_connectivity (cmesh);
  }
  return cmesh->boundary_connectivity->get_neighbors (ltreeid, bdy_dim, bdy_id);
}
