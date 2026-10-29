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

/** \file t8_cmesh_boundary_connectivity.hxx
 * The trees sharing a vertex, edge or face with a tree, together with the permutation of the boundary corners.
 * Vertices and edges are identified across all face connections of the cmesh, and additionally by the global
 * vertex ids if these were set (\ref t8_cmesh_set_global_vertices_of_tree). Edges are identified with an
 * orientation, so this also works for periodic meshes, in which both vertices of different edges may coincide.
 */

#ifndef T8_CMESH_BOUNDARY_CONNECTIVITY_HXX
#define T8_CMESH_BOUNDARY_CONNECTIVITY_HXX

#include <t8_cmesh/t8_cmesh.hxx>
#include <utility>
#include <vector>

/** The vertex, edge and face connectivity of the trees of a replicated cmesh. */
struct t8_cmesh_boundary_connectivity
{
 public:
  /** Build the connectivity of a committed, replicated cmesh.
   * \param [in] cmesh  A committed cmesh that is not partitioned.
   */
  explicit t8_cmesh_boundary_connectivity (const t8_cmesh_t cmesh);

  /** The trees sharing a boundary with a tree, see \ref t8_cmesh_get_boundary_neighbors. */
  std::vector<t8_cmesh_boundary_neighbor>
  get_neighbors (const t8_locidx_t ltreeid, const int bdy_dim, const int bdy_id) const;

 private:
  /** A tree and one of its vertices or edges. */
  using tree_boundary = std::pair<t8_locidx_t, int>;

  t8_cmesh_t cmesh; /**< The cmesh of this connectivity. */
  std::vector<t8_locidx_t> vertex_offset; /**< The index of the first vertex of each tree. */
  std::vector<t8_locidx_t> edge_offset;   /**< The index of the first edge of each tree. */
  std::vector<int> vertex_class;          /**< The class of identified vertices of each tree vertex. */
  std::vector<int> edge_class;            /**< The class of identified edges of each tree edge. */
  /** 1 if the edge runs opposite to the first edge of its class. */
  std::vector<int> edge_parity;
  std::vector<std::vector<tree_boundary>> vertex_class_members; /**< The tree vertices of each class. */
  std::vector<std::vector<tree_boundary>> edge_class_members;   /**< The tree edges of each class. */
};

#endif /* !T8_CMESH_BOUNDARY_CONNECTIVITY_HXX */
