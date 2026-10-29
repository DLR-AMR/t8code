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
 * Vertices and edges are identified across all face connections of the cmesh, edges with their direction, so this
 * also works for periodic meshes, in which both vertices of different edges may coincide. In a replicated cmesh,
 * vertices are additionally identified by the global vertex ids if these were set
 * (\ref t8_cmesh_set_global_vertices_of_tree).
 *
 * Each vertex and edge class is labeled by the smallest key 12 * global tree id + vertex or edge in it. In a
 * partitioned cmesh, each process labels the classes of its local trees and ghosts, and the processes exchange the
 * labels of trees that are present on several processes with the first owner of the tree until the labels agree.
 * The members of each class are then collected on the first owner of the tree of the label, which sends them to all
 * processes with a local tree in the class.
 */

#ifndef T8_CMESH_BOUNDARY_CONNECTIVITY_HXX
#define T8_CMESH_BOUNDARY_CONNECTIVITY_HXX

#include <t8_cmesh/t8_cmesh.hxx>
#include <vector>

/** The vertex, edge and face connectivity of the local trees of a cmesh. */
struct t8_cmesh_boundary_connectivity
{
 public:
  /** Build the connectivity of a committed cmesh.
   * \param [in] cmesh  A committed cmesh.
   * \param [in] comm   The communicator of \a cmesh, only used if \a cmesh is partitioned.
   * \note Collective if \a cmesh is partitioned.
   */
  t8_cmesh_boundary_connectivity (const t8_cmesh_t cmesh, sc_MPI_Comm comm);

  /** The trees sharing a boundary with a local tree, see \ref t8_cmesh_get_boundary_neighbors. */
  std::vector<t8_cmesh_boundary_neighbor>
  get_neighbors (const t8_locidx_t ltreeid, const int bdy_dim, const int bdy_id) const;

 private:
  /** A member of a vertex or edge class. */
  struct member
  {
    t8_gloidx_t gtreeid; /**< The global id of the tree. */
    t8_eclass_t eclass;  /**< The eclass of the tree. */
    int bdy_id;          /**< The vertex or edge of the tree. */
    int parity;          /**< For edges, 1 if the edge runs opposite to the first edge of the class. */
  };

  /** The classes of the vertices (bdy_dim 0) or edges (bdy_dim 1) of the local trees. */
  struct classes
  {
    std::vector<t8_locidx_t> offset;          /**< The index of the first vertex or edge of each local tree. */
    std::vector<int> class_of;                /**< The class of each vertex or edge of the local trees. */
    std::vector<int> parity;                  /**< The parity of each vertex or edge of the local trees. */
    std::vector<std::vector<member>> members; /**< The members of each class. */
  };

  /** Compute the classes of vertices or edges.
   * \param [in] bdy_dim  0 for vertices, 1 for edges.
   * \param [in] comm     The communicator, used if the cmesh is partitioned.
   */
  void
  build_classes (const int bdy_dim, sc_MPI_Comm comm);

  t8_cmesh_t cmesh;      /**< The cmesh of this connectivity. */
  classes vertices;      /**< The classes of the vertices. */
  classes edges;         /**< The classes of the edges. */
};

#endif /* !T8_CMESH_BOUNDARY_CONNECTIVITY_HXX */
