/*
  This file is part of t8code.
  t8code is a C library to manage a collection (a forest) of multiple
  connected adaptive space-trees of general element classes in parallel.

  Copyright (C) 2026 the developers

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

/** \file t8_dgm.hxx
 * Delaunay Graph Mapping (DGM) used to interpolate the boundary displacements of a mesh deformation to the inner nodes of each Delaunay element.
 */
#pragma once

#include <array>
#include <vector>
#include <cmath>
#include <limits>
#include <memory>
#include <unordered_map>
#include <t8_vector/t8_vec.hxx>
#include <t8_cmesh/t8_cmesh_mesh_deformation/t8_cmesh_mesh_deformation_method.hxx>

/**
 * Struct for mesh deformation using Delaunay Graph Mapping (DGM).
 * The boundary nodes are triangulated with a Delaunay triangulation (the Delaunay graph). Each inner node is
 * located in a Delaunay element and its displacement is interpolated from the element's boundary nodes.
 */
struct t8_dgm: public t8_mesh_deformation_method
{
 public:
  /** Constructor. */
  t8_dgm () {};

  /** Destructor. */
  ~t8_dgm () {};

  /**
   * Transfers the boundary node data to the internal data structure of the DGM.
   * \param [in] boundary_node_data Map of global vertex IDs to boundary nodes.
   */
  void
  set_boundary_nodes (const std::unordered_map<t8_gloidx_t, t8_mesh_deformation_node>& boundary_node_data) override
  {
    boundary_nodes.clear ();
    boundary_nodes.reserve (boundary_node_data.size ());
    for (const auto& [global_vertex_id, boundary_node] : boundary_node_data) {
      boundary_nodes.push_back (boundary_node);
      boundary_nodes.back ().global_id = global_vertex_id;
    }
  }

  /**
   * Builds the Delaunay graph of the boundary nodes.
   * This step is mandatory to later be able to interpolate the inner nodes.
   */
  void
  solve () override;

  /**
   * Interpolates the inner node.
   * \param [in, out] inner_node   The inner node which will be interpolated.
   */
  void
  interpolate (t8_mesh_deformation_node& inner_node) const override
  {
  }

  /**
   * Returns the Delaunay triangles of the boundary nodes.
   * Each triangle holds three indices into the boundary node list.
   */
  const std::vector<std::array<size_t, 3>>&
  get_delaunay_elements () const
  {
    return delaunay_elements;
  }

 private:
  /**
   * Computes the 2D Delaunay triangulation of the boundary nodes with the Bowyer-Watson algorithm
   * and stores the resulting triangles in \a delaunay_elements.
   */
  void
  create_delaunay_mesh (const std::vector<t8_mesh_deformation_node>& boundary_node_data);

  /** List of all registered support nodes. */
  std::vector<t8_mesh_deformation_node> boundary_nodes;
  /** The Delaunay triangles. Each entry holds three indices into \a boundary_nodes. */
  std::vector<std::array<size_t, 3>> delaunay_elements;
};
