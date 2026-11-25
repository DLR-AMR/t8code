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

/** \file t8_dg_rbf.cxx
 *  This file implements the Delaunay graph based RBF method (DG-RBF) for the mesh deformation.
 */
#include <t8.h>
#include <t8_cmesh/t8_cmesh_mesh_deformation/t8_dg_rbf.hxx>
#include <t8_cmesh/t8_cmesh_mesh_deformation/t8_dgm.hxx>
#include <t8_vector/t8_vec.hxx>
#include <gmsh.h>
#include <cmath>

/** Rotates a position about a center, counter-clockwise about the x-, then the y- and then the z-axis. */
static t8_3D_vec
t8_dg_rbf_rotate (const t8_3D_vec& position, const t8_3D_vec& center, const t8_3D_vec& angles)
{
  t8_3D_vec relative = { position[0] - center[0], position[1] - center[1], position[2] - center[2] };
  for (int axis = 0; axis < 3; ++axis) {
    const int first = (axis + 1) % 3;
    const int second = (axis + 2) % 3;
    const double rotated_first = std::cos (angles[axis]) * relative[first] - std::sin (angles[axis]) * relative[second];
    relative[second] = std::sin (angles[axis]) * relative[first] + std::cos (angles[axis]) * relative[second];
    relative[first] = rotated_first;
  }
  return { relative[0] + center[0], relative[1] + center[1], relative[2] + center[2] };
}

void
t8_dg_rbf::interpolate (t8_mesh_deformation_node& inner_node) const
{
  t8_mesh_deformation_node inner_node_tmp = inner_node;
  if (rotation_dg_rbf) {
    /* The rotation DG-RBF returns the interpolated rotation angles as displacement. */
    t8_mesh_deformation_node angle_node = inner_node;
    rotation_dg_rbf->interpolate (angle_node);
    inner_node_tmp.position = t8_dg_rbf_rotate (inner_node.position, rotation_center, angle_node.displacement);
  }
  // Implementation for finding the Delaunay element containing the inner node
  const double tolerance = 1e-12;
  double num_steps = dgm.delaunay_elements.size ();

  for (int step = 0; step < num_steps; ++step) {
    bool found = false;
    for (size_t ielement = 0; ielement < dgm.delaunay_elements[step].size (); ++ielement) {
      const std::array<double, 3> inner_node_pos
        = { inner_node_tmp.position[0], inner_node_tmp.position[1], inner_node_tmp.position[2] };
      const auto& element = dgm.delaunay_elements[step][ielement];
      const auto& point_a = boundary_nodes[element[0]];
      const auto& point_b = boundary_nodes[element[1]];
      const auto& point_c = boundary_nodes[element[2]];
      const auto& point_d = boundary_nodes[element[3]];
      // Positions of the triangle vertices at the current substep
      const auto pos_a = dgm.t8_dgm_position_at_step (point_a, step, num_steps);
      const auto pos_b = dgm.t8_dgm_position_at_step (point_b, step, num_steps);
      const auto pos_c = dgm.t8_dgm_position_at_step (point_c, step, num_steps);
      const auto pos_d = dgm.t8_dgm_position_at_step (point_d, step, num_steps);
      // Calculate the double area of the triangle and check if it is degenerate
      const double det = dgm.t8_dgm_signed_volume (pos_a, pos_b, pos_c, pos_d, dimension);
      if (std::abs (det) < volume_tolerance) {
        /* Degenerate triangle. */
        continue;
      }
      // Compute the barycentric coordinates of the inner node with respect to the triangle (a, b, c)
      const double lambda_a = dgm.t8_dgm_signed_volume (inner_node_pos, pos_b, pos_c, pos_d, dimension) / det;
      const double lambda_b = dgm.t8_dgm_signed_volume (pos_a, inner_node_pos, pos_c, pos_d, dimension) / det;
      double lambda_d = 0.0;
      if (dimension == 3) {
        lambda_d = dgm.t8_dgm_signed_volume (pos_a, pos_b, pos_c, inner_node_pos, dimension) / det;
      }
      const double lambda_c = 1.0 - lambda_a - lambda_b - lambda_d;

      if (lambda_a >= -tolerance && lambda_b >= -tolerance && lambda_c >= -tolerance && lambda_d >= -tolerance) {
        /* RBF is used for the interpolation */
        t8_mesh_deformation_node step_node = inner_node_tmp;
        local_rbfs[step][ielement]->interpolate (step_node);
        for (int coord = 0; coord < 3; ++coord) {
          inner_node_tmp.position[coord] += step_node.displacement[coord];
        }
        found = true;
        break;
      }
    }
    SC_CHECK_ABORTF (found, "Global vertex %" T8_GLOIDX_FORMAT " is not inside any Delaunay element in substep %d.\n",
                     inner_node.global_id, step);
  }
  // After the interpolation, update the inner node's displacement to the global displacement.
  for (int coord = 0; coord < 3; ++coord) {
    inner_node.displacement[coord] = inner_node_tmp.position[coord] - inner_node.position[coord];
  }
}

void
t8_dg_rbf::solve ()
{
  // Volumes below this tolerance are roundoff, e.g. slivers spanned by coplanar boundary nodes.
  double max_coordinate = 0.0;
  for (const auto& node : boundary_nodes) {
    max_coordinate = std::max (
      { max_coordinate, std::abs (node.position[0]), std::abs (node.position[1]), std::abs (node.position[2]) });
  }
  volume_tolerance = 1e-12 * std::pow (max_coordinate, dimension);
  dgm.volume_tolerance = volume_tolerance;

  if (rotation_enabled) {
    /* The angles are interpolated with a DG-RBF without substeps, which carries the angles as displacement.
     * Moving boundary nodes get the rotation angles of the body, the fixed ones zero. */
    std::vector<t8_mesh_deformation_node> angle_nodes = boundary_nodes;
    for (auto& node : angle_nodes) {
      const bool is_moving = t8_norm (node.displacement) > 1e-12 * max_coordinate;
      node.displacement = is_moving ? rotation_angles : t8_3D_vec { 0.0, 0.0, 0.0 };
    }
    rotation_dg_rbf = std::make_unique<t8_dg_rbf> (scale_factor_support_radius, false);
    rotation_dg_rbf->set_dimension (dimension);
    rotation_dg_rbf->set_boundary_nodes (angle_nodes);
    rotation_dg_rbf->solve ();

    /* The translation is the part of the displacement which remains after the rotation of the boundary node. */
    for (size_t inode = 0; inode < boundary_nodes.size (); ++inode) {
      auto& node = boundary_nodes[inode];
      const t8_3D_vec rotated = t8_dg_rbf_rotate (node.position, rotation_center, angle_nodes[inode].displacement);
      for (int coord = 0; coord < 3; ++coord) {
        node.displacement[coord] = node.position[coord] + node.displacement[coord] - rotated[coord];
      }
      node.position = rotated;
    }
    dgm.set_boundary_nodes (boundary_nodes);
  }

  bool no_flipped_elements = false;
  double substeps = 0;
  const bool gmsh_was_initialized = gmsh::isInitialized ();
  if (!gmsh_was_initialized) {
    gmsh::initialize ();
    gmsh::option::setNumber ("General.Terminal", 0);
  }
  int max_substeps;
  if (substepping_enabled) {
    max_substeps = 16;
  }
  else {
    max_substeps = 1;
  }
  while (!no_flipped_elements && substeps < max_substeps) {
    dgm.create_delaunay_mesh (boundary_nodes, substeps);
    no_flipped_elements = dgm.check_delaunay_elements (substeps);
    substeps++;
  }
  if (substepping_enabled) {
    SC_CHECK_ABORTF (no_flipped_elements, "Delaunay graph: Flipped elements remain after the maximum of %i substeps.\n",
                     max_substeps);
    t8_debugf ("Delaunay graph: %f substeps were required to avoid flipped elements.\n", substeps - 1);
  }
  else {
    if (!no_flipped_elements) {
      t8_errorf ("Delaunay graph: Flipped elements remain after one substep. Consider enabling substepping.\n");
    }
    else {
      t8_errorf ("Delaunay graph: No flipped elements after one substep.\n");
    }
  }
  /* Build one local RBF per Delaunay element and step. The support nodes are the element's vertices. */
  const auto& elements = dgm.delaunay_elements;
  const double num_steps = elements.size ();
  local_rbfs.clear ();
  local_rbfs.resize (elements.size ());
  for (size_t step = 0; step < elements.size (); ++step) {
    local_rbfs[step].reserve (elements[step].size ());
    for (const auto& element : elements[step]) {
      std::vector<t8_mesh_deformation_node> support_nodes (dimension + 1);
      for (int ivertex = 0; ivertex <= dimension; ++ivertex) {
        const auto& node = boundary_nodes[element[ivertex]];
        const auto position = t8_dgm::t8_dgm_position_at_step (node, step, num_steps);
        support_nodes[ivertex] = node;
        for (int coord = 0; coord < 3; ++coord) {
          support_nodes[ivertex].position[coord] = position[coord];
          support_nodes[ivertex].displacement[coord] = node.displacement[coord] / num_steps;
        }
      }
      auto rbf = std::make_unique<t8_rbf> (rbf_type, scale_factor_support_radius);
      rbf->set_boundary_nodes (support_nodes);
      rbf->solve ();
      local_rbfs[step].push_back (std::move (rbf));
    }
  }

  if (!gmsh_was_initialized) {
    gmsh::finalize ();
  }
}
