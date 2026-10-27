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

/** \file t8_cmesh_mesh_deformation.cxx
 *  This file implements the routines for CAD-based mesh deformation.
 */

#include <unordered_map>
#include <set>
#include <tuple>
#include <t8_cmesh/t8_cmesh_vertex_connectivity/t8_cmesh_vertex_connectivity.hxx>
#include <t8_cad/t8_cad_handle.hxx>
#include <t8_cmesh/t8_cmesh.h>
#include <t8_cmesh/t8_cmesh.hxx>
#include <t8_cmesh/t8_cmesh_mesh_deformation/t8_cmesh_mesh_deformation.hxx>
#include <t8_geometry/t8_geometry_handler.hxx>
#include <t8_geometry/t8_geometry_implementations/t8_geometry_cad.hxx>

/**
 * Calculates the mean euclidean distance from a given vertex based on its neighbor vertices.
 * The result will later be multiplied with the scaling factor chosen by the user.
 * \param [in] cmesh The mesh.
 * \param [in] tree_list The list of trees where the vertex is present.
 * \param [in] current_global_vertex_id The global ID of the current vertex.
 * \return The mean distance to the neighbor vertices.
 */
static double
t8_calculate_mean_neighbor_distance (const t8_cmesh_t cmesh, const std::vector<std::pair<t8_locidx_t, int>> &tree_list,
                                     const t8_gloidx_t current_global_vertex_id)
{
  double total_distance = 0.0;
  /* Get the first tree where the current boundary vertex (with ID current_global_vertex_id) exists. */
  const auto &first_tree = tree_list.front ();
  const t8_locidx_t first_tree_id = first_tree.first;
  const int first_local_index = first_tree.second;

  const double *first_tree_coords = (const double *) t8_cmesh_get_attribute (
    cmesh, t8_get_package_id (), T8_CMESH_VERTICES_ATTRIBUTE_KEY, first_tree_id);

  /* Check if the coordinates are available. */
  if (first_tree_coords == nullptr) {
    t8_errorf ("Error: Coordinates attribute missing for tree %d.\n", first_tree_id);
    SC_ABORT ("Vertex coordinates are missing.");
  }

  const double current_vertex_x = first_tree_coords[3 * first_local_index + 0];
  const double current_vertex_y = first_tree_coords[3 * first_local_index + 1];
  const double current_vertex_z = first_tree_coords[3 * first_local_index + 2];

  /* Set up a set to calculate the distance to every neighbors vertex just once. */
  std::set<std::tuple<double, double, double>> unique_neighbors;

  /* Iterate over the tree list in which the current vertex is present. */
  for (const auto &[tree_id, local_vertex_index_of_the_current_vertex] : tree_list) {

    /* Get the vertex coordinates array of the current tree. */
    const double *tree_vertex_coords
      = (const double *) t8_cmesh_get_attribute (cmesh, t8_get_package_id (), T8_CMESH_VERTICES_ATTRIBUTE_KEY, tree_id);

    /* Check if the coordinates are available. */
    if (tree_vertex_coords != nullptr) {

      const int num_vertices = t8_eclass_num_vertices[t8_cmesh_get_tree_class (cmesh, tree_id)];

      for (int neighbor_vertex = 0; neighbor_vertex < num_vertices; ++neighbor_vertex) {
        /* If the vertex is the current vertex itself, we do not calculate the distance. */
        if (neighbor_vertex != local_vertex_index_of_the_current_vertex) {

          /* Get the coordinates of the neighbor vertex. */
          const double x_coord = tree_vertex_coords[3 * neighbor_vertex + 0];
          const double y_coord = tree_vertex_coords[3 * neighbor_vertex + 1];
          const double z_coord = tree_vertex_coords[3 * neighbor_vertex + 2];

          /* Save the coordinates in the set, in which double coordinates will be filtered out. */
          unique_neighbors.insert (std::make_tuple (x_coord, y_coord, z_coord));
        }
      }
    }
  }

  /* Check if the set is empty. */
  if (unique_neighbors.empty ()) {
    t8_errorf ("Error: No neighbor vertices found for global vertex %" T8_GLOIDX_FORMAT
               " to calculate local support radius.\n",
               current_global_vertex_id);
    SC_ABORT ("Calculation of local support radius failed due to missing neighbor vertices.");
  }
  /* Calculate the distance to every neighbor node. */
  for (const auto &[x_coord, y_coord, z_coord] : unique_neighbors) {
    const double distance_x = current_vertex_x - x_coord;
    const double distance_y = current_vertex_y - y_coord;
    const double distance_z = current_vertex_z - z_coord;

    total_distance += std::sqrt (distance_x * distance_x + distance_y * distance_y + distance_z * distance_z);
  }
  /* Return the average distance. */
  return total_distance / static_cast<double> (unique_neighbors.size ());
}

std::vector<t8_mesh_deformation_node>
t8_cmesh_mesh_deformation::calculate_displacement_surface_vertices (const t8_cad_handle *cad)
{
  T8_ASSERT (t8_cmesh_is_committed (associated_cmesh));

  const int mesh_dimension = t8_cmesh_get_dimension (associated_cmesh);

  /* The boundary nodes with their displacements. */
  std::vector<t8_mesh_deformation_node> boundary_node_data;

  for (const auto &global_vertex : *(associated_cmesh->vertex_connectivity)) {

    /* Get the list of all trees associated with the vertex. */
    const auto &tree_list = global_vertex.second;
    const t8_gloidx_t global_vertex_id = global_vertex.first;

    /* Get the first tree and the local corner index of the vertex in the tree. */
    const auto &first_tree = tree_list.front ();
    const t8_locidx_t first_tree_id = first_tree.first;
    const int local_corner_index = first_tree.second;

    /* Get the first tree as a reference. */
    const int *first_tree_geom_attribute = static_cast<const int *> (t8_cmesh_get_attribute (
      associated_cmesh, t8_get_package_id (), T8_CMESH_NODE_GEOMETRY_ATTRIBUTE_KEY, first_tree_id));

    /* Check if the geometry attribute is available for this tree. */
    if (first_tree_geom_attribute == nullptr) {
      t8_errorf ("Error: Geometry attribute missing for tree %d.\n", first_tree_id);
      SC_ABORT ("Geometry attribute is missing.\n");
    }

    const int first_tree_entity_dim = first_tree_geom_attribute[2 * tree_list[0].second];
    const int first_tree_entity_tag = first_tree_geom_attribute[2 * tree_list[0].second + 1];

    /* Check if all trees sharing this vertex have consistent geometry attributes. */
#if T8_ENABLE_DEBUG
    /* Iterate over all trees and compare to the reference tree. */
    for (const auto &[tree_id, local_corner_index] : tree_list) {

      const int *geom_attribute = static_cast<const int *> (
        t8_cmesh_get_attribute (associated_cmesh, t8_get_package_id (), T8_CMESH_NODE_GEOMETRY_ATTRIBUTE_KEY, tree_id));

      const int entity_dim = geom_attribute[2 * local_corner_index];
      const int entity_tag = geom_attribute[2 * local_corner_index + 1];

      /* Check if the attribute of the vertex is the same in all trees. */
      if (!(entity_dim == first_tree_entity_dim && entity_tag == first_tree_entity_tag)) {
        t8_errorf ("Error: Inconsistent entity info for global vertex %" T8_GLOIDX_FORMAT
                   ": tree %d: dim=%d tag=%d, expected dim=%d tag=%d.\n",
                   global_vertex_id, tree_id, entity_dim, entity_tag, first_tree_entity_dim, first_tree_entity_tag);
        SC_ABORT ("Inconsistency in vertex info.\n");
      }
    }
#endif /* T8_ENABLE_DEBUG */

    /* Check if this vertex is a boundary node. */
    if (first_tree_entity_dim < mesh_dimension && first_tree_entity_dim >= 0) {

      /* Get the pointer to the array of (u,v)-parameters for the CAD geometry. */
      const double *uv_attribute = (const double *) t8_cmesh_get_attribute (
        associated_cmesh, t8_get_package_id (), T8_CMESH_NODE_PARAMETERS_ATTRIBUTE_KEY, first_tree_id);

      /* Check if the (u,v)-parameters are available. */
      if (uv_attribute == nullptr) {
        t8_errorf ("Error: (u,v)-parameters are missing for tree %d.\n", first_tree_id);
        SC_ABORT ("(u,v)-parameters are missing.\n");
      }
      /* Get the (u,v)-parameter of the vertex. */
      const double *uv_parameter = &uv_attribute[2 * local_corner_index];

      /* Get the pointer to the coordinate array as it was before the deformation. */
      const double *old_coords = (const double *) t8_cmesh_get_attribute (
        associated_cmesh, t8_get_package_id (), T8_CMESH_VERTICES_ATTRIBUTE_KEY, first_tree_id);

      /* Check if the coordinates are available. */
      if (old_coords == nullptr) {
        t8_errorf ("Error: Coordinates attribute missing for tree %d.\n", first_tree_id);
        SC_ABORT ("Vertex coordinates are missing.\n");
      }

      gp_Pnt new_coords;

      /* Find the new coordinates of the vertex in the CAD file, based on the geometry it's lying on. */
      switch (first_tree_entity_dim) {
      case 0: {
        new_coords = cad->get_cad_point (first_tree_entity_tag);
        break;
      }
      case 1: {
        const Handle_Geom_Curve curve = cad->get_cad_curve (first_tree_entity_tag);
        curve->D0 (uv_parameter[0], new_coords);
        break;
      }
      case 2: {
        const Handle_Geom_Surface surface = cad->get_cad_surface (first_tree_entity_tag);
        surface->D0 (uv_parameter[0], uv_parameter[1], new_coords);
        break;
      }
      default:
        SC_ABORT_NOT_REACHED ();
      }

      /* Get the old coordinates before the deformation. */
      const double old_x = old_coords[3 * local_corner_index + 0];
      const double old_y = old_coords[3 * local_corner_index + 1];
      const double old_z = old_coords[3 * local_corner_index + 2];

      t8_mesh_deformation_node node;

      node.position = { old_x, old_y, old_z };

      /* Calculate the displacement of the vertex which should be then done in the deformation. */
      node.displacement = { new_coords.X () - old_x, new_coords.Y () - old_y, new_coords.Z () - old_z };

      /* Calculate the mean neighbor distance of the node, which can be used to calculate the radius for a compactly supported RBF. */
      node.mean_neighbor_distance = t8_calculate_mean_neighbor_distance (associated_cmesh, tree_list, global_vertex_id);

      /* Set the global id of the node. */
      node.global_id = global_vertex_id;
      boundary_node_data.push_back (node);
    }
  }
  return boundary_node_data;
}

void
t8_cmesh_mesh_deformation::apply_vertex_displacements (const std::vector<t8_mesh_deformation_node> &boundary_node_data,
                                                       std::shared_ptr<t8_cad_handle> cad)
{
  T8_ASSERT (t8_cmesh_is_committed (associated_cmesh));

  SC_CHECK_ABORT (deformation_method != nullptr, "No deformation method set. Call set_deformation_method first.\n");

  const int mesh_dimension = t8_cmesh_get_dimension (associated_cmesh);

  std::unordered_map<t8_gloidx_t, t8_3D_vec> vertex_displacements;
  std::vector<t8_mesh_deformation_node> inner_nodes;

  /* The displacements of the boundary nodes are already known from the new CAD geometry. */
  for (const auto &boundary_node : boundary_node_data) {
    vertex_displacements[boundary_node.global_id] = boundary_node.displacement;
  }

  /* Iterate over all vertices and collect the inner nodes, whose displacement has to be interpolated. */
  for (const auto &global_vertex : *(associated_cmesh->vertex_connectivity)) {
    /* Get the global vertex ID of the current vertex. */
    const t8_gloidx_t global_vertex_id = global_vertex.first;

    /* Get the list of trees where this vertex exists. */
    const auto &tree_list = global_vertex.second;

    /* Get the first tree where this vertex exists as a reference. */
    const auto &first_tree = tree_list.front ();

    /* Get the data of the first tree. */
    const int *first_tree_geom_attribute = static_cast<const int *> (t8_cmesh_get_attribute (
      associated_cmesh, t8_get_package_id (), T8_CMESH_NODE_GEOMETRY_ATTRIBUTE_KEY, first_tree.first));

    /* Check if the geometry attribute is available for this tree. */
    if (first_tree_geom_attribute == nullptr) {
      t8_errorf ("Error: Geometry attribute missing for tree %d.\n", first_tree.first);
      SC_ABORT ("Geometry attribute is missing.");
    }

    const int first_tree_entity_dim = first_tree_geom_attribute[2 * tree_list[0].second];

    /* Boundary nodes were already handled above. */
    if (first_tree_entity_dim < mesh_dimension && first_tree_entity_dim >= 0) {
      continue;
    }

    /* If it is an inner node, we do not know the displacement right away and need to interpolate to get the new coordinates. */
    const double *vertex_coords = static_cast<const double *> (t8_cmesh_get_attribute (
      associated_cmesh, t8_get_package_id (), T8_CMESH_VERTICES_ATTRIBUTE_KEY, first_tree.first));
    /* Check if the coordinates are available. */
    if (vertex_coords == nullptr) {
      t8_errorf ("Error: Coordinates attribute missing for tree %d.\n", first_tree.first);
      SC_ABORT ("Vertex coordinates are missing.");
    }

    /* Save the initial coordinates of the inner node. */
    t8_mesh_deformation_node inner_node;
    inner_node.global_id = global_vertex_id;
    inner_node.position = { vertex_coords[3 * first_tree.second], vertex_coords[3 * first_tree.second + 1],
                            vertex_coords[3 * first_tree.second + 2] };
    inner_nodes.push_back (inner_node);
  }

  /* Calculate the weights of the deformation method from the known boundary displacements. */
  deformation_method->set_boundary_nodes (boundary_node_data);
  deformation_method->solve ();

  /* Interpolate the displacement of the inner nodes. */
  for (auto &inner_node : inner_nodes) {
    deformation_method->interpolate (inner_node);
    vertex_displacements[inner_node.global_id] = inner_node.displacement;
  }

  for (const auto &global_vertex : *(associated_cmesh->vertex_connectivity)) {
    /* Get the global vertex ID of the current vertex. */
    const t8_gloidx_t global_vertex_id = global_vertex.first;
    /* Get the list of trees where this vertex exists. */
    const auto &tree_list = global_vertex.second;
    /* Update the vertex coordinates in each tree. */
    for (const auto &[tree_id, local_vertex_index] : tree_list) {
      /* Get the vertex coordinates of the current tree. */
      double *tree_vertex_coords = (double *) t8_cmesh_get_attribute (associated_cmesh, t8_get_package_id (),
                                                                      T8_CMESH_VERTICES_ATTRIBUTE_KEY, tree_id);

      /* Check if the coordinates are available. */
      if (tree_vertex_coords != nullptr) {
        /* Update the coordinates of the vertex. */
        for (int coord_index = 0; coord_index < 3; ++coord_index) {
          tree_vertex_coords[3 * local_vertex_index + coord_index]
            += vertex_displacements.at (global_vertex_id)[coord_index];
        }
      }
    }
  }

  /* Update the cad geometry. */
  detail::t8_geometry_handler *geometry_handler = detail::t8_cmesh_get_geometry_handler (associated_cmesh);
  T8_ASSERT (geometry_handler != nullptr);

  for (auto geom = geometry_handler->begin (); geom != geometry_handler->end (); ++geom) {
    if (geom->second->t8_geom_get_type () == T8_GEOMETRY_TYPE_CAD) {
      t8_geometry_cad *cad_geom = static_cast<t8_geometry_cad *> (geom->second.get ());
      cad_geom->update_cad_handle (cad);
      break;
    }
  }
}
