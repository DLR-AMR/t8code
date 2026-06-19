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

static double
calculate_local_support_radius (t8_cmesh_t cmesh, const std::vector<std::pair<t8_locidx_t, int>> &tree_list,
                                int current_global_vertex_id)
{
  double total_distance = 0.0;

  /** Get the first tree where this vertex exists. */
  const auto &first_tree = tree_list.front ();
  const t8_locidx_t first_tree_id = first_tree.first;
  const int first_local_index = first_tree.second;

  const double *first_tree_coords = (const double *) t8_cmesh_get_attribute (
    cmesh, t8_get_package_id (), T8_CMESH_VERTICES_ATTRIBUTE_KEY, first_tree_id);

  /** Check if the coordinates are available. */
  if (first_tree_coords == nullptr) {
    t8_errorf ("Error: Coordinates attribute missing for tree %d\n.", first_tree_id);
    SC_ABORTF ("Vertex coordinates are missing.");
  }

  const double current_vertex_x = first_tree_coords[3 * first_local_index + 0];
  const double current_vertex_y = first_tree_coords[3 * first_local_index + 1];
  const double current_vertex_z = first_tree_coords[3 * first_local_index + 2];

  /** Set up a set to calculate the distance to every neighbors vertex just once. */
  std::set<std::tuple<double, double, double>> unique_neighbors;

  /** Iterate over the tree list in which the current vertex is present. */
  for (const auto &[tree_id, local_vertex_index_of_the_current_vertex] : tree_list) {

    /** Get the vertex coordinates array of the current tree. */
    const double *tree_vertex_coords
      = (const double *) t8_cmesh_get_attribute (cmesh, t8_get_package_id (), T8_CMESH_VERTICES_ATTRIBUTE_KEY, tree_id);

    /** Check if the coordinates are available. */
    if (tree_vertex_coords != nullptr) {

      int num_vertices = t8_eclass_num_vertices[t8_cmesh_get_tree_class (cmesh, tree_id)];

      for (int neighbor_vertex = 0; neighbor_vertex < num_vertices; ++neighbor_vertex) {
        /** If the vertex is the current vertex itself, we do not calculate the distance. */
        if (neighbor_vertex != local_vertex_index_of_the_current_vertex) {

          /** Get the coordinates of the neighbor vertex. */
          double x = tree_vertex_coords[3 * neighbor_vertex + 0];
          double y = tree_vertex_coords[3 * neighbor_vertex + 1];
          double z = tree_vertex_coords[3 * neighbor_vertex + 2];

          /** Save the coordinates in the set, in which double coordinates will be filtered out. */
          unique_neighbors.insert (std::make_tuple (x, y, z));
        }
      }
    }
  }

  /** Check if the set is empty. */
  if (unique_neighbors.empty ()) {
    t8_errorf ("Error: No neighbor vertices found for global vertex %d to calculate local support radius.\n",
               current_global_vertex_id);
    SC_ABORTF ("Calculation of local support radius failed due to missing neighbor vertices.");
  }
  /** Calculate the distance to every neighbor node. */
  for (const auto &[x, y, z] : unique_neighbors) {
    double distance_x = current_vertex_x - x;
    double distance_y = current_vertex_y - y;
    double distance_z = current_vertex_z - z;

    total_distance += std::sqrt (distance_x * distance_x + distance_y * distance_y + distance_z * distance_z);
  }
  /** To get the average, we divide through the amount of neighbor nodes. */
  return (total_distance / static_cast<double> (unique_neighbors.size ()));
}

std::unordered_map<t8_gloidx_t, t8_rbf_boundary_node>
t8_cmesh_mesh_deformation::calculate_displacement_surface_vertices (const t8_cad_handle *cad,
                                                                    const t8_rbf_function_type rbf_type,
                                                                    const double scale_factor_support_radius)
{
  T8_ASSERT (t8_cmesh_is_committed (associated_cmesh));

  const int mesh_dimension = t8_cmesh_get_dimension (associated_cmesh);

  /* Map from global vertex id -> displacement vector. */
  std::unordered_map<t8_gloidx_t, t8_rbf_boundary_node> boundary_node_data;

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
      t8_errorf ("Error: Geometry attribute missing for tree %d\n.", first_tree_id);
      SC_ABORTF ("Geometry attribute is missing.\n");
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
        t8_errorf (
          "Error: Inconsistent entity info for global vertex %li: tree %d: dim=%d tag=%d, expected dim=%d tag=%d\n",
          global_vertex_id, tree_id, entity_dim, entity_tag, first_tree_entity_dim, first_tree_entity_tag);
        SC_ABORTF ("Inconsistency in vertex info.\n");
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
        t8_errorf ("Error: (u,v)-parameters are missing for tree %d\n.", first_tree_id);
        SC_ABORT ("(u,v)-parameters are missing.\n");
      }
      /* Get the (u,v)-parameter of the vertex. */
      const double *uv_parameter = &uv_attribute[2 * local_corner_index];

      /* Get the pointer to the coordinate array as it was before the deformation. */
      const double *old_coords = (const double *) t8_cmesh_get_attribute (
        associated_cmesh, t8_get_package_id (), T8_CMESH_VERTICES_ATTRIBUTE_KEY, first_tree_id);

      /* Check if the coordinates are available. */
      if (old_coords == nullptr) {
        t8_errorf ("Error: Coordinates attribute missing for tree %d\n.", first_tree_id);
        SC_ABORTF ("Vertex coordinates are missing.\n");
      }

      gp_Pnt new_coords;

      /* Find the new coordinates of the vertex in the CAD file, based on the geometry its lying on. */
      switch (first_tree_entity_dim) {
      case 0: {
        new_coords = cad->get_cad_point (first_tree_entity_tag);
        break;
      }
      case 1: {
        Handle_Geom_Curve curve = cad->get_cad_curve (first_tree_entity_tag);
        curve->D0 (uv_parameter[0], new_coords);
        break;
      }
      case 2: {
        Handle_Geom_Surface surface = cad->get_cad_surface (first_tree_entity_tag);
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

      t8_rbf_boundary_node node;

      node.position = { old_x, old_y, old_z };

      /* Calculate the displacement of the vertex which should be then done in the deformation. */
      node.displacement = { new_coords.X () - old_x, new_coords.Y () - old_y, new_coords.Z () - old_z };

      node.weight.fill (0.0);

      if (rbf_type == T8_RBF_CP_C2) {
        node.local_support_radius = scale_factor_support_radius
                                    * calculate_local_support_radius (associated_cmesh, tree_list, global_vertex_id);
      }
      boundary_node_data[global_vertex_id] = node;
    }
  }
  return boundary_node_data;
}

void
t8_cmesh_mesh_deformation::apply_vertex_displacements (
  std::unordered_map<t8_gloidx_t, t8_rbf_boundary_node> &boundary_node_data, std::shared_ptr<t8_cad_handle> cad,
  const t8_rbf_function_type rbf_type)
{
  T8_ASSERT (t8_cmesh_is_committed (associated_cmesh));

  t8_rbf rbf_handler (rbf_type);

  rbf_handler.set_boundary_nodes (std::move (boundary_node_data));

  /** Calculate the weights. */
  rbf_handler.solve ();

  /** Iterate over all vertices in the displacement map. */
  for (const auto &global_vertex : *(associated_cmesh->vertex_connectivity)) {
    /** Get the global vertex ID of the current vertex. */
    t8_gloidx_t global_vertex_id = global_vertex.first;

    /** Get the list of trees where this vertex exists. */
    const auto &tree_list = global_vertex.second;

    /** Get the first tree where this vertex exists as a reference. */
    const auto &first_tree = tree_list.front ();

    /* Get the data of the first tree. */
    const int *first_tree_geom_attribute = static_cast<const int *> (t8_cmesh_get_attribute (
      associated_cmesh, t8_get_package_id (), T8_CMESH_NODE_GEOMETRY_ATTRIBUTE_KEY, first_tree.first));

    /* Check if the geometry attribute is available for this tree. */
    if (first_tree_geom_attribute == nullptr) {
      t8_errorf ("Error: Geometry attribute missing for tree %d\n.", first_tree.first);
      SC_ABORTF ("Geometry attribute is missing.");
    }

    const int first_tree_entity_dim = first_tree_geom_attribute[2 * tree_list[0].second];

    const int mesh_dimension = t8_cmesh_get_dimension (associated_cmesh);

    /** Contains the displacement whether its a boundary node or a inner node. */
    t8_3D_vec final_displacement;

    /* Check if this vertex is a boundary node. 
      If so, we can use the already known displacement extracted from the new CAD geometry input given.*/
    if (first_tree_entity_dim < mesh_dimension && first_tree_entity_dim >= 0) {

      final_displacement = rbf_handler.get_boundary_displacement (global_vertex_id);
    }
    /** If it is a inner node we don not know the displacement right away and need to interpolate to get the new coordinates. */
    else {

      double *vertex_coords = (double *) t8_cmesh_get_attribute (associated_cmesh, t8_get_package_id (),
                                                                 T8_CMESH_VERTICES_ATTRIBUTE_KEY, first_tree.first);
      /** Check if the coordinates are available. */
      if (vertex_coords == nullptr) {
        t8_errorf ("Error: Coordinates attribute missing for tree %d\n.", first_tree.first);
        SC_ABORTF ("Vertex coordinates are missing.");
      }

      /** Save the initial coordinates of the inner node. */
      t8_rbf_node inner_node;
      inner_node.global_id = global_vertex_id;
      inner_node.position = { vertex_coords[3 * first_tree.second], vertex_coords[3 * first_tree.second + 1],
                              vertex_coords[3 * first_tree.second + 2] };

      /** Calculate the displacement of the inner node. */
      rbf_handler.interpolate (inner_node);
      final_displacement = inner_node.displacement;
    }

    /** Update the vertex coordinates in each tree. */
    for (const auto &[tree_id, local_vertex_index] : tree_list) {
      /** Get the vertex coordinates of the current tree. */
      double *tree_vertex_coords = (double *) t8_cmesh_get_attribute (associated_cmesh, t8_get_package_id (),
                                                                      T8_CMESH_VERTICES_ATTRIBUTE_KEY, tree_id);

      /** Check if the coordinates are available. */
      if (tree_vertex_coords != nullptr) {
        /** Update the coordinates of the vertex. */
        for (int coord_index = 0; coord_index < 3; ++coord_index) {
          tree_vertex_coords[3 * local_vertex_index + coord_index] += final_displacement[coord_index];
        }
      }
    }
  }

  /** Update the cad geometry. */
  t8_geometry_handler *geometry_handler = associated_cmesh->geometry_handler;
  T8_ASSERT (geometry_handler != nullptr);

  for (auto geom = geometry_handler->begin (); geom != geometry_handler->end (); ++geom) {
    if (geom->second->t8_geom_get_type () == T8_GEOMETRY_TYPE_CAD) {
      t8_geometry_cad *cad_geom = static_cast<t8_geometry_cad *> (geom->second.get ());
      cad_geom->update_cad_handle (cad);
      break;
    }
  }
}
