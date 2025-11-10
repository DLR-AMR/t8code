
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

#include <unordered_map>
#include <t8_cmesh/t8_cmesh_vertex_connectivity/t8_cmesh_vertex_connectivity.hxx>
#include <t8_data/t8_cad.hxx>
#include <t8_cmesh.h>
#include <t8_cmesh.hxx>
#include <t8_cmesh/t8_cmesh_mesh_deformation/t8_cmesh_mesh_deformation.hxx>
#include <t8_geometry/t8_geometry_handler.hxx>
#include <t8_geometry/t8_geometry_implementations/t8_geometry_cad.hxx>
/**
 * \param [in] cmesh The coarse mesh structure.
 * \param [in] cad A pointer to the CAD-based geometry object.
 * \return A map from global vertex id to its displacement vector, so the difference between old and new coordinates of a vertex. 
 */

std::unordered_map<t8_gloidx_t, t8_3D_vec>
calculate_displacement_surface_vertices (t8_cmesh_t cmesh, const t8_cad *cad)
{
  T8_ASSERT (t8_cmesh_is_committed (cmesh));

  /* Get number of global vertices on this rank */
  //const t8_gloidx_t num_global_vertices = t8_cmesh_get_num_global_vertices (cmesh);

  const int mesh_dimension = t8_cmesh_get_dimension (cmesh);

  /* Map from global vertex id -> displacement vector. */
  std::unordered_map<t8_gloidx_t, t8_3D_vec> displacements;

  for (const auto &global_vertex : *(cmesh->vertex_connectivity)) {

    /* Get the list of all trees associated with the vertex. */
    const auto &tree_list = global_vertex.second;
    const t8_gloidx_t global_vertex_id = global_vertex.first;

    /* Get the first tree and the local corner index of the vertex in the tree. */
    const auto &first_tree = tree_list.front ();
    const t8_locidx_t first_tree_id = first_tree.first;
    const int local_corner_index = first_tree.second;

    /* Get the first tree as a reference. */
    const int *first_tree_geom_attribute = static_cast<const int *> (
      t8_cmesh_get_attribute (cmesh, t8_get_package_id (), T8_CMESH_NODE_GEOMETRY_ATTRIBUTE_KEY, first_tree_id));

    const int first_tree_entity_dim = first_tree_geom_attribute[2 * tree_list[0].second];
    const int first_tree_entity_tag = first_tree_geom_attribute[2 * tree_list[0].second + 1];

#if T8_ENABLE_DEBUG
    /* Iterate over all trees and compare to the reference tree. */
    for (const auto &[tree_id, local_corner_index] : tree_list) {

      const int *geom_attribute = static_cast<const int *> (
        t8_cmesh_get_attribute (cmesh, t8_get_package_id (), T8_CMESH_NODE_GEOMETRY_ATTRIBUTE_KEY, tree_id));

      const int entity_dim = geom_attribute[2 * local_corner_index];
      const int entity_tag = geom_attribute[2 * local_corner_index + 1];

      /* Check if the attribute of the vertex is the same in all trees. */
      if (!(entity_dim == first_tree_entity_dim && entity_tag == first_tree_entity_tag)) {
        t8_errorf ("Inconsistent entity info for global vertex %li: tree %d: dim=%d tag=%d, expected dim=%d tag=%d\n",
                   global_vertex_id, tree_id, entity_dim, entity_tag, first_tree_entity_dim, first_tree_entity_tag);
        SC_ABORTF ("Aborting due to inconsistent vertex info.\n");
      }
    }
#endif /*T8_ENABLE_DEBUG */

    /* Check if this vertex is a boundary node. */
    if (first_tree_entity_dim < mesh_dimension && first_tree_entity_dim >= 0) {

      /* Get the (u,v)-parameter of the vertex. */
      const double *uv_parameter = &((const double *) t8_cmesh_get_attribute (
        cmesh, t8_get_package_id (), T8_CMESH_NODE_PARAMETERS_ATTRIBUTE_KEY, first_tree_id))[2 * local_corner_index];

      const double *old_coords = (const double *) t8_cmesh_get_attribute (
        cmesh, t8_get_package_id (), T8_CMESH_VERTICES_ATTRIBUTE_KEY, first_tree_id);

      gp_Pnt new_coords;

      /* Find the new coordinates of the vertex in the cad file, based on the geometry its lying on. */
      switch (first_tree_entity_dim) {
      case 0: {
        new_coords = cad->t8_geom_get_cad_point (first_tree_entity_tag);
        break;
      }
      case 1: {
        Handle_Geom_Curve curve = cad->t8_geom_get_cad_curve (first_tree_entity_tag);
        curve->D0 (uv_parameter[0], new_coords);
        break;
      }
      case 2: {
        Handle_Geom_Surface surface = cad->t8_geom_get_cad_surface (first_tree_entity_tag);
        surface->D0 (uv_parameter[0], uv_parameter[1], new_coords);
        break;
      }
      default:
        SC_ABORT_NOT_REACHED ();
      }

      /* Get the old coordinates. */
      const double old_x = old_coords[3 * local_corner_index + 0];
      const double old_y = old_coords[3 * local_corner_index + 1];
      const double old_z = old_coords[3 * local_corner_index + 2];

      /* Calculate the displacement of the vertex which should be then done in the deformation. */
      displacements[global_vertex_id] = { new_coords.X () - old_x, new_coords.Y () - old_y, new_coords.Z () - old_z };
    }
  }
  return displacements;
}
/*
 * \param [in] cmesh The coarse mesh structure.
 * \param [in] a map from global vertex id to its displacement vector
 *         (difference between old and new coordinates) with length equal to the mesh dimension.
 */
void
apply_vertex_displacements (t8_cmesh_t cmesh, const std::unordered_map<t8_gloidx_t, t8_3D_vec> &displacements,
                            std::shared_ptr<t8_cad> cad)
{
  T8_ASSERT (t8_cmesh_is_committed (cmesh));

  /* Iterate over all vertices in the displacement map. */
  for (const auto &[global_vertex, displacement] : displacements) {

    /* Get the list of trees where this vertex exists. */
    const auto &tree_list = cmesh->vertex_connectivity->get_tree_list_of_vertex (global_vertex);

    /*Update the vertex coordinates in each tree. */
    for (const auto &[tree_id, local_vertex_index] : tree_list) {

      /* Get the vertex coordinates of the current tree. */
      double *tree_vertex_coords
        = (double *) t8_cmesh_get_attribute (cmesh, t8_get_package_id (), T8_CMESH_VERTICES_ATTRIBUTE_KEY, tree_id);

      /* Update the coordinates of the vertex. */
      for (int coord_index = 0; coord_index < 3; ++coord_index) {
        tree_vertex_coords[3 * local_vertex_index + coord_index] += displacement[coord_index];
      }
    }
  }

  /* Update the cad geometry. */
  t8_geometry_handler *geometry_handler = cmesh->geometry_handler;
  T8_ASSERT (geometry_handler != nullptr);

  for (auto geom = geometry_handler->begin (); geom != geometry_handler->end (); ++geom) {
    if (geom->second->t8_geom_get_type () == T8_GEOMETRY_TYPE_CAD) {
      t8_geometry_cad *cad_geom = static_cast<t8_geometry_cad *> (geom->second.get ());
      cad_geom->update_cad_manager (cad);
      break;
    }
  }
}
