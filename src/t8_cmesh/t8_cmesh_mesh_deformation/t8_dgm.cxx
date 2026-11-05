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

/** \file t8_dgm.cxx
 *  This file implements the Delaunay Graph Mapping for the mesh deformation.
 *  The Delaunay triangulation is computed with a simple 2D Bowyer-Watson algorithm.
 */

#include <t8.h>
#include <t8_cmesh/t8_cmesh_mesh_deformation/t8_dgm.hxx>
#include <t8_vector/t8_vec.hxx>
#include <map>

#if T8_ENABLE_GMSH
#include <gmsh.h>
#endif

/**
 * Computes the position of a boundary node after \a step of \a num_steps equally sized substeps.
 * \param [in] node       The boundary node.
 * \param [in] step       The number of substeps that are already performed.
 * \param [in] num_steps  The total number of substeps.
 * \return                The position of the node after \a step substeps.
 */
static inline std::array<double, 3>
t8_dgm_position_at_step (const t8_mesh_deformation_node& node, const double step, const double num_steps)
{
  return { node.position[0] + node.displacement[0] / num_steps * step,
           node.position[1] + node.displacement[1] / num_steps * step,
           node.position[2] + node.displacement[2] / num_steps * step };
}

/**
 * Computes the double signed area of the triangle (a, b, c) in the xy-plane.
 * \param [in] a, b, c  The vertices of the triangle.
 * \return              The double signed area of the triangle.
 */
static inline double
t8_dgm_signed_volume (const std::array<double, 3>& point_a, const std::array<double, 3>& point_b,
                      const std::array<double, 3>& point_c, const std::array<double, 3>& point_d, int dimension)
{
  if (dimension == 2) {
    return (point_b[1] - point_c[1]) * (point_a[0] - point_c[0])
           - (point_b[0] - point_c[0]) * (point_a[1] - point_c[1]);
  }
  else {
    /* Compute the signed volume of the tetrahedron (a, b, c, d) */
    const t8_3D_vec edge_1 = { point_b[0] - point_a[0], point_b[1] - point_a[1], point_b[2] - point_a[2] };
    const t8_3D_vec edge_2 = { point_c[0] - point_a[0], point_c[1] - point_a[1], point_c[2] - point_a[2] };
    const t8_3D_vec edge_3 = { point_d[0] - point_a[0], point_d[1] - point_a[1], point_d[2] - point_a[2] };
    t8_3D_vec cross_product;
    t8_cross_3D (edge_2, edge_3, cross_product);
    return t8_dot (edge_1, cross_product);
  }
}

void
t8_dgm::solve ()
{
  bool no_flipped_elements = false;
  double substeps = 0;
  while (!no_flipped_elements) {
    create_delaunay_mesh (boundary_nodes, substeps);
    no_flipped_elements = check_delaunay_elements (substeps);
    substeps++;
  }
  t8_debugf ("Delaunay graph mapping: %f substeps were required to avoid flipped elements.\n", substeps - 1);
}

void
t8_dgm::create_delaunay_mesh (const std::vector<t8_mesh_deformation_node>& boundary_nodes, double substeps)
{
  delaunay_elements.assign (static_cast<size_t> (substeps) + 1, {});

  for (double step = 0; step <= substeps; ++step) {
    t8_debugf ("Delaunay graph mapping: Creating Delaunay mesh for step %f of %f.\n", step, substeps);
#if T8_ENABLE_GMSH
    const size_t num_nodes = boundary_nodes.size ();

    /* gmsh triangulates in the (u, v) plane, nodes are added with full 3D coordinates. */
    std::vector<double> coordinates_2d;
    std::vector<double> coordinates_3d;
    std::vector<std::size_t> node_tags;
    for (size_t inode = 0; inode < num_nodes; ++inode) {
      const auto position = t8_dgm_position_at_step (boundary_nodes[inode], step, substeps + 1);
      coordinates_2d.push_back (position[0]);
      coordinates_2d.push_back (position[1]);
      coordinates_3d.push_back (position[0]);
      coordinates_3d.push_back (position[1]);
      coordinates_3d.push_back (position[2]);
      /* gmsh node tags start at 1. */
      node_tags.push_back (inode + 1);
    }

    gmsh::initialize ();
    gmsh::model::add ("t8_dgm");
    if (dimension == 2) {
      /* The returned triangle node tags refer to the input points and start at 1. */
      std::vector<std::size_t> triangles;
      gmsh::algorithm::triangulate (coordinates_2d, triangles);

      /* Store the triangles as 0-based indices into the boundary node list. */
      const size_t num_triangles = triangles.size () / 3;
      for (size_t itri = 0; itri < num_triangles; ++itri) {
        delaunay_elements[step].push_back (
          { triangles[3 * itri] - 1, triangles[3 * itri + 1] - 1, triangles[3 * itri + 2] - 1 });
      }

      /* Build a discrete surface holding the nodes and triangles and write it to a gmsh mesh file. */
      const int surface_tag = gmsh::model::addDiscreteEntity (2);
      gmsh::model::mesh::addNodes (2, surface_tag, node_tags, coordinates_3d);
      std::vector<std::size_t> element_tags (num_triangles);
      for (size_t itri = 0; itri < num_triangles; ++itri) {
        element_tags[itri] = itri + 1;
      }
      /* Element type 2 is the 3-node triangle. */
      gmsh::model::mesh::addElementsByType (surface_tag, 2, element_tags, triangles);
    }
    else if (dimension == 3) {
      /* The returned tetrahedron node tags refer to the input points and start at 1. */
      std::vector<std::size_t> tetrahedra;
      /* Steiner points are not boundary nodes, so the tetrahedra must only use the input points. */
      std::vector<double> steiner_points;
      gmsh::algorithm::tetrahedralize (coordinates_3d, tetrahedra, steiner_points);
      SC_CHECK_ABORTF (steiner_points.empty (),
                       "Delaunay graph mapping: gmsh added %zu Steiner points to the tetrahedralization. These are not "
                       "boundary nodes.\n",
                       steiner_points.size () / 3);

      /* Store the tetrahedra as 0-based indices into the boundary node list. */
      const size_t num_tetrahedra = tetrahedra.size () / 4;
      for (size_t itet = 0; itet < num_tetrahedra; ++itet) {
        delaunay_elements[step].push_back ({ tetrahedra[4 * itet] - 1, tetrahedra[4 * itet + 1] - 1,
                                             tetrahedra[4 * itet + 2] - 1, tetrahedra[4 * itet + 3] - 1 });
      }

      /* Build a discrete volume holding the nodes and tetrahedra and write it to a gmsh mesh file. */
      const int volume_tag = gmsh::model::addDiscreteEntity (3);
      gmsh::model::mesh::addNodes (3, volume_tag, node_tags, coordinates_3d);
      std::vector<std::size_t> element_tags (num_tetrahedra);
      for (size_t itet = 0; itet < num_tetrahedra; ++itet) {
        element_tags[itet] = itet + 1;
      }
      /* Element type 4 is the 4-node tetrahedron. */
      gmsh::model::mesh::addElementsByType (volume_tag, 4, element_tags, tetrahedra);
    }
    else {
      SC_ABORTF ("Delaunay graph mapping is only implemented for 2D and 3D meshes.\n");
    }
    gmsh::write ("t8_dgm.msh");
    gmsh::finalize ();
#else
    SC_ABORTF ("t8code was not configured with gmsh. The Delaunay graph mapping requires gmsh.\n");
#endif
  }
}
bool
t8_dgm::check_delaunay_elements (double substeps) const
{

  for (double step = 0; step <= substeps; ++step) {
    for (const auto& element : delaunay_elements[step]) {
      const auto& point_a = boundary_nodes[element[0]];
      const auto& point_b = boundary_nodes[element[1]];
      const auto& point_c = boundary_nodes[element[2]];
      const auto& point_d = boundary_nodes[element[3]];
      // Calculate the double area of the triangle before and after the current substep
      const double det_before = t8_dgm_signed_volume (t8_dgm_position_at_step (point_a, step, substeps + 1),
                                                      t8_dgm_position_at_step (point_b, step, substeps + 1),
                                                      t8_dgm_position_at_step (point_c, step, substeps + 1),
                                                      t8_dgm_position_at_step (point_d, step, substeps + 1), dimension);
      // Check the determinant after the next substep to see if the triangle has flipped (step + 1).
      const double det_after
        = t8_dgm_signed_volume (t8_dgm_position_at_step (point_a, step + 1, substeps + 1),
                                t8_dgm_position_at_step (point_b, step + 1, substeps + 1),
                                t8_dgm_position_at_step (point_c, step + 1, substeps + 1),
                                t8_dgm_position_at_step (point_d, step + 1, substeps + 1), dimension);

      if (std::abs (det_after) < std::numeric_limits<double>::min ()) {
        t8_debugf ("The element with vertices %zu, %zu, %zu is degenerate after the deformation.\n", element[0],
                   element[1], element[2]);
        return false;
      }
      if (det_before * det_after < 0.0) {
        t8_debugf ("Element %zu, %zu, %zu flipped in step %.0f of %.0f: det_before = %e, det_after = %e\n", element[0],
                   element[1], element[2], step, substeps + 1, det_before, det_after);

        return false;
      }
    }
  }
  t8_debugf ("All Delaunay triangles are valid after the deformation.\n");
  return true;
}

void
t8_dgm::interpolate (t8_mesh_deformation_node& inner_node) const
{
  t8_mesh_deformation_node inner_node_tmp = inner_node;
  // Implementation for finding the Delaunay element containing the inner node
  const double tolerance = 1e-12;
  bool found = false;
  const double num_steps = delaunay_elements.size ();
  for (double step = 0; step < num_steps; ++step) {
    for (const auto& element : delaunay_elements[step]) {
      const std::array<double, 3> inner_node_pos
        = { inner_node_tmp.position[0], inner_node_tmp.position[1], inner_node_tmp.position[2] };
      const auto& point_a = boundary_nodes[element[0]];
      const auto& point_b = boundary_nodes[element[1]];
      const auto& point_c = boundary_nodes[element[2]];
      const auto& point_d = boundary_nodes[element[3]];
      // Positions of the triangle vertices at the current substep
      const auto pos_a = t8_dgm_position_at_step (point_a, step, num_steps);
      const auto pos_b = t8_dgm_position_at_step (point_b, step, num_steps);
      const auto pos_c = t8_dgm_position_at_step (point_c, step, num_steps);
      const auto pos_d = t8_dgm_position_at_step (point_d, step, num_steps);
      // Calculate the double area of the triangle and check if it is degenerate
      const double det = t8_dgm_signed_volume (pos_a, pos_b, pos_c, pos_d, dimension);
      if (std::abs (det) < std::numeric_limits<double>::min ()) {
        /* Degenerate triangle. */
        continue;
      }
      // Compute the barycentric coordinates of the inner node with respect to the triangle (a, b, c)
      const double lambda_a = t8_dgm_signed_volume (inner_node_pos, pos_b, pos_c, pos_d, dimension) / det;
      const double lambda_b = t8_dgm_signed_volume (pos_a, inner_node_pos, pos_c, pos_d, dimension) / det;
      double lambda_d = 0.0;
      if (dimension == 3) {
        lambda_d = t8_dgm_signed_volume (pos_a, pos_b, pos_c, inner_node_pos, dimension) / det;
      }
      const double lambda_c = 1.0 - lambda_a - lambda_b - lambda_d;

      if (lambda_a >= -tolerance && lambda_b >= -tolerance && lambda_c >= -tolerance && lambda_d >= -tolerance) {
        /* The barycentric coordinates are the interpolation weights, so move the node directly. */
        for (int coord = 0; coord < 3; ++coord) {
          inner_node_tmp.displacement[coord]
            = lambda_a * point_a.displacement[coord] / num_steps + lambda_b * point_b.displacement[coord] / num_steps
              + lambda_c * point_c.displacement[coord] / num_steps + lambda_d * point_d.displacement[coord] / num_steps;
          // The new position is the old position plus the displacement.
          inner_node_tmp.position[coord] += inner_node_tmp.displacement[coord];
        }
        found = true;
        break;
      }
    }
  }
  // After the interpolation, update the inner node's displacement to the global displacement.
  for (int coord = 0; coord < 3; ++coord) {
    inner_node.displacement[coord] = inner_node_tmp.position[coord] - inner_node.position[coord];
  }
}
