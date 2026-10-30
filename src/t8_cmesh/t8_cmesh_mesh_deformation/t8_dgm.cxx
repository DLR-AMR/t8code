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
#include <map>

#if T8_ENABLE_GMSH
#include <gmsh.h>
#endif

void
t8_dgm::solve ()
{
  create_delaunay_mesh (boundary_nodes);
}

void
t8_dgm::create_delaunay_mesh (const std::vector<t8_mesh_deformation_node>& boundary_node_data)
{
#if T8_ENABLE_GMSH
  const size_t num_nodes = boundary_node_data.size ();

  /* gmsh triangulates in the (u, v) plane, nodes are added with full 3D coordinates. */
  std::vector<double> coordinates_2d;
  std::vector<double> coordinates_3d;
  std::vector<std::size_t> node_tags;
  for (size_t inode = 0; inode < num_nodes; ++inode) {
    const auto& position = boundary_node_data[inode].position;
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

  /* The returned triangle node tags refer to the input points and start at 1. */
  std::vector<std::size_t> triangles;
  gmsh::algorithm::triangulate (coordinates_2d, triangles);

  /* Store the triangles as 0-based indices into the boundary node list. */
  const size_t num_triangles = triangles.size () / 3;
  delaunay_elements.clear ();
  for (size_t itri = 0; itri < num_triangles; ++itri) {
    delaunay_elements.push_back ({ triangles[3 * itri] - 1, triangles[3 * itri + 1] - 1, triangles[3 * itri + 2] - 1 });
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

  gmsh::write ("t8_dgm.msh");
  gmsh::finalize ();
#else
  SC_ABORTF ("t8code was not configured with gmsh. The Delaunay graph mapping requires gmsh.\n");
#endif
}

void
t8_dgm::interpolate (t8_mesh_deformation_node& inner_node) const
{
  // Implementation for finding the Delaunay element containing the inner node
  const double tolerance = 1e-12;
  const double inner_node_x = inner_node.position[0];
  const double inner_node_y = inner_node.position[1];
  bool found = false;

  for (const auto& element : delaunay_elements) {
    const auto& point_a = boundary_nodes[element[0]];
    const auto& point_b = boundary_nodes[element[1]];
    const auto& point_c = boundary_nodes[element[2]];
    // Calculate the double area of the triangle and check if it is degenerate
    const double det = (point_b.position[1] - point_c.position[1]) * (point_a.position[0] - point_c.position[0])
                       - (point_b.position[0] - point_c.position[0]) * (point_a.position[1] - point_c.position[1]);
    if (std::abs (det) < std::numeric_limits<double>::min ()) {
      /* Degenerate triangle. */
      continue;
    }
    // Compute the barycentric coordinates of the point (px, py) with respect to the triangle (a, b, c)
    const double lambda_a = ((point_b.position[1] - point_c.position[1]) * (inner_node_x - point_c.position[0])
                             - (point_b.position[0] - point_c.position[0]) * (inner_node_y - point_c.position[1]))
                            / det;
    const double lambda_b = ((point_c.position[1] - point_a.position[1]) * (inner_node_x - point_a.position[0])
                             - (point_c.position[0] - point_a.position[0]) * (inner_node_y - point_a.position[1]))
                            / det;
    const double lambda_c = 1.0 - lambda_a - lambda_b;

    if (lambda_a >= -tolerance && lambda_b >= -tolerance && lambda_c >= -tolerance) {
      /* The barycentric coordinates are the interpolation weights, so move the node directly. */
      for (int coord = 0; coord < 3; ++coord) {
        inner_node.displacement[coord] = lambda_a * point_a.displacement[coord] + lambda_b * point_b.displacement[coord]
                                         + lambda_c * point_c.displacement[coord];
      }
      found = true;
      break;
    }
  }
}
