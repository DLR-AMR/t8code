/*
  This file is part of t8code.
  t8code is a C library to manage a collection (a forest) of multiple
  connected adaptive space-trees of general element classes in parallel.

  Copyright (C) 2015 the developers

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

/** \file t8_cmesh_cad.cxx
 *
 * TODO: document this file
 */

#include <t8_cmesh/t8_cmesh.hxx>
#include <t8_cmesh/t8_cmesh_cad.hxx>
#include <t8_geometry/t8_geometry_implementations/t8_geometry_linear.hxx>

#if T8_ENABLE_OCC
#include <t8_geometry/t8_geometry_implementations/t8_geometry_cad.hxx>
#include <gp_Pnt.hxx>
#include <BRepBuilderAPI_MakeEdge.hxx>
#include <BRepBuilderAPI_MakeFace.hxx>
#include <BRepPrimAPI_MakePrism.hxx>
#include <BRepAlgoAPI_Fuse.hxx>
#include <TopoDS_Face.hxx>
#include <TopoDS_Edge.hxx>
#include <TopoDS.hxx>
#include <gp_Ax2.hxx>
#include <gp_Dir.hxx>
#include <gp_Circ.hxx>
#include <gp_Vec.hxx>
#include <BRep_Tool.hxx>
#include <Geom_Surface.hxx>
#endif

t8_cmesh_t
t8_cmesh_new_hollow_cylinder (sc_MPI_Comm comm, int num_tangential_trees, int num_axial_trees, int num_radial_trees,
                              int with_cad_geometry)
{
  /* Create the cmesh and the geometry */
  t8_cmesh_t cmesh;
  t8_cmesh_init (&cmesh);
  t8_cmesh_set_profiling (cmesh, 1);

  if (with_cad_geometry) {
#if T8_ENABLE_OCC
    /* Create the two cad cylinder surfaces */
    const double radius_inner = 0.25;
    const double radius_outer = 0.5;
    const gp_Pnt origin (0, 0, 0);
    const gp_Dir z_dir (0, 0, 1);
    const gp_Ax2 axis (origin, z_dir);
    const gp_Vec height (0, 0, 1);
    const gp_Circ circle_outer (axis, radius_outer);
    const gp_Circ circle_inner (axis, radius_inner);
    BRepBuilderAPI_MakeEdge make_outer_edge (circle_outer);
    const TopoDS_Edge edge_outer = make_outer_edge.Edge ();
    const TopoDS_Face face_outer = TopoDS::Face (BRepPrimAPI_MakePrism (edge_outer, height));
    const Handle_Geom_Surface cylinder_outer = BRep_Tool::Surface (face_outer);
    BRepBuilderAPI_MakeEdge make_inner_edge (circle_inner);
    const TopoDS_Edge edge_inner = make_inner_edge.Edge ();
    const TopoDS_Face face_inner = TopoDS::Face (BRepPrimAPI_MakePrism (edge_inner, height));
    const Handle_Geom_Surface cylinder_inner = BRep_Tool::Surface (face_inner);

    /* Fill a shape with cylinders and register an cad geometry with this shape. */
    TopoDS_Shape shape;
    shape = BRepBuilderAPI_MakeFace (cylinder_outer, 1e-6).Face ();
    shape = BRepAlgoAPI_Fuse (shape, BRepBuilderAPI_MakeFace (cylinder_inner, 1e-6).Face ());

    t8_cmesh_register_geometry<t8_geometry_cad> (cmesh, shape, "cad surface");

#else  /* !T8_ENABLE_OCC */
    SC_ABORTF ("OCC not linked");
#endif /* T8_ENABLE_OCC */
  }
  else {
    t8_cmesh_register_geometry<t8_geometry_linear> (cmesh);
  }

#if T8_ENABLE_OCC
  /* Save the indices of the cylinders inside the shape for later usage. 
   * The indices start with 1 and are in the same order as we put in the cylinders. */
  int cylinder_outer_index = 1, cylinder_inner_index = 2;
#endif /* T8_ENABLE_OCC */

  /* Start the construction of the actual cylindrical cmesh. We are going to use three loops
   * to iterate over the three dimensions of cylinder coordinates. */
  const double radius_outer = 0.5, radius_inner = 0.25;
  const double dr = (radius_outer - radius_inner) / num_radial_trees;
  const double dphi = 2.0 * M_PI / num_tangential_trees;
  const double dh = 1.0 / num_axial_trees;
  /* Allocate memory for saving the node coordinates 
   * and in case of usage of the cad geometry, the node parameters */
  double *vertices;
  vertices = T8_ALLOC (double, num_tangential_trees *num_axial_trees *num_radial_trees * 24);
#if T8_ENABLE_OCC
  double *parameters;
  parameters = T8_ALLOC (double, num_tangential_trees *num_axial_trees * 8);
#endif /* T8_ENABLE_OCC */

  /* Compute vertex coordinates and parameters */
  for (int i_tangential_trees = 0; i_tangential_trees < num_tangential_trees; ++i_tangential_trees) {
    for (int i_axial_trees = 0; i_axial_trees < num_axial_trees; ++i_axial_trees) {
      for (int i_radial_trees = 0; i_radial_trees < num_radial_trees; ++i_radial_trees) {
        t8_cmesh_set_tree_class (
          cmesh, (i_tangential_trees * num_axial_trees + i_axial_trees) * num_radial_trees + i_radial_trees,
          T8_ECLASS_HEX);
        const int current_tree_vertices
          = ((i_tangential_trees * num_axial_trees + i_axial_trees) * num_radial_trees + i_radial_trees) * 24;
        vertices[current_tree_vertices + 0]
          = cos ((i_tangential_trees + 1) * dphi) * (radius_inner + (i_radial_trees + 1) * dr);
        vertices[current_tree_vertices + 1]
          = sin ((i_tangential_trees + 1) * dphi) * (radius_inner + (i_radial_trees + 1) * dr);
        vertices[current_tree_vertices + 2] = -0.5 + i_axial_trees * dh;
        vertices[current_tree_vertices + 3]
          = cos ((i_tangential_trees + 1) * dphi) * (radius_inner + i_radial_trees * dr);
        vertices[current_tree_vertices + 4]
          = sin ((i_tangential_trees + 1) * dphi) * (radius_inner + i_radial_trees * dr);
        vertices[current_tree_vertices + 5] = -0.5 + i_axial_trees * dh;
        vertices[current_tree_vertices + 6]
          = cos (i_tangential_trees * dphi) * (radius_inner + (i_radial_trees + 1) * dr);
        vertices[current_tree_vertices + 7]
          = sin (i_tangential_trees * dphi) * (radius_inner + (i_radial_trees + 1) * dr);
        vertices[current_tree_vertices + 8] = -0.5 + i_axial_trees * dh;
        vertices[current_tree_vertices + 9] = cos (i_tangential_trees * dphi) * (radius_inner + i_radial_trees * dr);
        vertices[current_tree_vertices + 10] = sin (i_tangential_trees * dphi) * (radius_inner + i_radial_trees * dr);
        vertices[current_tree_vertices + 11] = -0.5 + i_axial_trees * dh;
        vertices[current_tree_vertices + 12]
          = cos ((i_tangential_trees + 1) * dphi) * (radius_inner + (i_radial_trees + 1) * dr);
        vertices[current_tree_vertices + 13]
          = sin ((i_tangential_trees + 1) * dphi) * (radius_inner + (i_radial_trees + 1) * dr);
        vertices[current_tree_vertices + 14] = -0.5 + (i_axial_trees + 1) * dh;
        vertices[current_tree_vertices + 15]
          = cos ((i_tangential_trees + 1) * dphi) * (radius_inner + i_radial_trees * dr);
        vertices[current_tree_vertices + 16]
          = sin ((i_tangential_trees + 1) * dphi) * (radius_inner + i_radial_trees * dr);
        vertices[current_tree_vertices + 17] = -0.5 + (i_axial_trees + 1) * dh;
        vertices[current_tree_vertices + 18]
          = cos (i_tangential_trees * dphi) * (radius_inner + (i_radial_trees + 1) * dr);
        vertices[current_tree_vertices + 19]
          = sin (i_tangential_trees * dphi) * (radius_inner + (i_radial_trees + 1) * dr);
        vertices[current_tree_vertices + 20] = -0.5 + (i_axial_trees + 1) * dh;
        vertices[current_tree_vertices + 21] = cos (i_tangential_trees * dphi) * (radius_inner + i_radial_trees * dr);
        vertices[current_tree_vertices + 22] = sin (i_tangential_trees * dphi) * (radius_inner + i_radial_trees * dr);
        vertices[current_tree_vertices + 23] = -0.5 + (i_axial_trees + 1) * dh;

        t8_cmesh_set_tree_vertices (
          cmesh, (i_tangential_trees * num_axial_trees + i_axial_trees) * num_radial_trees + i_radial_trees,
          vertices + ((i_tangential_trees * num_axial_trees + i_axial_trees) * num_radial_trees + i_radial_trees) * 24,
          24);

        /* Assign parameters if cad is enabled */
#if T8_ENABLE_OCC
        if (with_cad_geometry) {
          /* Calculate parameters if cell lies on boundary */
          const int current_tree_parameters = (i_tangential_trees * num_axial_trees + i_axial_trees) * 8;
          if (i_radial_trees == 0 || i_radial_trees == num_radial_trees - 1) {
            parameters[current_tree_parameters + 0] = (i_tangential_trees + 1) * dphi;
            parameters[current_tree_parameters + 1] = 0.5 - i_axial_trees * dh;
            parameters[current_tree_parameters + 2] = i_tangential_trees * dphi;
            parameters[current_tree_parameters + 3] = 0.5 - i_axial_trees * dh;
            parameters[current_tree_parameters + 4] = (i_tangential_trees + 1) * dphi;
            parameters[current_tree_parameters + 5] = 0.5 + -(i_axial_trees + 1) * dh;
            parameters[current_tree_parameters + 6] = i_tangential_trees * dphi;
            parameters[current_tree_parameters + 7] = 0.5 - (i_axial_trees + 1) * dh;
          }
          int edges[24] = { 0 };

          const int current_tree
            = (i_tangential_trees * num_axial_trees + i_axial_trees) * num_radial_trees + i_radial_trees;
          /* If geometry on both sides of cell */
          if (num_radial_trees == 1) {
            /* Assign cad geometries to the corresponding faces */
            int faces[6] = { 0 };
            faces[0] = cylinder_outer_index;
            faces[1] = cylinder_inner_index;

            /* Assign attributes to cmesh cells */
            t8_cmesh_set_attribute (cmesh, current_tree, t8_get_package_id (), T8_CMESH_CAD_FACE_ATTRIBUTE_KEY, faces,
                                    6 * sizeof (int), 0);
            t8_cmesh_set_attribute (cmesh, current_tree, t8_get_package_id (), T8_CMESH_CAD_EDGE_ATTRIBUTE_KEY, edges,
                                    24 * sizeof (int), 0);
            t8_cmesh_set_attribute (cmesh, current_tree, t8_get_package_id (),
                                    T8_CMESH_CAD_FACE_PARAMETERS_ATTRIBUTE_KEY + 0,
                                    parameters + current_tree_parameters, 8 * sizeof (double), 1);
            t8_cmesh_set_attribute (cmesh, current_tree, t8_get_package_id (),
                                    T8_CMESH_CAD_FACE_PARAMETERS_ATTRIBUTE_KEY + 1,
                                    parameters + current_tree_parameters, 8 * sizeof (double), 1);
          }
          /* If geometry only on face 1 */
          else if (i_radial_trees == 0) {
            /* Assign cad geometries to the corresponding faces */
            int faces[6] = { 0 };
            faces[1] = cylinder_inner_index;

            /* Assign attributes to cmesh cells */
            t8_cmesh_set_attribute (cmesh, current_tree, t8_get_package_id (), T8_CMESH_CAD_FACE_ATTRIBUTE_KEY, faces,
                                    6 * sizeof (int), 0);
            t8_cmesh_set_attribute (cmesh, current_tree, t8_get_package_id (), T8_CMESH_CAD_EDGE_ATTRIBUTE_KEY, edges,
                                    24 * sizeof (int), 0);
            t8_cmesh_set_attribute (cmesh, current_tree, t8_get_package_id (),
                                    T8_CMESH_CAD_FACE_PARAMETERS_ATTRIBUTE_KEY + 1,
                                    parameters + current_tree_parameters, 8 * sizeof (double), 1);
          }
          /* If geometry only on face 0 */
          else if (i_radial_trees == num_radial_trees - 1) {
            /* Assign cad geometries to the corresponding faces */
            int faces[6] = { 0 };
            faces[0] = cylinder_outer_index;

            /* Assign attributes to cmesh cells */
            t8_cmesh_set_attribute (cmesh, current_tree, t8_get_package_id (), T8_CMESH_CAD_FACE_ATTRIBUTE_KEY, faces,
                                    6 * sizeof (int), 0);
            t8_cmesh_set_attribute (cmesh, current_tree, t8_get_package_id (), T8_CMESH_CAD_EDGE_ATTRIBUTE_KEY, edges,
                                    24 * sizeof (int), 0);
            t8_cmesh_set_attribute (cmesh, current_tree, t8_get_package_id (),
                                    T8_CMESH_CAD_FACE_PARAMETERS_ATTRIBUTE_KEY + 0,
                                    parameters + current_tree_parameters, 8 * sizeof (double), 1);
          }
          /* If there is no geometry */
          else {
            /* Assign cad geometries to the corresponding faces */
            int faces[6] = { 0 };

            /* Assign attributes to cmesh cells */
            t8_cmesh_set_attribute (cmesh, current_tree, t8_get_package_id (), T8_CMESH_CAD_FACE_ATTRIBUTE_KEY, faces,
                                    6 * sizeof (int), 0);
            t8_cmesh_set_attribute (cmesh, current_tree, t8_get_package_id (), T8_CMESH_CAD_EDGE_ATTRIBUTE_KEY, edges,
                                    24 * sizeof (int), 0);
          }
        }
#endif /* T8_ENABLE_OCC */
        /* Join radial neighbors */
        if (i_radial_trees > 0) {
          t8_cmesh_set_join (
            cmesh, (i_tangential_trees * num_axial_trees + i_axial_trees) * num_radial_trees + i_radial_trees - 1,
            (i_tangential_trees * num_axial_trees + i_axial_trees) * num_radial_trees + i_radial_trees, 0, 1, 0);
        }
      }

      /* Join axial neighbors */
      if (i_axial_trees > 0) {
        for (int i_radial_trees = 0; i_radial_trees < num_radial_trees; ++i_radial_trees) {
          t8_cmesh_set_join (
            cmesh, (i_tangential_trees * num_axial_trees + i_axial_trees - 1) * num_radial_trees + i_radial_trees,
            (i_tangential_trees * num_axial_trees + i_axial_trees) * num_radial_trees + i_radial_trees, 5, 4, 0);
        }
      }
    }
  }

  /* Join tangential neighbors of seam */
  for (int i_axial_trees = 0; i_axial_trees < num_axial_trees; ++i_axial_trees) {
    for (int i_radial_trees = 0; i_radial_trees < num_radial_trees; ++i_radial_trees) {
      t8_cmesh_set_join (
        cmesh, (num_axial_trees + i_axial_trees) * num_radial_trees + i_radial_trees,
        ((num_tangential_trees - 1) * num_axial_trees + i_axial_trees) * num_radial_trees + i_radial_trees, 2, 3, 0);
    }
  }

  /* Join all other tangential neighbors */
  for (int i_tangential_trees = 1; i_tangential_trees < num_tangential_trees; ++i_tangential_trees) {
    for (int i_axial_trees = 0; i_axial_trees < num_axial_trees; ++i_axial_trees) {
      for (int i_radial_trees = 0; i_radial_trees < num_radial_trees; ++i_radial_trees) {
        t8_cmesh_set_join (
          cmesh, ((i_tangential_trees - 1) * num_axial_trees + i_axial_trees) * num_radial_trees + i_radial_trees,
          (i_tangential_trees * num_axial_trees + i_axial_trees) * num_radial_trees + i_radial_trees, 2, 3, 0);
      }
    }
  }

  /* Commit the cmesh and free allocated memory. */
  t8_cmesh_commit (cmesh, comm);
  T8_FREE (vertices);
#if T8_ENABLE_OCC
  T8_FREE (parameters);
#endif /* T8_ENABLE_OCC */
  return cmesh;
}
