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

/** \file t8_cmesh_occ.cxx
 *
 * TODO: document this file
 */

#include <t8_cmesh/t8_cmesh_occ.hxx>
#include <t8_geometry/t8_geometry_implementations/t8_geometry_linear.h>
#include <t8_geometry/t8_geometry_implementations/t8_geometry_occ.hxx>
#include <t8_cmesh_vtk.h>

#if T8_WITH_OCC
#include <gp_Pnt.hxx>
#include <BRepBuilderAPI_MakeEdge.hxx>
#include <BRepPrimAPI_MakePrism.hxx>
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
t8_cmesh_new_hollow_cylinder (sc_MPI_Comm comm, int num_tangential_trees, 
                              int num_axial_trees, int num_radial_trees, 
                              int with_occ_geometry)
{  
  t8_cmesh_t cmesh;
  t8_cmesh_init (&cmesh);
  t8_cmesh_set_profiling(cmesh, 1);
  t8_geometry_occ *geometry_occ = new t8_geometry_occ (3, "occ surface dim=3");
  t8_geometry_c *geometry_linear = t8_geometry_linear_new (3);
  int cylinder_outer_index, cylinder_inner_index;
  
  if (with_occ_geometry)
  {
    #if T8_WITH_OCC
    /* Create occ cylinder surfaces */
    const double radius_inner = 0.25;
    const double radius_outer = 0.5;
    const gp_Pnt origin(0, 0, 0);
    const gp_Dir z_dir(0, 0, 1);
    const gp_Ax2 axis(origin, z_dir);
    const gp_Vec height(0, 0, 1);
    const gp_Circ circle_outer(axis, radius_outer);
    const gp_Circ circle_inner(axis, radius_inner);
    BRepBuilderAPI_MakeEdge make_outer_edge(circle_outer);
    const TopoDS_Edge edge_outer = make_outer_edge.Edge();
    const TopoDS_Face face_outer = TopoDS::Face(BRepPrimAPI_MakePrism(edge_outer, height));
    const Handle_Geom_Surface cylinder_outer = BRep_Tool::Surface(face_outer);
    BRepBuilderAPI_MakeEdge make_inner_edge(circle_inner);
    const TopoDS_Edge edge_inner = make_inner_edge.Edge();
    const TopoDS_Face face_inner = TopoDS::Face(BRepPrimAPI_MakePrism(edge_inner, height));
    const Handle_Geom_Surface cylinder_inner = BRep_Tool::Surface(face_inner);
    cylinder_outer_index = geometry_occ->t8_geom_push_occ_surface(cylinder_outer);
    cylinder_inner_index = geometry_occ->t8_geom_push_occ_surface(cylinder_inner);
    #else /* !T8_WITH_OCC */
    SC_ABORTF("OCC not linked");
    #endif /* T8_WITH_OCC */ 
  }
  
  double *vertices, *parameters;
  const double radius_outer = 0.5, radius_inner = 0.25;
  const double dr = (radius_outer - radius_inner) / num_radial_trees;
  const double dphi = 2.0 * M_PI / num_tangential_trees;
  const double dh = 1.0 / num_axial_trees;
  vertices = T8_ALLOC(double, num_tangential_trees * num_axial_trees * num_radial_trees * 24);
  parameters = T8_ALLOC(double, num_tangential_trees * num_axial_trees * 8);
  
  /* Compute vertex coordinates and parameters */
  for (int i_tangential_trees = 0; i_tangential_trees < num_tangential_trees; ++i_tangential_trees)
  {
    for (int i_axial_trees = 0; i_axial_trees < num_axial_trees; ++i_axial_trees)
    {
      for (int i_radial_trees = 0; i_radial_trees < num_radial_trees; ++i_radial_trees)
      {
        t8_cmesh_set_tree_class (cmesh, (i_tangential_trees * num_axial_trees + i_axial_trees) * num_radial_trees + i_radial_trees, T8_ECLASS_HEX);
        vertices[((i_tangential_trees * num_axial_trees + i_axial_trees) * num_radial_trees + i_radial_trees) * 24 + 0] = cos((i_tangential_trees + 1) * dphi) * (radius_inner + (i_radial_trees + 1) * dr);
        vertices[((i_tangential_trees * num_axial_trees + i_axial_trees) * num_radial_trees + i_radial_trees) * 24 + 1] = sin((i_tangential_trees + 1) * dphi) * (radius_inner + (i_radial_trees + 1) * dr);
        vertices[((i_tangential_trees * num_axial_trees + i_axial_trees) * num_radial_trees + i_radial_trees) * 24 + 2] = -0.5 + i_axial_trees * dh;
        vertices[((i_tangential_trees * num_axial_trees + i_axial_trees) * num_radial_trees + i_radial_trees) * 24 + 3] = cos((i_tangential_trees + 1) * dphi) * (radius_inner + i_radial_trees * dr);
        vertices[((i_tangential_trees * num_axial_trees + i_axial_trees) * num_radial_trees + i_radial_trees) * 24 + 4] = sin((i_tangential_trees + 1) * dphi) * (radius_inner + i_radial_trees * dr);
        vertices[((i_tangential_trees * num_axial_trees + i_axial_trees) * num_radial_trees + i_radial_trees) * 24 + 5] = -0.5 + i_axial_trees * dh;
        vertices[((i_tangential_trees * num_axial_trees + i_axial_trees) * num_radial_trees + i_radial_trees) * 24 + 6] = cos(i_tangential_trees * dphi) * (radius_inner + (i_radial_trees + 1) * dr);
        vertices[((i_tangential_trees * num_axial_trees + i_axial_trees) * num_radial_trees + i_radial_trees) * 24 + 7] = sin(i_tangential_trees * dphi) * (radius_inner + (i_radial_trees + 1) * dr);
        vertices[((i_tangential_trees * num_axial_trees + i_axial_trees) * num_radial_trees + i_radial_trees) * 24 + 8] = -0.5 + i_axial_trees * dh;
        vertices[((i_tangential_trees * num_axial_trees + i_axial_trees) * num_radial_trees + i_radial_trees) * 24 + 9] = cos(i_tangential_trees * dphi) * (radius_inner + i_radial_trees * dr);
        vertices[((i_tangential_trees * num_axial_trees + i_axial_trees) * num_radial_trees + i_radial_trees) * 24 + 10] = sin(i_tangential_trees * dphi) * (radius_inner + i_radial_trees * dr);
        vertices[((i_tangential_trees * num_axial_trees + i_axial_trees) * num_radial_trees + i_radial_trees) * 24 + 11] = -0.5 + i_axial_trees * dh;
        vertices[((i_tangential_trees * num_axial_trees + i_axial_trees) * num_radial_trees + i_radial_trees) * 24 + 12] = cos((i_tangential_trees + 1) * dphi) * (radius_inner + (i_radial_trees + 1) * dr);
        vertices[((i_tangential_trees * num_axial_trees + i_axial_trees) * num_radial_trees + i_radial_trees) * 24 + 13] = sin((i_tangential_trees + 1) * dphi) * (radius_inner + (i_radial_trees + 1) * dr);
        vertices[((i_tangential_trees * num_axial_trees + i_axial_trees) * num_radial_trees + i_radial_trees) * 24 + 14] = -0.5 + (i_axial_trees + 1) * dh;
        vertices[((i_tangential_trees * num_axial_trees + i_axial_trees) * num_radial_trees + i_radial_trees) * 24 + 15] = cos((i_tangential_trees + 1) * dphi) * (radius_inner + i_radial_trees * dr);
        vertices[((i_tangential_trees * num_axial_trees + i_axial_trees) * num_radial_trees + i_radial_trees) * 24 + 16] = sin((i_tangential_trees + 1) * dphi) * (radius_inner + i_radial_trees * dr);
        vertices[((i_tangential_trees * num_axial_trees + i_axial_trees) * num_radial_trees + i_radial_trees) * 24 + 17] = -0.5 + (i_axial_trees + 1) * dh;
        vertices[((i_tangential_trees * num_axial_trees + i_axial_trees) * num_radial_trees + i_radial_trees) * 24 + 18] = cos(i_tangential_trees * dphi) * (radius_inner + (i_radial_trees + 1) * dr);
        vertices[((i_tangential_trees * num_axial_trees + i_axial_trees) * num_radial_trees + i_radial_trees) * 24 + 19] = sin(i_tangential_trees * dphi) * (radius_inner + (i_radial_trees + 1) * dr);
        vertices[((i_tangential_trees * num_axial_trees + i_axial_trees) * num_radial_trees + i_radial_trees) * 24 + 20] = -0.5 + (i_axial_trees + 1) * dh;
        vertices[((i_tangential_trees * num_axial_trees + i_axial_trees) * num_radial_trees + i_radial_trees) * 24 + 21] = cos(i_tangential_trees * dphi) * (radius_inner + i_radial_trees * dr);
        vertices[((i_tangential_trees * num_axial_trees + i_axial_trees) * num_radial_trees + i_radial_trees) * 24 + 22] = sin(i_tangential_trees * dphi) * (radius_inner + i_radial_trees * dr);
        vertices[((i_tangential_trees * num_axial_trees + i_axial_trees) * num_radial_trees + i_radial_trees) * 24 + 23] = -0.5 + (i_axial_trees + 1) * dh;
        t8_cmesh_set_tree_vertices (cmesh, (i_tangential_trees * num_axial_trees + i_axial_trees) * num_radial_trees + i_radial_trees, vertices + ((i_tangential_trees * num_axial_trees + i_axial_trees) * num_radial_trees + i_radial_trees) * 24, 24);

        /* Assign parameters if occ is enabled */
        if (with_occ_geometry)
        {
          /* Calculate parameters if cell lies on boundary */
          if (i_radial_trees == 0 || i_radial_trees == num_radial_trees - 1)
          {
            parameters[(i_tangential_trees * num_axial_trees + i_axial_trees) * 8 + 0] = (i_tangential_trees + 1) * dphi;
            parameters[(i_tangential_trees * num_axial_trees + i_axial_trees) * 8 + 1] = 0.5 - i_axial_trees * dh;
            parameters[(i_tangential_trees * num_axial_trees + i_axial_trees) * 8 + 2] = i_tangential_trees * dphi;
            parameters[(i_tangential_trees * num_axial_trees + i_axial_trees) * 8 + 3] = 0.5 - i_axial_trees * dh;
            parameters[(i_tangential_trees * num_axial_trees + i_axial_trees) * 8 + 4] = (i_tangential_trees + 1) * dphi;
            parameters[(i_tangential_trees * num_axial_trees + i_axial_trees) * 8 + 5] = 0.5 + -(i_axial_trees + 1) * dh;
            parameters[(i_tangential_trees * num_axial_trees + i_axial_trees) * 8 + 6] = i_tangential_trees * dphi;
            parameters[(i_tangential_trees * num_axial_trees + i_axial_trees) * 8 + 7] = 0.5 -(i_axial_trees + 1) * dh;
          }

          int edges[24] = {-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
                           -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1};

          /* If geometry on both sides of cell */
          if (num_radial_trees == 1)
          {
            #if T8_WITH_OCC
            /* Assign occ geometries to the corresponding faces */
            int faces[6] = {-1, -1, -1, -1, -1, -1};
            faces[0] = cylinder_outer_index;
            faces[1] = cylinder_inner_index;
            

            /* Assign attributes to cmesh cells */
            t8_cmesh_set_attribute (cmesh, (i_tangential_trees * num_axial_trees + i_axial_trees) * num_radial_trees + i_radial_trees, t8_get_package_id(), 
                                    T8_CMESH_OCC_FACE_ATTRIBUTE_KEY, 
                                    faces, 6 * sizeof(int), 0);
            t8_cmesh_set_attribute (cmesh, (i_tangential_trees * num_axial_trees + i_axial_trees) * num_radial_trees + i_radial_trees, t8_get_package_id(), 
                                    T8_CMESH_OCC_EDGE_ATTRIBUTE_KEY, 
                                    edges, 24 * sizeof(int), 0);
            t8_cmesh_set_attribute (cmesh, (i_tangential_trees * num_axial_trees + i_axial_trees) * num_radial_trees + i_radial_trees, t8_get_package_id(), 
                                    T8_CMESH_OCC_FACE_PARAMETERS_ATTRIBUTE_KEY + 0, 
                                    parameters + (i_tangential_trees * num_axial_trees + i_axial_trees) * 8, 8 * sizeof(double), 0);
            t8_cmesh_set_attribute (cmesh, (i_tangential_trees * num_axial_trees + i_axial_trees) * num_radial_trees + i_radial_trees, t8_get_package_id(), 
                                    T8_CMESH_OCC_FACE_PARAMETERS_ATTRIBUTE_KEY + 1, 
                                    parameters + (i_tangential_trees * num_axial_trees + i_axial_trees) * 8, 8 * sizeof(double), 0);
            #endif /* T8_WITH_OCC */
          }
          /* If geometry only on face 1 */
          else if (i_radial_trees == 0)
          {
            #if T8_WITH_OCC
            /* Assign occ geometries to the corresponding faces */
            int faces[6] = {-1, -1, -1, -1, -1, -1};
            faces[1] = cylinder_inner_index;

            /* Assign attributes to cmesh cells */
            t8_cmesh_set_attribute (cmesh, (i_tangential_trees * num_axial_trees + i_axial_trees) * num_radial_trees + i_radial_trees, t8_get_package_id(), 
                                    T8_CMESH_OCC_FACE_ATTRIBUTE_KEY, 
                                    faces, 6 * sizeof(int), 0);
            t8_cmesh_set_attribute (cmesh, (i_tangential_trees * num_axial_trees + i_axial_trees) * num_radial_trees + i_radial_trees, t8_get_package_id(), 
                                    T8_CMESH_OCC_EDGE_ATTRIBUTE_KEY, 
                                    edges, 24 * sizeof(int), 0);
            t8_cmesh_set_attribute (cmesh, (i_tangential_trees * num_axial_trees + i_axial_trees) * num_radial_trees + i_radial_trees, t8_get_package_id(), 
                                    T8_CMESH_OCC_FACE_PARAMETERS_ATTRIBUTE_KEY + 1, 
                                    parameters + (i_tangential_trees * num_axial_trees + i_axial_trees) * 8, 8 * sizeof(double), 0);
            #endif /* T8_WITH_OCC */
          }
          /* If geometry only on face 0 */
          else if (i_radial_trees == num_radial_trees - 1)
          {
            #if T8_WITH_OCC
            /* Assign occ geometries to the corresponding faces */
            int faces[6] = {-1, -1, -1, -1, -1, -1};
            faces[0] = cylinder_outer_index;

            /* Assign attributes to cmesh cells */
            t8_cmesh_set_attribute (cmesh, (i_tangential_trees * num_axial_trees + i_axial_trees) * num_radial_trees + i_radial_trees, t8_get_package_id(), 
                                    T8_CMESH_OCC_FACE_ATTRIBUTE_KEY, 
                                    faces, 6 * sizeof(int), 0);
            t8_cmesh_set_attribute (cmesh, (i_tangential_trees * num_axial_trees + i_axial_trees) * num_radial_trees + i_radial_trees, t8_get_package_id(), 
                                    T8_CMESH_OCC_EDGE_ATTRIBUTE_KEY, 
                                    edges, 24 * sizeof(int), 0);
            t8_cmesh_set_attribute (cmesh, (i_tangential_trees * num_axial_trees + i_axial_trees) * num_radial_trees + i_radial_trees, t8_get_package_id(), 
                                    T8_CMESH_OCC_FACE_PARAMETERS_ATTRIBUTE_KEY + 0, 
                                    parameters + (i_tangential_trees * num_axial_trees + i_axial_trees) * 8, 8 * sizeof(double), 0);
            #endif /* T8_WITH_OCC */
          }
          /* If there is no geometry */
          else
          {
            #if T8_WITH_OCC
            /* Assign occ geometries to the corresponding faces */
            int faces[6] = {-1, -1, -1, -1, -1, -1};

            /* Assign attributes to cmesh cells */
            t8_cmesh_set_attribute (cmesh, (i_tangential_trees * num_axial_trees + i_axial_trees) * num_radial_trees + i_radial_trees, t8_get_package_id(), 
                                    T8_CMESH_OCC_FACE_ATTRIBUTE_KEY, 
                                    faces, 6 * sizeof(int), 0);
            t8_cmesh_set_attribute (cmesh, (i_tangential_trees * num_axial_trees + i_axial_trees) * num_radial_trees + i_radial_trees, t8_get_package_id(), 
                                    T8_CMESH_OCC_EDGE_ATTRIBUTE_KEY, 
                                    edges, 24 * sizeof(int), 0);
            #endif /* T8_WITH_OCC */
          }
        }
        /* Join radial neighbors */
        if (i_radial_trees > 0)
        {
          t8_cmesh_set_join (cmesh, 
                            (i_tangential_trees * num_axial_trees + i_axial_trees) * num_radial_trees + i_radial_trees - 1, 
                            (i_tangential_trees * num_axial_trees + i_axial_trees) * num_radial_trees + i_radial_trees, 
                            0, 1, 0);
        }
      }

      /* Join axial neighbors */
      if (i_axial_trees > 0)
      {
        for (int i_radial_trees = 0; i_radial_trees < num_radial_trees; ++i_radial_trees)
        {
          t8_cmesh_set_join (cmesh, 
                            (i_tangential_trees * num_axial_trees + i_axial_trees - 1) * num_radial_trees + i_radial_trees, 
                            (i_tangential_trees * num_axial_trees + i_axial_trees) * num_radial_trees + i_radial_trees, 
                            5, 4, 0);
        }
      }
    }
  }
    
  /* Join tangential neighbors of seam */
  for (int i_axial_trees = 0; i_axial_trees < num_axial_trees; ++i_axial_trees)
  {
    for (int i_radial_trees = 0; i_radial_trees < num_radial_trees; ++i_radial_trees)
    {
      t8_cmesh_set_join (cmesh, 
                        (num_axial_trees + i_axial_trees) * num_radial_trees + i_radial_trees, 
                        ((num_tangential_trees - 1) * num_axial_trees + i_axial_trees) * num_radial_trees + i_radial_trees, 
                        2, 3, 0);
    }
  }
  
  /* Join all other tangential neighbors */
  for (int i_tangential_trees = 1; i_tangential_trees < num_tangential_trees; ++i_tangential_trees)
  {
    for (int i_axial_trees = 0; i_axial_trees < num_axial_trees; ++i_axial_trees)
    {
      for (int i_radial_trees = 0; i_radial_trees < num_radial_trees; ++i_radial_trees)
      {
        t8_cmesh_set_join (cmesh, 
                          ((i_tangential_trees - 1) * num_axial_trees + i_axial_trees) * num_radial_trees + i_radial_trees, 
                          (i_tangential_trees * num_axial_trees + i_axial_trees) * num_radial_trees + i_radial_trees,
                          2, 3, 0);
      }
    }
  }
  
  if (with_occ_geometry)
  {
    t8_cmesh_register_geometry (cmesh, geometry_occ);
  }
  else
  {
    t8_cmesh_register_geometry (cmesh, geometry_linear);
  }
  t8_cmesh_commit (cmesh, comm);
  T8_FREE(vertices);
  T8_FREE(parameters);
  return cmesh;
}
