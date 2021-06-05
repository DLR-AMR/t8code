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

t8_cmesh_t
t8_cmesh_new_hollow_cylinder (sc_MPI_Comm comm, int num_tangential_trees, 
                              int num_axial_trees,  int with_occ_geometry, 
                              int do_partition)
{  
  t8_cmesh_t cmesh;
  t8_cmesh_init (&cmesh);
  t8_geometry_c *geometry;
  
  if (with_occ_geometry)
  {
    #if T8_WITH_OCC
    /* Create occ cylinder surfaces */
    double radius_inner = 0.25;
    double radius_outer = 0.5;
    gp_Pnt origin(0, 0, 0);
    gp_Dir z_dir(0, 0, 1);
    gp_Ax2 axis(origin, z_dir);
    gp_Vec height(0, 0, 1);
    gp_Circ circle_outer(axis, radius_outer);
    gp_Circ circle_inner(axis, radius_inner);
    BRepBuilderAPI_MakeEdge make_outer_edge(circle_outer);
    TopoDS_Edge edge_outer = make_outer_edge.Edge();
    TopoDS_Face face_outer = TopoDS::Face(BRepPrimAPI_MakePrism(edge_outer, height));
    Handle_Geom_Surface cylinder_outer = BRep_Tool::Surface(face_outer);
    BRepBuilderAPI_MakeEdge make_inner_edge(circle_inner);
    TopoDS_Edge edge_inner = make_inner_edge.Edge();
    TopoDS_Face face_inner = TopoDS::Face(BRepPrimAPI_MakePrism(edge_inner, height));
    Handle_Geom_Surface cylinder_inner = BRep_Tool::Surface(face_inner);
    t8_global_occ_surface[0] = cylinder_outer;
    t8_global_occ_surface[1] = cylinder_inner;
    geometry = new t8_geometry_occ (3, "occ surface dim=3", NULL);
    #else /* !T8_WITH_OCC */
    SC_ABORTF("OCC not linked");
    #endif /* T8_WITH_OCC */ 
  }
  else
  {
    geometry = t8_geometry_linear_new (3);
  }
  
  double *vertices, *parameters;
  double radius_outer = 0.5, radius_inner = 0.25;
  vertices = T8_ALLOC(double, num_tangential_trees * num_axial_trees * 24);
  parameters = T8_ALLOC(double, num_tangential_trees * num_axial_trees * 8);
  int faces[6] = {0, 1, -1, -1, -1, -1};
  for (int i = 0; i < num_tangential_trees; ++i)
  {
    for (int j = 0; j < num_axial_trees; ++j)
    {
      t8_cmesh_set_tree_class (cmesh, i * num_axial_trees + j, T8_ECLASS_HEX);
      vertices[(i * num_axial_trees + j) * 24 + 0] = cos((i + 1) * 2 * M_PI / num_tangential_trees) * radius_outer;
      vertices[(i * num_axial_trees + j) * 24 + 1] = sin((i + 1) * 2 * M_PI / num_tangential_trees) * radius_outer;
      vertices[(i * num_axial_trees + j) * 24 + 2] = j / num_axial_trees;
      vertices[(i * num_axial_trees + j) * 24 + 3] = cos((i + 1) * 2 * M_PI / num_tangential_trees) * radius_inner;
      vertices[(i * num_axial_trees + j) * 24 + 4] = sin((i + 1) * 2 * M_PI / num_tangential_trees) * radius_inner;
      vertices[(i * num_axial_trees + j) * 24 + 5] = j / num_axial_trees;
      vertices[(i * num_axial_trees + j) * 24 + 6] = cos(i * 2 * M_PI / num_tangential_trees) * radius_outer;
      vertices[(i * num_axial_trees + j) * 24 + 7] = sin((i) * 2 * M_PI / num_tangential_trees) * radius_outer;
      vertices[(i * num_axial_trees + j) * 24 + 8] = j / num_axial_trees;
      vertices[(i * num_axial_trees + j) * 24 + 9] = cos(i * 2 * M_PI / num_tangential_trees) * radius_inner;
      vertices[(i * num_axial_trees + j) * 24 + 10] = sin(i * 2 * M_PI / num_tangential_trees) * radius_inner;
      vertices[(i * num_axial_trees + j) * 24 + 11] = j / num_axial_trees;
      vertices[(i * num_axial_trees + j) * 24 + 12] = cos((i + 1) * 2 * M_PI / num_tangential_trees) * radius_outer;
      vertices[(i * num_axial_trees + j) * 24 + 13] = sin((i + 1) * 2 * M_PI / num_tangential_trees) * radius_outer;
      vertices[(i * num_axial_trees + j) * 24 + 14] = (j + 1) / num_axial_trees;
      vertices[(i * num_axial_trees + j) * 24 + 15] = cos((i + 1) * 2 * M_PI / num_tangential_trees) * radius_inner;
      vertices[(i * num_axial_trees + j) * 24 + 16] = sin((i + 1) * 2 * M_PI / num_tangential_trees) * radius_inner;
      vertices[(i * num_axial_trees + j) * 24 + 17] = (j + 1) / num_axial_trees;
      vertices[(i * num_axial_trees + j) * 24 + 18] = cos(i * 2 * M_PI / num_tangential_trees) * radius_outer;
      vertices[(i * num_axial_trees + j) * 24 + 19] = sin((i) * 2 * M_PI / num_tangential_trees) * radius_outer;
      vertices[(i * num_axial_trees + j) * 24 + 20] = (j + 1) / num_axial_trees;
      vertices[(i * num_axial_trees + j) * 24 + 21] = cos(i * 2 * M_PI / num_tangential_trees) * radius_inner;
      vertices[(i * num_axial_trees + j) * 24 + 22] = sin(i * 2 * M_PI / num_tangential_trees) * radius_inner;
      vertices[(i * num_axial_trees + j) * 24 + 23] = (j + 1) / num_axial_trees;
      t8_cmesh_set_tree_vertices (cmesh, i * num_axial_trees + j, vertices + (i * num_axial_trees + j) * 24, 24);
      
      if (with_occ_geometry)
      {
        #if T8_WITH_OCC
        parameters[(i * num_axial_trees + j) * 8 + 0] = (i + 1) * 2 * M_PI / num_tangential_trees;
        parameters[(i * num_axial_trees + j) * 8 + 1] = j / num_axial_trees;
        parameters[(i * num_axial_trees + j) * 8 + 2] = i * 2 * M_PI / num_tangential_trees;
        parameters[(i * num_axial_trees + j) * 8 + 3] = j / num_axial_trees;
        parameters[(i * num_axial_trees + j) * 8 + 4] = (i + 1) * 2 * M_PI / num_tangential_trees;
        parameters[(i * num_axial_trees + j) * 8 + 5] = -(j + 1) / num_axial_trees;
        parameters[(i * num_axial_trees + j) * 8 + 6] = i * 2 * M_PI / num_tangential_trees;
        parameters[(i * num_axial_trees + j) * 8 + 7] = -(j + 1) / num_axial_trees;
        t8_cmesh_set_attribute (cmesh, i * num_axial_trees + j, t8_get_package_id(), T8_CMESH_OCC_SURFACE_ATTRIBUTE_KEY, 
                                faces, 6 * sizeof(int), 1);
        t8_cmesh_set_attribute (cmesh, i * num_axial_trees + j, t8_get_package_id(), T8_CMESH_OCC_SURFACE_PARAMETERS_ATTRIBUTE_KEY + 0, 
                                parameters + (i * num_axial_trees + j) * 8, 8 * sizeof(double), 1);
        t8_cmesh_set_attribute (cmesh, i * num_axial_trees + j, t8_get_package_id(), T8_CMESH_OCC_SURFACE_PARAMETERS_ATTRIBUTE_KEY + 1, 
                                parameters + (i * num_axial_trees + j) * 8, 8 * sizeof(double), 1);
        #endif /* T8_WITH_OCC */
      }
    }
    /* Join axial neighbours */
    for (int j = 0; j < num_axial_trees - 1; ++j)
    {
      t8_cmesh_set_join(cmesh, i * num_axial_trees + j, i * num_axial_trees + j + 1, 5, 4, 0);
    }
  }
    
  /* Join tangential neighbours of seam */
  for (int j = 0; j < num_axial_trees; ++j)
    {
      t8_cmesh_set_join(cmesh, num_tangential_trees - 1, 0, 2, 3, 0);
    }
  
  /* Join all other tangential neighbours */
  for (int i = 0; i < num_tangential_trees - 1; ++i)
  {
    for (int j = 0; j < num_axial_trees; ++j)
    {
      t8_cmesh_set_join(cmesh, i * num_axial_trees + j, (i + 1) * num_axial_trees + j, 2, 3, 0);
    }
  }
  
  t8_cmesh_register_geometry (cmesh, geometry);
 
  if (do_partition) 
  {
    int                 mpirank, mpisize, mpiret;
    int                 first_tree, last_tree;
    mpiret = sc_MPI_Comm_rank (comm, &mpirank);
    SC_CHECK_MPI (mpiret);
    mpiret = sc_MPI_Comm_size (comm, &mpisize);
    SC_CHECK_MPI (mpiret);
    first_tree = (mpirank * num_axial_trees * num_tangential_trees) / mpisize;
    last_tree = ((mpirank + 1) * num_axial_trees * num_tangential_trees) / mpisize - 1;
    t8_cmesh_set_partition_range (cmesh, 3, first_tree, last_tree);
  }
  t8_cmesh_commit (cmesh, comm);
  T8_FREE(vertices);
  if (with_occ_geometry) {
    #if T8_WITH_OCC
    T8_FREE(parameters);
    #endif /* T8_WITH_OCC */
  }
  return cmesh;
}