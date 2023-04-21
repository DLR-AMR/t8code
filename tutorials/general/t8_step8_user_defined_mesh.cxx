/*
  This file is part of t8code.
  t8code is a C library to manage a collection (a forest) of multiple
  connected adaptive space-trees of general element types in parallel.

  Copyright (C) 2023 the developers

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

/* See also: https://github.com/DLR-AMR/t8code/wiki/Step-6-Computing-stencils
 *
 * This is step8 of the t8code tutorials using the C++ interface of t8code.
 * In the following we will create two user defined meshes.
 * The first example is given by a periodic two dimensional mesh using linear
 * geometry consisting of four triangles and and two quads.
 * The second example is given by a non-periodic three dimensional mesh 
 * with linear geometry constructed using one tetrahedron, two prisms and one 
 * hexaedron.
 *
 * How you can experiment here:
 *   - Look at the paraview output files of the different meshes.
 *   - Change the element types of the mesh.
 *   - Change the face connections between the different elements.
 *   - Create an own mesh.
 *  */

#include <t8.h>                 /* General t8code header, always include this. */
#include <t8_cmesh.h>           /* cmesh definition and basic interface. */
#include <t8_forest.h>          /* forest definition and basic interface. */
#include <t8_schemes/t8_default/t8_default_cxx.hxx>     /* default refinement scheme. */
#include <t8_cmesh_vtk_writer.h> /* write file in vtu file */
#include <t8_geometry/t8_geometry_implementations/t8_geometry_linear.h> /* linear geometry of the cmesh */

T8_EXTERN_C_BEGIN ();

/* Different steps of creating a cmesh
 * 1. Defining an array with all vertices
 *    The vertices are ordered in a listing for each cell.
 *    Thus, there can be duplicates in the List.
 *    Example: double vertices[numberOfValues] = {
 *               //point values for cell 1
 *               x_1,y_1,z_1         //(x,y,z) of first point of cell 1
 *               x_2,y_2,z_2         //(x,y,z) of second point of cell 1
 *                   .
 *                   .
 *                   .
 *               x_n,y_n,z_n         //(x,y,z) of nth point (last point) of cell 1
 *
 *               //point values for cell 2
 *               x_1,y_1,z_1         //(x,y,z) of first point of cell 2
 *               x_2,y_2,z_2         //(x,y,z) of second point of cell 2
 *                   .
 *                   .
 *                   .
 *               x_m,y_m,z_m         //(x,y,z) of nth point (last point) of cell 2
 *
 *                   .
 *                   .
 *                   .
 *
 *               //point values for the last cell 
 *               x_1,y_1,z_1         //(x,y,z) of first point of the last cell
 *               x_2,y_2,z_2         //(x,y,z) of second point of the last cell
 *                   .
 *                   .
 *                   .
 *               x_o,y_o,z_o         //(x,y,z) of nth point (last point) of the last cell
 *             };
 *
 * 2. Initialization of the mesh
 *    Example: t8_cmesh_t          cmesh;
 *             t8_cmesh_init (&cmesh);
 * 
 *
 * 3. Definition of the geometry
 *             t8_geometry_c      *linear_geom = [defineTheGeometry];
 *
 * 4. Defitition of the classes of the different trees - each tree is defined by one cell
 *    Example: //Class of the first tree
 *             t8_cmesh_set_tree_class (cmesh, 0, T8_ECLASS_[TYPE]);
 *             //Class of the second tree
 *             t8_cmesh_set_tree_class (cmesh, 1, T8_ECLASS_[TYPE]);
 *                   .
 *                   .
 *                   .
 *             //Class of the last tree
 *             t8_cmesh_set_tree_class (cmesh, x, T8_ECLASS_[TYPE]);
 *
 * 5. Classification of the vertices for each tree
 *             // Vertices of the first tree
 *             t8_cmesh_set_tree_vertices (cmesh, 0, [pointerToVerticesOfTreeOne], [numberOfVerticesTreeOne]);
 *             // Vertices of the second tree
 *             t8_cmesh_set_tree_vertices (cmesh, 1, [pointerToVerticesOfTreeTwo] , [numberOfVerticesTreeTwo]);
 *                   .
 *                   .
 *                   .
 *             // Vertices of the last tree
 *             t8_cmesh_set_tree_vertices (cmesh, x, [pointerToVerticesOfTree(x+1)] , [numberOfVerticesTree(x+1)]);
 *
 * 6. Definition of the face neighboors between the different trees
 *             // List of all face neighboor connections
 *             t8_cmesh_set_join (cmesh, [treeId1], [treeId2], [faceIdInTree1], [faceIdInTree2], [orientation]);
 *             t8_cmesh_set_join (cmesh, [treeId1], [treeId2], [faceIdInTree1], [faceIdInTree2], [orientation]);
 *                   .
 *                   .
 *                   .
 *             t8_cmesh_set_join (cmesh, [treeId1], [treeId2], [faceIdInTree1], [faceIdInTree2], [orientation]);
 *
 * 7. Commit the mesh
 *             t8_cmesh_commit (cmesh, comm);
 *  */



T8_EXTERN_C_END ();
