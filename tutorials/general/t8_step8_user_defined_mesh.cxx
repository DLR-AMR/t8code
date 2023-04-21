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



T8_EXTERN_C_END ();
