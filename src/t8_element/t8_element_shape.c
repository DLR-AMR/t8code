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

#include <t8_element/t8_element_shape.h>

const int t8_element_shape_max_num_corner[T8_ECLASS_MAX_DIM + 1] = { 0, 2, 4, 8 };

/** The number of codimension-one boundaries of an element class. */
int
t8_element_shape_num_faces (int element_shape)
{
  return t8_eclass_num_faces[element_shape];
}

/** For each dimension the maximum possible number of faces of an element_shape of that dimension. */
int
t8_element_shape_max_num_faces (int element_shape)
{
  return t8_eclass_max_num_faces[element_shape];
}

/** The number of vertices of an element class. */
int
t8_element_shape_num_vertices (int element_shape)
{
  return t8_eclass_num_vertices[element_shape];
}

/** The vtk cell type for the element_shape */
int
t8_element_shape_vtk_type (int element_shape)
{
  return t8_eclass_vtk_type[element_shape];
}

/** Maps the t8code corner number of the element to the vtk corner number 
 * \param [in] element_shape  The shape of the element.
 * \param [in] index          The index of the corner in z-order (t8code numeration).
 * \return                    The corresponding vtk index. 
*/
int
t8_element_shape_t8_to_vtk_corner_number (int element_shape, int index)
{
  return t8_eclass_t8_to_vtk_corner_number[element_shape][index];
}

/** Maps the vtk corner number of the element to the t8code corner number
 * \param [in] element_shape  The shape of the element.
 * \param [in] index          The index of the corner in vtk ordering. 
 * \return                    The corresponding t8code index.
 */
int
t8_element_shape_vtk_to_t8_corner_number (int element_shape, int index)
{
  return t8_eclass_vtk_to_t8_corner_number[element_shape][index];
}

/** For each element_shape, the name of this class as a string */
const char*
t8_element_shape_to_string (int element_shape)
{
  return t8_eclass_to_string[element_shape];
}

/* Currently t8_element_shape equals t8_eclass, if they will differ, this function has to be adapted. */
int
t8_element_shape_compare (t8_element_shape_t element_shape1, t8_element_shape_t element_shape2)
{
  return t8_eclass_compare ((t8_eclass_t) element_shape1, (t8_eclass_t) element_shape2);
}
