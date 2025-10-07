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

/** \file t8_element_shape.h
 * TODO: comment
 */

#ifndef T8_ELEMENT_SHAPE_H
#define T8_ELEMENT_SHAPE_H

#include <t8.h>
#include <t8_eclass.h>

T8_EXTERN_C_BEGIN ();

/** Type definition for the geometric shape of an element.
 * Currently the possible shapes are the same as the possible element classes.
 * I.e. T8_ECLASS_VERTEX, T8_ECLASS_TET, etc... */
typedef t8_eclass_t t8_element_shape_t;

/** The MPI datatype used for t8_element_shape_t */
#define T8_MPI_ELEMENT_SHAPE_TYPE (T8_ASSERT (sizeof (int) == sizeof (t8_element_shape_t)), sc_MPI_INT)

/** The maximum number of boundary faces an element class can have. */
#define T8_ELEMENT_SHAPE_MAX_FACES 6
/** The maximum number of corners a 3-dimensional element class can have. */
#define T8_ELEMENT_SHAPE_MAX_CORNERS 8

/** Maximum possible number of corner nodes of an element in a specific dimension */
extern const int t8_element_shape_max_num_corner[T8_ECLASS_MAX_DIM + 1];

/** The number of codimension-one boundaries of an element class. */
int
t8_element_shape_num_faces (int element_shape);

/** For each dimension the maximum possible number of faces of an element_shape of that dimension. */
int
t8_element_shape_max_num_faces (int element_shape);

/** The number of vertices of an element class. */
int
t8_element_shape_num_vertices (int element_shape);

/** The vtk cell type for the element_shape */
int
t8_element_shape_vtk_type (int element_shape);

/** Maps the t8code corner number of the element to the vtk corner number 
 * \param [in] element_shape  The shape of the element.
 * \param [in] index          The index of the corner in z-order (t8code numeration).
 * \return                    The corresponding vtk index. 
*/
int
t8_element_shape_t8_to_vtk_corner_number (int element_shape, int index);

/** Maps the vtk corner number of the element to the t8code corner number
 * \param [in] element_shape  The shape of the element.
 * \param [in] index          The index of the corner in vtk ordering. 
 * \return                    The corresponding t8code index.
 */
int
t8_element_shape_t8_corner_number (int element_shape, int index);

/** For each element_shape, the name of this class as a string */
const char*
t8_element_shape_to_string (int element_shape);

/** Compare two element_shapes of the same dimension
 *  as necessary for face neighbor orientation.
 *  The implemented order is Triangle < Square in 2D and
 *  Tet < Hex < Prism < Pyramid in 3D.
 *  \param [in] element_shape1 The first element_shape to compare.
 *  \param [in] element_shape2 The second element_shape to compare.
 *  \return 0 if the element_shapes are equal, 1 if element_shape1 > element_shape2
 *            and -1 if element_shape1 < element_shape2
 */
int
t8_element_shape_compare (t8_element_shape_t element_shape1, t8_element_shape_t element_shape2);

T8_EXTERN_C_END ();

#endif /* !T8_ELEMENT_SHAPE_H */
