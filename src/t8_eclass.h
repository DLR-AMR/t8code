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

/** \file t8_eclass.h
 * We define all possible element classes that occur in hybrid meshes.
 *
 * Notable examples are triangles, tetrahedra, quadrilaterals and hexahedra.
 * We cover all dimensions between zero and three, so it is in principal
 * possible to build a topological complex out of these element classes.
 */

#ifndef T8_ECLASS_H
#define T8_ECLASS_H

#include <t8.h>

T8_EXTERN_C_BEGIN ();

/** This enumeration contains all possible element classes. */
typedef enum t8_eclass
{
  T8_ECLASS_ZERO = 0,
  /** The vertex is the only zero-dimensional element class. */
  T8_ECLASS_VERTEX = T8_ECLASS_ZERO,
  /** The line is the only one-dimensional element class. */
  T8_ECLASS_LINE,
  /** The quadrilateral is one of two element classes in two dimensions. */
  T8_ECLASS_QUAD,
  /** The element class for a triangle. */
  T8_ECLASS_TRIANGLE,
  /** The hexahedron is one three-dimensional element class. */
  T8_ECLASS_HEX,
  /** The tetrahedron is another three-dimensional element class. */
  T8_ECLASS_TET,
  /** The prism has five sides: two opposing triangles joined by three quadrilaterals. */
  T8_ECLASS_PRISM,
  /** The pyramid has a quadrilateral as base and four triangles as sides. */
  T8_ECLASS_PYRAMID,
  /** This is no element class but can be used as the number of element classes. */
  T8_ECLASS_COUNT
}
t8_eclass_t;

/** The MPI datatype used for t8_eclass_t */
#define T8_MPI_ECLASS_TYPE (T8_ASSERT (sizeof (int) == sizeof (t8_eclass_t)),\
  sc_MPI_INT)

/** The maximum number of boundary faces an element class can have. */
#define T8_ECLASS_MAX_FACES 6
/** The maximum number of cornes a 2-dimensional element class can have. */
#define T8_ECLASS_MAX_CORNERS_2D 4
/** The maximum number of cornes an element class can have. */
#define T8_ECLASS_MAX_CORNERS 8
/** The maximal possible dimension for an eclass */
#define T8_ECLASS_MAX_DIM 3

/** Map each of the element classes to its dimension. */
extern const int    t8_eclass_to_dimension[T8_ECLASS_COUNT];

/** The number of codimension-one boundaries of an element class. */
extern const int    t8_eclass_num_faces[T8_ECLASS_COUNT];

/** For each dimension the maximum possible number of faces of an eclass of that dimension. */
extern const int    t8_eclass_max_num_faces[T8_ECLASS_MAX_DIM + 1];

/** For each eclass and each face f the entry i gives the vertex number
 * of f's i-th vertex within all vertices of the tree. */
extern const int
     t8_face_vertex_to_tree_vertex[T8_ECLASS_COUNT][T8_ECLASS_MAX_FACES]
  [T8_ECLASS_MAX_CORNERS_2D];

/** Each face is either 0 or 1 oriented, depending on the order of its vertices.
 * We say a face is 0 oriented, if its normal vector points inwards,
 * 1 oriented otherwise.
 * The normal vector is computed as the cross product of v_1 - v_0 and v_2 - v_0.
 * v_i being the i-th vertex.
 * The faces of an eclass of dimension 2 or lower are all 0 oriented.
 */
extern const int
     t8_eclass_face_orientation[T8_ECLASS_COUNT][T8_ECLASS_MAX_FACES];

/** The number of vertices of an element class. */
extern const int    t8_eclass_num_vertices[T8_ECLASS_COUNT];

/** The vtk cell type for the eclass */
extern const int    t8_eclass_vtk_type[T8_ECLASS_COUNT];

/** Map the t8code corner number to the vtk corner number */
extern const int
     t8_eclass_vtk_corner_number[T8_ECLASS_COUNT][T8_ECLASS_MAX_CORNERS];

/** The shape of an element according its eclass*/
extern const t8_eclass_t
    t8_eclass_shape[T8_ECLASS_COUNT];

/** For each of the element classes, list the type of the faces. */
extern const int
     t8_eclass_face_types[T8_ECLASS_COUNT][T8_ECLASS_MAX_FACES];

/** For each of the element classes, count the boundary points. */
extern const int
     t8_eclass_boundary_count[T8_ECLASS_COUNT][T8_ECLASS_COUNT];

/** For each eclass, the name of this class as a string */
extern const char  *t8_eclass_to_string[T8_ECLASS_COUNT];

/** Query the element class and count of boundary points.
 * \param [in] theclass         We query a point of this element class.
 * \param [in] min_dim          Ignore boundary points of lesser dimension.
 *                              The ignored points get a count value of 0.
 * \param [out] per_eclass      Array of length T8_ECLASS_COUNT to be filled
 *                              with the count of the boundary objects,
 *                              counted per each of the element classes.
 * \return                      The count over all boundary points.
 */
int                 t8_eclass_count_boundary (t8_eclass_t theclass,
                                              int min_dim, int *per_eclass);

/** Compare two eclasses of the same dimension
 *  as necessary for face neighbor orientation.
 *  The implemented order is Triangle < Square in 2D and
 *  Tet < Hex < Prism < Pyramid in 3D.
 *  \param [in] eclass1 The first eclass to compare.
 *  \param [in] eclass2 The second eclass to compare.
 *  \return 0 if the eclasses are equal, 1 if eclass1 > eclass2
 *            and -1 if eclass1 < eclass2
 */
int                 t8_eclass_compare (t8_eclass_t eclass1,
                                       t8_eclass_t eclass2);

T8_EXTERN_C_END ();

#endif /* !T8_ELEMENT_H */
