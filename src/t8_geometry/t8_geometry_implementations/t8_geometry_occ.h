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

/** \file t8_geometry_occ.h
 * This header provides the C interface to create a occ geometry.
 */

#ifndef T8_GEOMETRY_OCC_H
#define T8_GEOMETRY_OCC_H

#include <t8.h>
#include <t8_geometry/t8_geometry.h>

/* The vertices of each edge of a hexahedron. Used in the occ geometry. */
extern const int    t8_edge_vertex_to_tree_vertex[T8_ECLASS_MAX_EDGES][2];

/* The faces connected to each edge. */
extern const int    t8_edge_to_face[T8_ECLASS_MAX_EDGES][2];

/* The edges connected to a face */
extern const int    t8_face_to_edge[T8_ECLASS_MAX_FACES][T8_ECLASS_MAX_EDGES];

/** This typedef holds virtual functions for a particular geometry.
 * We need it so that we can use t8_geometry_occ_c pointers in .c files
 * without them seeing the actual C++ code (and then not compiling)
 */
typedef struct t8_geometry_occ t8_geometry_occ_c;

T8_EXTERN_C_BEGIN ();

/** Create a new occ geometry of a given dimension.
 * \param [in] dimension  0 <= \a dimension <= 3. The dimension.
 * \param [in] fileprefix Prefix of a .brep file from which to extract an occ geometry.
 * \param [in] name       The name to give this geometry.
 * \return                A pointer to an allocated t8_geometry_occ struct, as
 *                        if the t8_geometry_occ (int dimension, const char *fileprefix, const char *name) 
 *                        constructor was called.
 */
t8_geometry_occ_c  *t8_geometry_occ_new (int dimension,
                                         const char *fileprefix,
                                         const char *name_in);

/** Destroy a occ geometry that was created with \ref t8_geometry_occ_new.
 * \param [in,out] geom A occ geometry. Set to NULL on output.
 */
void                t8_geometry_occ_destroy (t8_geometry_occ_c ** geom);

T8_EXTERN_C_END ();

#endif /* !T8_GEOMETRY_OCC_H! */
