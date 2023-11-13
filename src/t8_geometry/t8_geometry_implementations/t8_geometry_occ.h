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
#include <t8_geometry/t8_geometry_with_vertices.h>

/** This typedef holds virtual functions for a particular geometry.
 * We need it so that we can use t8_geometry_occ_c pointers in .c files
 * without them seeing the actual C++ code (and then not compiling)
 */
typedef struct t8_geometry_occ t8_geometry_occ_c;

T8_EXTERN_C_BEGIN ();

/**
 * Create a new occ geometry with a given dimension. The geometry
 * is currently viable with quad/hex and triangle trees. Tets will be supported soon.
 * The geometry uses as many vertices as the tree type has, as well as
 * additional geometry information, which is extracted from a .brep file.
 * The vertices are saved via the \ref t8_cmesh_set_tree_vertices function.
 * Since the internals of this geometry are finely tuned to the .brep file
 * it is recommended to only use it with the \ref t8_cmesh_readmshfile function.
 * \param [in] dim        0 <= \a dimension <= 3. The dimension.
 * \param [in] fileprefix Prefix of a .brep file from which to extract an occ geometry.
 * \param [in] name       The name to give this geometry.
 * \return                A pointer to an allocated t8_geometry_occ struct, as
 *                        if the \ref t8_geometry_occ (int dim, const *char fileprefix, 
 *                        const char *name) 
 *                        constructor was called.
 */
t8_geometry_occ_c *
t8_geometry_occ_new (int dim, const char *fileprefix, const char *name_in);

/** Destroy a occ geometry that was created with \ref t8_geometry_occ_new.
 * \param [in,out] geom A occ geometry. Set to NULL on output.
 */
void
t8_geometry_occ_destroy (t8_geometry_occ_c **geom);

#if T8_ENABLE_DEBUG
/** Query whether a given geometry is \ref t8_geometry_occ_c.
 * \param [in] geometry   A geometry.
 * \return     True (non-zero) if and only if the geometry is of type \ref t8_geometry_occ_c.
 * \note       This function is currently only available in debug mode.
 */
int
t8_geom_is_occ (const t8_geometry_c *geometry);
#endif

T8_EXTERN_C_END ();

#endif /* !T8_GEOMETRY_OCC_H! */
