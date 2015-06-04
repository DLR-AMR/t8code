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

/** \file t8_geometry.h
 * We define the geometry transformation for a forest of trees in this file.
 * This is mainly used for vtk output.
 */

#ifndef T8_GEOMETRY_H
#define T8_GEOMETRY_H

#include <t8.h>

T8_EXTERN_C_BEGIN ();

typedef struct t8_geometry *t8_geometry_t;

/** Forward transformation from the reference unit square to physical space.
 * Note that the two-dimensional connectivities have 3D vertex coordinates
 * that can be used in the transformation if so desired.
 * The physical space "xyz" is user-defined, currently used for VTK output.
 * \param [in] geom The underlying t8_geometry_t struct.
 * \param [in] which_tree The tree_id of the coarse tree to be considered.
 * \param [in] abc  The reference coordinates in [0,1]^d to be transformed.
 * \param [out] xyz The physical coordinates that abc get mapped to.
 */
typedef void        (*t8_geometry_X_t) (t8_geometry_t geom,
                                        t8_topidx_t which_tree,
                                        const double abc[3], double xyz[3]);

/** Destructor prototype for a user-allocated \a t8_geometry_t.
 * It is invoked by t8_geometry_reset.  If the user chooses to
 * reserve the structure statically, simply don't call t8_geometry_reset.
 * \param [in] geom The geometry to be reseted.
 */
typedef void        (*t8_geometry_reset_t) (t8_geometry_t * pgeom);

/* TODO: Comment these out */

void                t8_geometry_init (t8_geometry_t * pgeom);

void                t8_geometry_set_name (t8_geometry_t geom,
                                          const char *name);

void                t8_geometry_set_user (t8_geometry_t geom, void *user);

void                t8_geometry_set_transformation (t8_geometry_t geom,
                                                    t8_geometry_X_t X);

void                t8_geometry_set_reset (t8_geometry_t geom,
                                           t8_geometry_reset_t reset);

void                t8_geometry_ref (t8_geometry_t geom);

void                t8_geometry_unref (t8_geometry_t * pgeom);

void                t8_geometry_reset (t8_geometry_t * pgeom);

/** Create a geometry that maps the unit square to itself via the identity mapping.
 * This function exists to provide the minimal example of a t8_geometry_t.
 * It should not be used for coarse meshes with more than one tree.
 */
t8_geometry_t       t8_geometry_new_identity (void);

T8_EXTERN_C_END ();

#endif /* !T8_GEOMETRY_H! */
