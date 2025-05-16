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
 * Typedef for the t8_geometry class in order to be usable as a pointer
 * from .c files.
 */

#ifndef T8_GEOMETRY_H
#define T8_GEOMETRY_H

#include <t8.h>
#include <t8_refcount.h>

/** This enumeration contains all possible geometries. */
typedef enum t8_geometry_type {
  /** The zero geometry maps all points to zero. */
  T8_GEOMETRY_TYPE_ZERO = 0,
  /** The linear geometry uses linear interpolations to interpolate between the tree vertices. */
  T8_GEOMETRY_TYPE_LINEAR,
  /** The linear, axis aligned geometry uses only 2 vertices, since it is axis aligned. */
  T8_GEOMETRY_TYPE_LINEAR_AXIS_ALIGNED,
  /** The Lagrange geometry uses a mapping with Lagrange polynomials to approximate curved elements . */
  T8_GEOMETRY_TYPE_LAGRANGE,
  /** The analytic geometry uses a user-defined analytic function to map into the physical domain. */
  T8_GEOMETRY_TYPE_ANALYTIC,
  /** The opencascade geometry uses CAD shapes to map trees exactly to the underlying CAD model. */
  T8_GEOMETRY_TYPE_CAD,
  /** This is no geometry type but can be used as the number of geometry types. */
  T8_GEOMETRY_TYPE_COUNT,
  /** This is no geometry type but is used as error type to describe invalid geometries */
  T8_GEOMETRY_TYPE_INVALID
    /** This is no geometry type but is used for every geometry, where no type is defined */
    T8_GEOMETRY_TYPE_UNDEFINED
} t8_geometry_type_t;

/** This typedef holds virtual functions for a particular geometry.
 * We need it so that we can use t8_geometry_c pointers in .c files
 * without them seeing the actual C++ code (and then not compiling)
 */
typedef struct t8_geometry t8_geometry_c;

/** This typedef holds virtual functions for the geometry handler.
 * We need it so that we can use t8_geometry_handler_c pointers in .c files
 * without them seeing the actual C++ code (and then not compiling)
 * TODO: Delete this when the cmesh is a proper cpp class.
 */
typedef struct t8_geometry_handler t8_geometry_handler_c;

/* The t8_geometry_c type must be know to cmesh.h, thus we
 * include it after the typedef. */
#include <t8_cmesh.h>

T8_EXTERN_C_BEGIN ();

/**
 * Evaluates the geometry of a tree at a given reference point.
 * \param [in]  cmesh      The cmesh
 * \param [in]  gtreeid    The global id of the tree
 * \param [in]  ref_coords The reference coordinates at which to evaluate the geometry
 * \param [in]  num_coords The number of reference coordinates
 * \param [out] out_coords The evaluated coordinates
 */
void
t8_geometry_evaluate (t8_cmesh_t cmesh, t8_gloidx_t gtreeid, const double *ref_coords, const size_t num_coords,
                      double *out_coords);

/** Evaluates the jacobian of a tree at a given reference point.
 * \param[in]  cmesh      The cmesh
 * \param[in]  gtreeid    The global id of the tree
 * \param[in]  ref_coords The reference coordinates at which to evaluate the jacobian
 * \param[in]  num_coords The number of reference coordinates
 * \param[out] jacobian   The jacobian at the reference coordinates
 */
void
t8_geometry_jacobian (t8_cmesh_t cmesh, t8_gloidx_t gtreeid, const double *ref_coords, const size_t num_coords,
                      double *jacobian);

/** This function returns the geometry type of a tree.
 * \param[in] cmesh       The cmesh
 * \param[in] gtreeid     The global id of the tree
 * \return                The geometry type of the tree with id \ref gtreeid
 */
t8_geometry_type_t
t8_geometry_get_type (t8_cmesh_t cmesh, t8_gloidx_t gtreeid);

/**
 * Check if a tree has a negative volume
 * 
 * \param[in] cmesh       The cmesh to check
 * \param[in] gtreeid     The global id of the tree
 * \return                True if the tree with id \ref gtreeid has a negative volume. False otherwise.  
 */
int
t8_geometry_tree_negative_volume (const t8_cmesh_t cmesh, const t8_gloidx_t gtreeid);

T8_EXTERN_C_END ();

#endif /* !T8_GEOMETRY_H */
