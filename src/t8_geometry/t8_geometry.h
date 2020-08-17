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

/** This typedef holds virtual functions for a particular geometry.
 * We need it so that we can use t8_geometry_c pointers in .c files
 * without them seeing the actual C++ code (and then not compiling)
 */
typedef struct t8_geometry t8_geometry_c;

/* The t8_geometry_c type must be know to cmesh.h, thus we
 * include it after the typedef. */
#include <t8_cmesh.h>

typedef struct t8_geometry_handler
{
  sc_array_t          registered_geometries;
                                        /**< Stores all geometries that are handled by this geometry_handler. */
  t8_geometry_c      *active_geometry;
                                  /**< Points to the currently loaded geometry (the geometry that was used last and is likely to be used next). */
  t8_gloidx_t         active_tree;
                              /**< The global tree id of the last tree for which geometry was used. */
  int                 is_committed;
                               /**< If true, no new geometries can be registered. */
  t8_refcount_t       rc;
                     /**< The reference count of the geometry handler. */
} t8_geometry_handler_t;

T8_EXTERN_C_BEGIN ();

/** 
 * Initialize a geometry handler. This will allocate memory.
 * \param [in,out]  pgeom_handler On input a pointer to unallocated memory.
 *                                On output this will point to an initialized geometry handler.
 */
void                t8_geom_handler_init (t8_geometry_handler_t **
                                          pgeom_handler);

/** 
 * Increase the reference counter of a geometry handler.
 * \param [in] geom_handler An initialize geometry handler.
 */
void                t8_geom_handler_ref (t8_geometry_handler_t *
                                         geom_handler);

/**
 * Decrease the reference count of a geometry handler, if the refcount
 * reaches 0, the handler will get destroyed.
 * \param [in,out] pgeom_handler Pointer to an initialized geometry handler.
 *                               If the refcount reaches 0, will point to NULL afterwards.
 */
void                t8_geom_handler_unref (t8_geometry_handler_t **
                                           pgeom_handler);

/**
 * Destroy a geometry handler and free its memory. 
 * This is only valid if its reference count is 1.
 * If you are unsure, call \ref t8_geom_handler_unref instead.
 * \param [in,out] pgeom_handler Pointer to an initialized geometry handler.
 *                               Will point to NULL afterwards.
 */
void                t8_geom_handler_destroy (t8_geometry_handler_t **
                                             pgeom_handler);

/**
 * Add a geometry to the geometry handler.
 * \param [in,out] geom_handler An initialized but not committed geometry handler.
 * \param [in]     geometry     The geometry to add.
 */
void                t8_geom_handler_register_geometry (t8_geometry_handler_t *
                                                       geom_handler,
                                                       const t8_geometry_c *
                                                       geometry);

/**
 * Commit a geometry handler. This specifies that no geometries will
 * be added to it and makes it ready to be used.
 * \param [in,out] geom_handler An initialized but not committed geometry handler.
 */
void                t8_geom_handler_commit (t8_geometry_handler_t *
                                            geom_handler);

/* Check if a geometry handler was committed. */
int                 t8_geom_handler_is_committed (const t8_geometry_handler_t
                                                  * geom_handler);

/**
 * Given a geometries name find that geometry in the geometry handler
 * and return it.
 * \param [in] geom_handler A committed geometry handler.
 * \param [in] name         The name of a geometry.
 * \return                  A pointer to the geomery or NULL if it was not found.
 */
const t8_geometry_c *t8_geom_handler_find_geometry (const
                                                    t8_geometry_handler_t *
                                                    geom_handler,
                                                    const char *name);

void                t8_geometry_evaluate (t8_cmesh_t cmesh,
                                          t8_gloidx_t gtreeid,
                                          const double *ref_coords,
                                          double *out_coords);

void                t8_geometry_jacobian (t8_cmesh_t cmesh,
                                          t8_gloidx_t gtreeid,
                                          const double *ref_coords,
                                          double *jacobian);

T8_EXTERN_C_END ();

#endif /* !T8_GEOMETRY_H! */
