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

/** \file t8_forest_ghost.h
 * We define the ghost routine to create a layer of halo elements
 * for a forest of trees in this file.
 */

/* TODO: begin documenting this file: make doxygen 2>&1 | grep t8_forest_ghost */

#ifndef T8_FOREST_GHOST_H
#define T8_FOREST_GHOST_H

#include <t8.h>
#include <t8_forest/t8_forest_types.h>

T8_EXTERN_C_BEGIN ();
/* TODO: comment */
void                t8_forest_ghost_init (t8_forest_ghost_t * pghost);

/** Increase the reference count of a ghost structure.
 * \param [in,out]  ghost     On input, this ghost structure must exist with
 *                            positive reference count.
 */
void                t8_forest_ghost_ref (t8_forest_ghost_t ghost);

/** Descrease the reference count of a ghost structure.
 * If the counter reaches zero, the ghost structure is destroyed.
 * See also \ref t8_forest_ghost_destroy, which is to be preferred when it is
 * known that the last reference to a cmesh is deleted.
 * \param [in,out]  pghost      On inputt, the ghost structure pointed to must
 *                              exist with positive reference count.
 *                              If the reference count reaches zero, the ghost
 *                              structure is destroyed and this pointer is set
 *                              to NULL.
 *                              Otherwise, the pointer is not changed.
 */
void                t8_forest_ghost_unref (t8_forest_ghost_t * pghost);

/** Verify that a ghost structure has only one reference left and destroy it.
 * This function is preferred over \ref t8_ghost_unref when it is known
 * that the last reference is to be deleted.
 * \param [in,out]  pghost     This ghost structure must have a reference count of one.
 *                             It can be in any state (committed or not).
 *                             Then it effectively calls \ref t8_forest_ghost_unref.
 */
void                t8_forest_ghost_destroy (t8_forest_ghost_t * pghost);

/* TODO: Document */
void                t8_forest_ghost_create (t8_forest_t forest);

T8_EXTERN_C_END ();

#endif /* !T8_FOREST_GHOST_H! */
