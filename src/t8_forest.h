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

/** \file t8_forest.h
 * We define the forest of tet-trees in this file.
 */

/* TODO: begin documenting this file: make doxygen 2>&1 | grep t8_forest */

#ifndef T8_FOREST_H
#define T8_FOREST_H

#include <t8_cmesh.h>
#include <t8_element.h>

typedef struct t8_forest *t8_forest_t;

T8_EXTERN_C_BEGIN ();

/** Create a new forest with reference count one.
 * This forest needs to be specialized with the t8_forest_set_* calls.
 * Then in needs to be set up with \see t8_forest_construct.
 * \param [in,out] pforest      On input, this pointer must be non-NULL.
 *                              On return, this pointer set to the new forest.
 */
void                t8_forest_new (t8_forest_t * pforest);

void                t8_forest_set_mpicomm (t8_forest_t forest,
                                           sc_MPI_Comm mpicomm, int do_dup);

/** Set the cmesh associated to a forest.
 *  If the cmesh was not referenced yet, the forest takes ownership of the
 *  cmesh and it will be destroyed when the forest is destroyed.
 *  To reference cmesh, call t8_cmesh_ref before passing it to t8_forest_set_mesh.
 *  \param [in,out] forest      The forest whose cmesh variable will be set.
 *  \param [in]     cmesh       The cmesh to be set.
 */
void                t8_forest_set_cmesh (t8_forest_t forest,
                                         t8_cmesh_t cmesh);

/** Set the scheme associated to a forest.
 *  If the scheme was not referenced yet, the forest takes ownership of the
 *  scheme and it will be destroyed when the forest is destroyed.
 *  To reference scheme, call t8_scheme_ref before passing it to t8_forest_set_mesh.
 *  \param [in,out] forest      The forest whose scheme variable will be set.
 *  \param [in]     scheme       The scheme to be set.
 */
void                t8_forest_set_scheme (t8_forest_t forest,
                                          t8_scheme_t * scheme);

void                t8_forest_set_level (t8_forest_t forest, int level);

/* TODO: by default we take ownership of the 'from' forest.
 *       This means that we call forest_unref (from) in forest_construct.
 *       The caller can keep ownership by calling forest_ref (from) before
 *       passing from into this function.
 */
void                t8_forest_set_copy (t8_forest_t forest,
                                        const t8_forest_t * from);
void                t8_forest_set_adapt (t8_forest_t forest,
                                         const t8_forest_t * from);
void                t8_forest_set_partition (t8_forest_t forest,
                                             const t8_forest_t * from);

void                t8_forest_set_balance (t8_forest_t forest,
                                           int do_balance);
void                t8_forest_set_ghost (t8_forest_t forest, int do_ghost);

/* TODO: use assertions and document that the forest_set (..., from) and
 *       set_load are mutually exclusive. */
void                t8_forest_set_load (t8_forest_t forest,
                                        const char *filename);

/** After allocating and adding properties to a forest, finish its construction.
 * \param [in,out] forest       Must be created with \see t8_forest_new and
 *                              specialized with t8_forest_set_* calls first.
 */
void                t8_forest_construct (t8_forest_t forest);

void                t8_forest_save (t8_forest_t forest);
void                t8_forest_write_vtk (t8_forest_t forest,
                                         const char *filename);

void                t8_forest_iterate (t8_forest_t forest);

/** Increase the reference counter of a forest.
 * \param [in,out] forest       On input, this forest must exist with positive
 *                              reference count.  It may be in any state.
 */
void                t8_forest_ref (t8_forest_t forest);

/** Decrease the reference counter of a forest.
 * If the counter reaches zero, this forest is destroyed.
 * \param [in,out] pforest      On input, the forest pointed to must exist
 *                              with positive reference count.  It may be in
 *                              any state.  If the reference count reaches
 *                              zero, the forest is destroyed and this pointer
 *                              set to NULL.
 *                              Otherwise, the pointer is not changed and
 *                              the forest is not modified in other ways.
 */
void                t8_forest_unref (t8_forest_t * pforest);

T8_EXTERN_C_END ();

#endif /* !T8_FOREST_H */
