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

#ifndef T8_FOREST_H
#define T8_FOREST_H

#include <t8_cmesh.h>
#include <t8_element.h>

typedef struct t8_forest *t8_forest_t;

T8_EXTERN_C_BEGIN ();

void                t8_forest_new (t8_forest_t * pforest);

void                t8_forest_set_mpicomm (t8_forest_t forest,
                                           sc_MPI_Comm mpicomm, int do_dup);
void                t8_forest_set_dimension (t8_forest_t forest,
                                             int dimension);
void                t8_forest_set_level (t8_forest_t forest, int level);
void                t8_forest_set_cmesh (t8_forest_t forest,
                                         t8_cmesh_t cmesh);
void                t8_forest_set_scheme (t8_forest_t forest,
                                          t8_scheme_t * scheme);

void                t8_forest_set_copy (t8_forest_t forest,
                                        const t8_forest_t * from);
void                t8_forest_set_adapt (t8_forest_t forest,
                                         const t8_forest_t * from);
void                t8_forest_set_partition (t8_forest_t forest,
                                             const t8_forest_t * from);
void                t8_forest_set_destroy_from (t8_forest_t forest,
                                                int destroy_from);

void                t8_forest_set_balance (t8_forest_t forest,
                                           int do_balance);
void                t8_forest_set_add_ghost (t8_forest_t forest,
                                             int do_ghost);

void                t8_forest_set_load (t8_forest_t forest,
                                        const char *filename);

void                t8_forest_construct (t8_forest_t forest);

void                t8_forest_save (t8_forest_t forest);
void                t8_forest_write_vtk (t8_forest_t forest,
                                         const char *filename);

void                t8_forest_iterate (t8_forest_t forest);

void                t8_forest_destroy (t8_forest_t * pforest);

T8_EXTERN_C_END ();

#endif /* !T8_FOREST_H */
