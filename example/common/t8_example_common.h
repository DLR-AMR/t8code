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

/** file t8_vtk.h
 * This header file collects macros that are needed for
 * the forest and cmesh vtk routines.
 * \see t8_forest_vtk.h \see t8_cmesh_vtk.h
 */

#ifndef T8_EXAMPLE_COMMON_H
#define T8_EXAMPLE_COMMON_H

#include <t8.h>

T8_EXTERN_C_BEGIN ();

/* function declarations */

/** Adapt a forest such that always the second child of the first
 * tree is refined and no other elements. This results in a highly
 * imbalanced forest.
 */
int                 t8_common_adapt_balance (t8_forest_t forest,
                                             t8_forest_t forest_from,
                                             t8_locidx_t which_tree,
                                             t8_eclass_scheme_c * ts,
                                             int num_elements,
                                             t8_element_t * elements[]);

T8_EXTERN_C_END ();

#endif /* !T8_EXAMPLE_COMMON_H */
