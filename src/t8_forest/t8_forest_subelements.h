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

/** \file t8_forest_eliminate_hanging_nodes.h
 * We define the eliminate_hanging_nodes routine to transform a 1:2 balanced, nonconformal forest
 * into a conformal forest. The routine relies on a 2D quad-scheme that has been balanced, such that 
 * there is a 1:2 balance between all elements.
 */

/* TODO: begin documenting this file: make doxygen 2>&1 | grep t8_forest_balance */

#ifndef T8_FOREST_SUBELEMENTS_H
#define T8_FOREST_SUBELEMENTS_H

#include <t8.h>
#include <t8_forest/t8_forest_types.h>

T8_EXTERN_C_BEGIN ();

int                 t8_forest_subelements_adapt (t8_forest_t forest, t8_forest_t forest_from,
                                                 t8_locidx_t ltree_id, t8_locidx_t lelement_id,
                                                 t8_eclass_scheme_c * ts,
                                                 int num_elements, t8_element_t * elements[]);

void                t8_forest_subelements (t8_forest_t forest);

/* Check whether the local elements of a forest are balanced. */
int                 t8_forest_subelements_used (t8_forest_t forest);

T8_EXTERN_C_END ();

#endif /* !T8_FOREST_BALANCE_H! */
