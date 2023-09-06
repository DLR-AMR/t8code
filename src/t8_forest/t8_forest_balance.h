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

/** \file t8_forest_balance.h
 * We define the balance routine to establish a 1:2 balance among the elements
 * in a forest.
 */

/* TODO: begin documenting this file: make doxygen 2>&1 | grep t8_forest_balance */

#ifndef T8_FOREST_BALANCE_H
#define T8_FOREST_BALANCE_H

#include <t8.h>
#include <t8_forest/t8_forest_types.h>

T8_EXTERN_C_BEGIN ();

/* TODO: document
 * only temporary and will be replaced in future */
void
t8_forest_balance (t8_forest_t forest, int repartition);

/* Check whether the local elements of a forest are balanced. */
int
t8_forest_is_balanced (t8_forest_t forest);

T8_EXTERN_C_END ();

#endif /* !T8_FOREST_BALANCE_H! */
