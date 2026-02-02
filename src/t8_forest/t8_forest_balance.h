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

#ifndef T8_FOREST_BALANCE_H
#define T8_FOREST_BALANCE_H

#include <t8.h>
#include <t8_forest/t8_forest_general.h>

/* We want to export the whole implementation to be callable from "C" */
T8_EXTERN_C_BEGIN ();

/**
 * Balance the forest.
 * 
 * This function adjust the forest such that it satisfies a 2:1 balance condition: The levels of neighboring 
 * elements may differ by at most one. Such a balanced forest simplifies many algorithms, making balancing 
 * an important procedure for both developers and external users.
 * 
 * \param[in,out] forest      The forest to be balanced.
 * \param[in]     repartition Switch deciding whether the forest is repartitioned after balancing.
 */
void
t8_forest_balance (t8_forest_t forest, int repartition);

/**
 * Check whether the local elements of a forest are balanced. 
 *
 * \param[in,out] forest The forest to be checked.
 * \return               1 if the local elements are balanced, 0 otherwise.
 */
int
t8_forest_is_balanced (t8_forest_t forest);

T8_EXTERN_C_END ();

#endif /* !T8_FOREST_BALANCE_H */
