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

/** \file t8_forest_cxx.h
 * We define the forest routines that need access to the
 * c++ element interface.
 */

/* TODO: begin documenting this file: make doxygen 2>&1 | grep t8_forest_cxx */

#ifndef T8_FOREST_CXX_H
#define T8_FOREST_CXX_H

#include <t8.h>
#include <t8_forest.h>

T8_EXTERN_C_BEGIN ();

/* TODO: document with doxygen */

/* For each tree in a forest compute its first and last descendant */
void                t8_forest_compute_desc (t8_forest_t forest);

/* Create the elements on this process given a uniform partition
 * of the coarse mesh. */
void                t8_forest_populate (t8_forest_t forest);

/* return nonzero if the first tree of a forest is shared with a smaller
 * process.
 * This is the case if and only if the first descendant of the first tree that we store is
 * not the first possible descendant of that tree.
 */
int                 t8_forest_first_tree_shared (t8_forest_t forest);

/* Allocate memory for trees and set their values as in from.
 * For each tree allocate enough element memory to fit the elements of from.
 * If copy_elements is true, copy the elements of from into the element memory.
 */
void                t8_forest_copy_trees (t8_forest_t forest,
                                          t8_forest_t from,
                                          int copy_elements);

T8_EXTERN_C_END ();

#endif /* !T8_FOREST_CXX_H! */
