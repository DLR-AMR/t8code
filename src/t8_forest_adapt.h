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

/** \file t8_forest_adapt.h
 * We define the adapt routine to refine and corsen a forest of trees in this file.
 */

/* TODO: begin documenting this file: make doxygen 2>&1 | grep t8_forest_adapt */

#ifndef T8_FOREST_ADAPT_H
#define T8_FOREST_ADAPT_H

#include <t8.h>
#include <t8_forest.h>

/* TODO: There is no user_data yet */
/** Callback function prototype to replace one set of elements with another.
 *
 * This is used by the adapt routine when the elements of an existing, valid
 * forest are changed.  The callback allows the user to make changes to newly
 * initialized elements before the elements that they replace are destroyed.
 *
 * \param [in] forest      the forest
 * \param [in] which_tree  the tree containing \a outgoing and \a incoming
 * \param [in] ts          the eclass scheme of the tree
 * \param [in] num_outgoing The number of outgoing elements.
 * \param [in] outgoing     The outgoing elements: after the callback, the
 *                          user_data will be destroyed. (at the current state there is no user data)
 * \param [in] num_incoming The number of incoming elements.
 * \param [in,out] incoming The incoming elements: prior to the callback,
 *                          the user_data is allocated, and the forest_init_t callback,
 *                          if it has been provided, will be called.
 *
 * If an element is being refined, num_outgoing will be 1 and num_incoming will
 * be the number of children, and vice versa if a family is being coarsened.
 */
typedef void        (*t8_forest_replace_t) (t8_forest_t * forest,
                                            t8_topidx_t which_tree,
                                            t8_eclass_scheme_t * ts,
                                            int num_outgoing,
                                            t8_element_t * outgoing[],
                                            int num_incoming,
                                            t8_element_t * incoming[]);

#endif /* !T8_FOREST_ADAPT_H! */
