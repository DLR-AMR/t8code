/*
  This file is part of t8code.
  t8code is a C library to manage a collection (a forest) of multiple
  connected adaptive space-trees of general element classes in parallel.

  Copyright (C) 2026 the developers

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

/** \file t8_forest_ghost.hxx
 * C++ interface of the ghost exchange of element data, in addition to \ref t8_forest_ghost_exchange_data in
 * \ref t8_forest_general.h.
 */

#ifndef T8_FOREST_GHOST_HXX
#define T8_FOREST_GHOST_HXX

#include <t8_forest/t8_forest_general.h>

struct t8_abstract_vector_handler;

/** \brief Exchange the element data of the ghost elements of a forest.
 *
 * Before the call, \a handler holds one item per local leaf element of \a forest, possibly followed by further
 * entries, e.g. outdated items of the ghost elements, which are dropped. Afterwards it holds the items of the local
 * leaf elements followed by those of the ghost elements, received from their owners.
 * Items with a bytewise data handler, e.g. trivially copyable types with their default t8_data_handler, are
 * exchanged without packing. All other items are serialized by the data handler of \a handler, e.g. one passed to
 * the t8_vector_handler constructor. The send buffer is kept in the ghost layer and reused by all exchanges of
 * \a forest. Without ghost layer, the handler only keeps the items of the local leaf elements.
 * If profiling is enabled, the time of the whole exchange is stored as ghost wait time.
 * Collective over the communicator of the forest.
 *
 * \code
 * t8_vector_handler handler (std::move (element_data));
 * t8_forest_ghost_exchange_data (forest, handler);
 * element_data = handler.take_data ();
 * \endcode
 *
 * \param[in] forest The committed forest.
 * \param[in,out] handler The element data, e.g. a t8_vector_handler.
 */
void
t8_forest_ghost_exchange_data (t8_forest_t forest, t8_abstract_vector_handler &handler);

#endif /* !T8_FOREST_GHOST_HXX */
