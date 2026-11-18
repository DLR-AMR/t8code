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

/** \file t8_forest_partition.hxx
 * C++ interface of the forest partition, in addition to \ref t8_forest_partition.h.
 */

#ifndef T8_FOREST_PARTITION_HXX
#define T8_FOREST_PARTITION_HXX

#include <t8_forest/t8_forest_partition.h>

struct t8_abstract_vector_handler;

/** \brief Repartition element data accordingly to a partitioned forest.
 *
 * The data is repartitioned in place: Before the call, \a handler holds one item per local leaf element of
 * \a forest_from, possibly followed by entries for its ghost elements, which are dropped. Afterwards it holds one
 * item per local leaf element of \a forest_to. New storage has room for the ghost elements of \a forest_to, so
 * that appending their entries, e.g. for a ghost exchange, does not reallocate.
 * Items with a bytewise data handler, e.g. trivially copyable types with their default t8_data_handler, are
 * exchanged without packing or intermediate buffers. All other items are serialized by the data handler of
 * \a handler, e.g. one passed to the t8_vector_handler constructor.
 * Collective over the communicator of the forests.
 *
 * \code
 * t8_vector_handler handler (std::move (element_data));
 * t8_forest_partition_data (forest_from, forest_to, handler);
 * element_data = handler.take_data ();
 * \endcode
 *
 * \param[in] forest_from The committed forest before the partitioning step.
 * \param[in] forest_to The committed forest that was partitioned from \a forest_from.
 * \param[in,out] handler The element data, e.g. a t8_vector_handler.
 *
 * \note If the element offsets of a forest have not been created yet, they are created.
 */
void
t8_forest_partition_data (t8_forest_t forest_from, t8_forest_t forest_to, t8_abstract_vector_handler &handler);

#endif /* !T8_FOREST_PARTITION_HXX */
