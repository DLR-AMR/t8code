/*
  This file is part of t8code.
  t8code is a C library to manage a collection (a forest) of multiple
  connected adaptive space-trees of general element classes in parallel.

  Copyright (C) 2025 the developers

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

/** \file t8_forest_partition_for_coarsening.h
 * adjust the partition of the receiving forest so that all elements from a family
 * belong to the same process before applying the partition
 */

#ifndef T8_FOREST_PFC_H
#define T8_FOREST_PFC_H

#include <t8.h>
#include <t8_forest/t8_forest_general.h>

T8_EXTERN_C_BEGIN ();

/** \brief Correct the first element of this rank in the new partition if it splits an element family.
 *
 *  The default partitioning distributes the elements into equally-sized partitions. For coarsening, however,
 *  all elements of a family have to be on the same process in order to be coarsened into their parent element.
 *  If a full family has elements on both sides of the first element of this rank, the first element is moved
 *  behind the family. The price to be paid is a slight deviation from the optimal balance of elements among
 *  processors.
 *
 * \param [in] forest   The forest to partition, with forest->set_from committed.
 * \param [in] first    The global index of the first element of this rank in the new partition.
 * \return The corrected global index of the first element of this rank.
 * \note This function is collective.
*/
t8_gloidx_t
t8_forest_pfc_correct_first (t8_forest_t forest, const t8_gloidx_t first);

T8_EXTERN_C_END ();
#endif /* T8_FOREST_PFC_H */
