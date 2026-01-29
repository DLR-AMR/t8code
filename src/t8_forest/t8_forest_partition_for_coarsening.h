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

/** \brief Correct the partitioning if element families are split accorss process boundaries.
 *
 *  The default partitioning distributes the elements into equally-sized partitions. For coarsening, however,
 *  all elements of a family have to be on the same process in order to be coarsened into their parent element.
 *  This function corrects the partitioning such that no families are split across process boundaries.
 *  The price to be paid is a slight deviation from the optimal balance of elements among processors.
 *
 * \param [in,out] forest   the forest. On input, it has been partitioned into equally-sized element partitions.
 *                                      On output, the partitioning has been adjusted such that no element families
 *                                      are split across the process boundaries.
*/
void
t8_forest_pfc_correction_offsets (t8_forest_t forest);

T8_EXTERN_C_END ();
#endif /* T8_FOREST_PFC_H */
