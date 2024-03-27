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

/** \file t8_forest_partition.h
 * We define the partition routine to partition a forest of trees in this file.
 */

/* TODO: begin documenting this file: make doxygen 2>&1 | grep t8_forest_partition */

#ifndef T8_FOREST_PARTITION_H
#define T8_FOREST_PARTITION_H

#include <t8.h>
#include <t8_forest/t8_forest_general.h>

T8_EXTERN_C_BEGIN ();
/* TODO: document */
void
t8_forest_partition (t8_forest_t forest);

/** Create the element_offset array of a partitioned forest.
 * \param [in,out]  forest The forest.
 * \a forest must be committed before calling this function.
 */
void
t8_forest_partition_create_offsets (t8_forest_t forest);

/** If \ref t8_forest_partition_create_offsets was already called,
 * compute for a given rank the next greater rank that is not empty.
 * \param [in]      forest The forest.
 * \param [in]      rank   An MPI rank.
 * \return                 A rank q > \a rank such that the forest has
 *                         elements on \a q. If such a \a q does not exist,
 *                         returns mpisize.
 */
int
t8_forest_partition_next_nonempty_rank (t8_forest_t forest, int rank);

/** Create the array of global_first_descendant ids of a partitioned forest.
 * \param [in,out]  forest The forest.
 * \a forest must be committed before calling this function.
 */
void
t8_forest_partition_create_first_desc (t8_forest_t forest);

/** Create the array tree offsets of a partitioned forest.
 * This arrays stores at position p the global id of the first tree of this process.
 * Or if this tree is shared, it stores -(global_id) - 1.
 * \param [in,out]  forest  The forest.
 * \a forest must be committed before calling this function.
 */
void
t8_forest_partition_create_tree_offsets (t8_forest_t forest);

/* TODO: document */
/* data_in has length forest_from->num_local_elements
 * data_out   --  --  forest_to->num_local_elements
 */
void
t8_forest_partition_data (t8_forest_t forest_from, t8_forest_t forest_to, const sc_array_t *data_in,
                          sc_array_t *data_out);

/** Test if the last descendant of the last element of current rank has
 * a smaller linear id than the stored first descendant of rank+1.
 * If this is not the case, elements overlap.
 * \param [in]  forest  The forest.
 * \note \a forest must be committed before calling this function.
 */
void
t8_forest_partition_test_boundary_element (const t8_forest_t forest);

T8_EXTERN_C_END ();

#endif /* !T8_FOREST_PARTITION_H */
