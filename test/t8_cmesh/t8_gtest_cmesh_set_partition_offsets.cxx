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

#include <gtest/gtest.h>
#include <t8_cmesh.h>
#include <t8_schemes/t8_default/t8_default_cxx.hxx>
#include "t8_cmesh/t8_cmesh_trees.h"
#include "t8_cmesh/t8_cmesh_partition.h"
#include <t8_cmesh/t8_cmesh_testcases.h>

/* We create a cmesh, partition it and repartition it several times.
 * At the end we result in the same partition as at the beginning and we
 * compare this cmesh with the initial one. If they are equal the test is
 * passed.
 */

TEST (t8_cmesh_set_partition_offsets, test_set_offsets)
{
  t8_cmesh_t cmesh;
  const t8_gloidx_t num_trees = 1;
  const int main_process = 0;
  /* Build a valid offset array. For this test it is onlt necessary that 
   * the array corresponds to any valid partition.
   * We use the offset_concentrate function to build an offset array for a partition
   * that concentrates all trees at one process. */
  t8_shmem_init (sc_MPI_COMM_WORLD);
  t8_shmem_array_t shmem_array = t8_cmesh_offset_concentrate (main_process, sc_MPI_COMM_WORLD, num_trees);

  /* Initialize the cmesh */
  t8_cmesh_init (&cmesh);

  /* Set the partition offsets */
  t8_cmesh_set_partition_offsets (cmesh, shmem_array);
  /* Destroy the cmesh */
  t8_cmesh_unref (&cmesh);
}
