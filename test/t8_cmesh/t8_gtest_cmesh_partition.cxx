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
#include <t8_schemes/t8_default/t8_default.hxx>
#include "t8_cmesh/t8_cmesh_trees.h"
#include "t8_cmesh/t8_cmesh_partition.h"
#include <test/t8_gtest_macros.hxx>
#include <test/t8_cmesh_generator/t8_cmesh_example_sets.hxx>
#include <test/t8_gtest_schemes.hxx>

/* We create a cmesh, partition it and repartition it several times.
 * At the end we result in the same partition as at the beginning and we
 * compare this cmesh with the initial one. If they are equal the test is
 * passed.
 */

class t8_cmesh_partition_class:
  public testing::TestWithParam<std::tuple<std::tuple<int, t8_eclass_t>, cmesh_example_base *>> {
 protected:
  void
  SetUp () override
  {
    const int scheme_id = std::get<0> (std::get<0> (GetParam ()));
    scheme = create_from_scheme_id (scheme_id);
    eclass = std::get<1> (std::get<0> (GetParam ()));
    size_t found = std::get<1> (GetParam ())->name.find (std::string ("empty"));
    if (found != std::string::npos) {
      /* Tests not working for empty cmeshes */
      GTEST_SKIP ();
    }

    cmesh_original = std::get<1> (GetParam ())->cmesh_create ();
  }
  t8_cmesh_t cmesh_original;
  const t8_scheme *scheme;
  t8_eclass_t eclass;
};

static void
test_cmesh_committed (t8_cmesh_t cmesh)
{
  ASSERT_TRUE (t8_cmesh_is_committed (cmesh)) << "Cmesh commit failed.";
  ASSERT_TRUE (t8_cmesh_trees_is_face_consistent (cmesh, cmesh->trees)) << "Cmesh face consistency failed.";
}

TEST_P (t8_cmesh_partition_class, test_cmesh_partition_concentrate)
{

  const int level = 11;
  int mpisize;
  int mpiret;
  int mpirank;
  t8_cmesh_t cmesh_partition;
  t8_cmesh_t cmesh_partition_new1;
  t8_cmesh_t cmesh_partition_new2;
  t8_shmem_array_t offset_concentrate;

  test_cmesh_committed (cmesh_original);

  mpiret = sc_MPI_Comm_size (sc_MPI_COMM_WORLD, &mpisize);
  SC_CHECK_MPI (mpiret);
  /* Set up the partitioned cmesh */
  for (int i = 0; i < 2; i++) {
    t8_cmesh_init (&cmesh_partition);
    t8_cmesh_set_derive (cmesh_partition, cmesh_original);
    /* Uniform partition according to level */
    t8_cmesh_set_partition_uniform (cmesh_partition, level, scheme);
    t8_cmesh_commit (cmesh_partition, sc_MPI_COMM_WORLD);

    test_cmesh_committed (cmesh_partition);
    cmesh_original = cmesh_partition;
  }

  mpiret = sc_MPI_Comm_rank (sc_MPI_COMM_WORLD, &mpirank);
  SC_CHECK_MPI (mpiret);
  /* Since we want to repartition the cmesh_partition in each step,
   * we need to ref it. This ensure that we can still work with it after
   * another cmesh is derived from it. 
   */
  t8_cmesh_ref (cmesh_partition);

  cmesh_partition_new1 = cmesh_partition;

  /* We repartition the cmesh to be concentrated on each rank once */
  for (int irank = 0; irank < mpisize; irank++) {
    t8_cmesh_init (&cmesh_partition_new2);
    t8_cmesh_set_derive (cmesh_partition_new2, cmesh_partition_new1);
    /* Create an offset array where each tree resides on irank */
    offset_concentrate
      = t8_cmesh_offset_concentrate (irank, sc_MPI_COMM_WORLD, t8_cmesh_get_num_trees (cmesh_partition));
    /* Set the new cmesh to be partitioned according to that offset */
    t8_cmesh_set_partition_offsets (cmesh_partition_new2, offset_concentrate);
    /* Commit the cmesh and test if successful */
    t8_cmesh_commit (cmesh_partition_new2, sc_MPI_COMM_WORLD);
    test_cmesh_committed (cmesh_partition_new2);

    /* Switch the rolls of the cmeshes */
    cmesh_partition_new1 = cmesh_partition_new2;
    cmesh_partition_new2 = NULL;
  }
  /* We partition the resulting cmesh according to a uniform level refinement.
   * This cmesh should now be equal to the initial cmesh. 
   */
  for (int i = 0; i < 2; i++) {
    t8_cmesh_init (&cmesh_partition_new2);
    t8_cmesh_set_derive (cmesh_partition_new2, cmesh_partition_new1);
    t8_cmesh_set_partition_uniform (cmesh_partition_new2, level, t8_scheme_new_default ());
    t8_cmesh_commit (cmesh_partition_new2, sc_MPI_COMM_WORLD);
    cmesh_partition_new1 = cmesh_partition_new2;
  }
  ASSERT_TRUE (t8_cmesh_is_equal (cmesh_partition_new2, cmesh_partition)) << "Cmesh equality check failed.";
  t8_cmesh_destroy (&cmesh_partition_new2);
  t8_cmesh_destroy (&cmesh_partition);
}

/* Test all cmeshes over all different inputs we get through their id */
INSTANTIATE_TEST_SUITE_P (t8_gtest_cmesh_partition, t8_cmesh_partition_class, AllCmeshsParam,
                          pretty_print_base_example);
