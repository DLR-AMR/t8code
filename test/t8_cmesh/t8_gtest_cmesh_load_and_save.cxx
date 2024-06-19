/*
This file is part of t8code.
t8code is a C library to manage a collection (a forest) of multiple
connected adaptive space-trees of general element classes in parallel.

Copyright (C) 2024 the developers

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
#include <test/t8_gtest_macros.hxx>
#include <t8_cmesh.h>
#include <t8_schemes/t8_default/t8_default_cxx.hxx>
#include "t8_cmesh/t8_cmesh_trees.h"
#include "t8_cmesh/t8_cmesh_partition.h"

#include <test/t8_gtest_macros.hxx>
#include <test/t8_cmesh_generator/t8_cmesh_example_sets.hxx>

#define T8_TEST_LOAD_AND_SAVE_NUM_TREES 100

/* The tests that do commit the cmesh iterate over eclasses and the number of 
 * tress, hence they have a TestWithParam with eclass and int. */
class cmesh_save_and_load: public testing::TestWithParam<cmesh_example_base *> {
 protected:
  void
  SetUp () override
  {
    std::string name;
    GetParam ()->param_to_string (name);
    size_t found = name.find (std::string ("t8_cmesh_new_bigmesh_"));
    if (found != std::string::npos) {
      skipped = true;
      /* Test not working for cmeshes without geometry*/
      GTEST_SKIP ();
    }
    found = GetParam ()->name.find (std::string ("empty"));
    if (found != std::string::npos) {
      skipped = true;

      /* Tests not working for empty cmeshes */
      GTEST_SKIP ();
    }

    cmesh_original = GetParam ()->cmesh_create ();

    t8_cmesh_init (&cmesh);
    t8_cmesh_ref (cmesh_original);
    t8_cmesh_set_derive (cmesh, cmesh_original);

    t8_shmem_init (comm);

    int mpiret = sc_MPI_Comm_size (comm, &mpisize);
    SC_CHECK_MPI (mpiret);
    mpiret = sc_MPI_Comm_rank (comm, &mpirank);
    SC_CHECK_MPI (mpiret);

    const t8_gloidx_t num_trees = t8_cmesh_get_num_trees (cmesh_original);

    t8_shmem_array_t shmem_array = t8_cmesh_offset_concentrate (mpisize - 1, comm, num_trees);

    /* Set the partition offsets */
    t8_cmesh_set_partition_offsets (cmesh, shmem_array);

    /* Commit the cmesh */
    t8_cmesh_commit (cmesh, comm);

    ASSERT_TRUE (t8_cmesh_is_committed (cmesh));
    const int save_return = t8_cmesh_save (cmesh, fileprefix);

    ASSERT_EQ (save_return, 1);

    snprintf (filename, BUFSIZ, "%s_%04d.cmesh", fileprefix, mpirank);
  }

  void
  TearDown () override
  {
    if (!skipped) {
      /* Destroy the cmesh */
      t8_cmesh_unref (&cmesh);
      t8_cmesh_unref (&cmesh_original);
    }

    //const int removed_file = remove (filename);
    //ASSERT_TRUE(!removed_file);
  }

  t8_cmesh_t cmesh;
  t8_cmesh_t cmesh_original;
  int mpisize = 0;
  int mpirank = 0;
  sc_MPI_Comm comm = sc_MPI_COMM_WORLD;
  const char fileprefix[15] = "cmesh_tmp_save";
  char filename[BUFSIZ];
  bool skipped = false;
};

TEST_P (cmesh_save_and_load, test_load)
{
  t8_cmesh_t compare_cmesh = t8_cmesh_load_and_distribute (fileprefix, mpisize, comm, T8_LOAD_SIMPLE, 1);

  ASSERT_EQ (t8_cmesh_is_committed (compare_cmesh), 1);

  t8_cmesh_t cmesh_part_uniform;

  t8_cmesh_init (&cmesh_part_uniform);

  t8_cmesh_set_derive (cmesh_part_uniform, compare_cmesh);

  t8_cmesh_set_partition_uniform (cmesh_part_uniform, 0, t8_scheme_new_default_cxx ());

  t8_cmesh_commit (cmesh_part_uniform, comm);
  ASSERT_TRUE (t8_cmesh_is_committed (cmesh_part_uniform));

  t8_cmesh_unref (&cmesh_part_uniform);
}

/* Make a test suite that iterates over all classes and a tree count from 0 to the maximum. */
INSTANTIATE_TEST_SUITE_P (t8_cmesh_test_load_and_save, cmesh_save_and_load, AllCmeshsParam, pretty_print_base_example);