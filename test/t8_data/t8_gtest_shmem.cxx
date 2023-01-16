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
#include <t8_data/t8_shmem.h>

/* This Test test the shared memory for mpi. It includes tests for initalization and finalizing, the type, the communicator,
 * the element size and count of the shared memory  */

/* *INDENT-OFF* */
#define T8_TEST_SHMEM_NUM_COMMS 2

class shmem : public testing::TestWithParam<int>{
protected:
  void SetUp() override {
    icomm = GetParam();
    comm = test_comms[icomm];
  }
  int icomm;
  sc_MPI_Comm         comm;
  sc_MPI_Comm         test_comms[T8_TEST_SHMEM_NUM_COMMS] = { sc_MPI_COMM_WORLD, sc_MPI_COMM_SELF };
  const char         *test_comms_names[T8_TEST_SHMEM_NUM_COMMS] =
    { "comm world", "comm self" };
};

TEST_P(shmem, test_shmem_init_finalize){
  sc_MPI_Comm         intranode;
  sc_MPI_Comm         internode;
  int                 intrarank;
  int                 interrank;
  int                 intrasize;
  int                 intersize;
  int                 mpiret;

  /* setup shared memory usage */
  t8_shmem_init (comm);

  /* Get intranode and internode comm */
  sc_mpi_comm_get_node_comms (comm, &intranode, &internode);

#if T8_ENABLE_MPI
  /* Check that they are not NULL */
  EXPECT_TRUE(intranode != sc_MPI_COMM_NULL);
  EXPECT_TRUE(internode != sc_MPI_COMM_NULL);
#endif

  /* Compute ranks and size */
  mpiret = sc_MPI_Comm_size (intranode, &intrasize);
  SC_CHECK_MPI (mpiret);
  mpiret = sc_MPI_Comm_rank (intranode, &intrarank);
  SC_CHECK_MPI (mpiret);
  mpiret = sc_MPI_Comm_size (internode, &intersize);
  SC_CHECK_MPI (mpiret);
  mpiret = sc_MPI_Comm_rank (internode, &interrank);
  SC_CHECK_MPI (mpiret);

  /* finalize shared mem usage */
  t8_shmem_finalize (comm);

  /* Get intranode and internode comm */
  sc_mpi_comm_get_node_comms (comm, &intranode, &internode);
  /* Check that they are NULL */
  EXPECT_TRUE(intranode == sc_MPI_COMM_NULL);
  EXPECT_TRUE(internode == sc_MPI_COMM_NULL);
}

TEST_P(shmem, test_sc_shmem_alloc){
  int                 mpirank;
  int                 mpisize;
  int                 mpiret;
  sc_MPI_Comm         intranode;
  sc_MPI_Comm         internode;

  mpiret = sc_MPI_Comm_rank (comm, &mpirank);
  SC_CHECK_MPI (mpiret);
  mpiret = sc_MPI_Comm_size (comm, &mpisize);
  SC_CHECK_MPI (mpiret);

  /* Seed random number generator. */
  srand (0);

  /* Checking shared memory type */
  for (int shmem_type_int = SC_SHMEM_BASIC;
       shmem_type_int < SC_SHMEM_NUM_TYPES; ++shmem_type_int) {
    const sc_shmem_type_t shmem_type = (sc_shmem_type_t) shmem_type_int;
    int                 intrasize, intrarank;

    /* setup shared memory usage */
    t8_shmem_init (comm);
    t8_shmem_set_type (comm, shmem_type);

#if T8_ENABLE_MPI
    const sc_shmem_type_t control_shmem_type = sc_shmem_get_type (comm);
    EXPECT_TRUE(shmem_type == control_shmem_type);
#endif

    sc_mpi_comm_get_node_comms (comm, &intranode, &internode);
#if T8_ENABLE_MPI
    EXPECT_TRUE(intranode != sc_MPI_COMM_NULL);
    EXPECT_TRUE(internode != sc_MPI_COMM_NULL);
#endif

    /* Get the size and rank on the shared memory region (intranode) */
    mpiret = sc_MPI_Comm_size (intranode, &intrasize);
    SC_CHECK_MPI (mpiret);
    mpiret = sc_MPI_Comm_rank (intranode, &intrarank);
    SC_CHECK_MPI (mpiret);

    /* Allocate one integer */
    int                *shared_int =
      (int *) sc_shmem_malloc (t8_get_package_id (), sizeof (int), intrasize,
                               comm);

    /* Write something into the array */
    const int           random_number = rand ();
    if (sc_shmem_write_start (shared_int, comm)) {
      if (
#if defined(SC_ENABLE_MPIWINSHARED)
           shmem_type == SC_SHMEM_WINDOW
           || shmem_type == SC_SHMEM_WINDOW_PRESCAN ||
#endif
#if defined(__bgq__)
           shmem_type == SC_SHMEM_BGQ || shmem_type == SC_SHMEM_BGQ_PRESCAN ||
#endif
           0                    /* If neither SC_ENABLE_MPIWINSHARED or __bgq__ are defined, we neet a condition in this if
                                 * but do not want to execute any code. Hence, 0. */
        ) {
        EXPECT_EQ(intrarank, 0);
      }
      shared_int[0] = random_number;
    }
    sc_shmem_write_end (shared_int, comm);

    EXPECT_EQ(shared_int[0], random_number);

    /* Free the memory */
    sc_shmem_free (t8_get_package_id (), shared_int, comm);
    /* Finalize shmem usage so that we can use a different shmem type
     * in the next loop iteration. */
    t8_shmem_finalize (comm);
  }
}

TEST_P(shmem, test_shmem_array){
  const int           array_length = 100;
  const int           element_size = sizeof (t8_gloidx_t);
  int                 mpirank, mpisize;
  int                 mpiret;

  mpiret = sc_MPI_Comm_rank (comm, &mpirank);
  SC_CHECK_MPI (mpiret);
  mpiret = sc_MPI_Comm_size (comm, &mpisize);
  SC_CHECK_MPI (mpiret);

  /* Seed random number generator. */
  srand (0);

  /* Checking shared memory type */
  for (int shmem_type_int = SC_SHMEM_BASIC;
       shmem_type_int < SC_SHMEM_NUM_TYPES; ++shmem_type_int) {
    const sc_shmem_type_t shmem_type = (sc_shmem_type_t) shmem_type_int;
    
    /* setup shared memory usage */
    t8_shmem_init (comm);
    t8_shmem_set_type (comm, shmem_type);

#if T8_ENABLE_MPI
    const sc_shmem_type_t control_shmem_type = sc_shmem_get_type (comm);
    EXPECT_EQ(shmem_type, control_shmem_type);
#endif

    /* Allocate one integer */
    t8_shmem_array_t    shmem_array;
    t8_shmem_array_init (&shmem_array, element_size, array_length, comm);

    sc_MPI_Comm         check_comm = t8_shmem_array_get_comm (shmem_array);
    /* Check communicator of shared memory array. */
    EXPECT_EQ(comm, check_comm);

    /* Check element count of shared memory array. */
    const size_t        check_count =
      t8_shmem_array_get_elem_count (shmem_array);
    EXPECT_TRUE(check_count == array_length);

    /* Check element size of shared memory array. */
    const size_t        check_size =
      t8_shmem_array_get_elem_size (shmem_array);
    EXPECT_TRUE(check_size == element_size);

    /* Write into array */
    /* In the first half we use the t8_shmem_array_set_gloidx function,
     * we then use a standard store to write 
     * and at the end we use t8_shmem_array_index_for_writing.
     */
    if (t8_shmem_array_start_writing (shmem_array)) {
      /* Double check that array_length is big enough so that each of the three cases
       * is covered. */
      EXPECT_TRUE(0 < array_length / 3
                      && array_length / 3 < (int) (2. / 3 * array_length)
                      && (int) (2. / 3 * array_length) < array_length);

      t8_gloidx_t        *array =
        t8_shmem_array_get_gloidx_array_for_writing (shmem_array);
      for (int i = 0; i < array_length / 3; ++i) {
        t8_shmem_array_set_gloidx (shmem_array, i, i);
      }
      for (int i = array_length / 3; i < 2. / 3 * array_length; ++i) {
        array[i] = i;
      }
      for (int i = 2. / 3 * array_length; i < array_length; ++i) {
        t8_gloidx_t        *index =
          (t8_gloidx_t *) t8_shmem_array_index_for_writing (shmem_array, i);
        *index = i;
      }
    }
    t8_shmem_array_end_writing (shmem_array);

    /* Check value at each position */
    for (int i = 0; i < array_length; ++i) {
      t8_gloidx_t         value = t8_shmem_array_get_gloidx (shmem_array, i);
      EXPECT_TRUE(value == i);
    }

    /* Copy */
    t8_shmem_array_t    copy_array;
    t8_shmem_array_init (&copy_array, element_size, array_length, comm);
    t8_shmem_array_copy (copy_array, shmem_array);
    /* Check equality of arrays after copying. */
    EXPECT_TRUE(t8_shmem_array_is_equal (copy_array, shmem_array));

    t8_shmem_array_destroy (&shmem_array);
    t8_shmem_array_destroy (&copy_array);

    t8_shmem_finalize (comm);
  }
}

INSTANTIATE_TEST_SUITE_P(t8_gtest_shmem, shmem, testing::Range(0, T8_TEST_SHMEM_NUM_COMMS));
/* *INDENT-ON* */
