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
#include <t8_gtest_memory_macros.hxx>

/* This Test test the shared memory for mpi. It includes tests for initialization and finalizing, the type, the communicator,
 * the element size and count of the shared memory  */

#define T8_TEST_SHMEM_NUM_COMMS 2

class shmem: public testing::TestWithParam<std::tuple<int, sc_MPI_Comm, int>> {
 protected:
  void
  SetUp () override
  {
    icomm = std::get<0> (GetParam ());
    comm = std::get<1> (GetParam ());
    shmem_type_int = std::get<2> (GetParam ());
    shmem_type = (sc_shmem_type_t) shmem_type;
  }
  int icomm;
  sc_MPI_Comm comm;
  int shmem_type_int;
  sc_shmem_type_t shmem_type;

  const char *test_comms_names[T8_TEST_SHMEM_NUM_COMMS] = { "comm world", "comm self" };
};

TEST_P (shmem, test_shmem_init_finalize)
{
  sc_MPI_Comm intranode;
  sc_MPI_Comm internode;
  int intrarank;
  int interrank;
  int intrasize;
  int intersize;
  int mpiret;

  /* setup shared memory usage */
  const int intrasize_from_init = t8_shmem_init (comm);
  ASSERT_GT (intrasize_from_init, 0) << "Error in t8_shmem_init. No intranode communicator set.";

  /* Get intranode and internode comm */
  sc_mpi_comm_get_node_comms (comm, &intranode, &internode);

#if T8_ENABLE_MPI
  /* Check that they are not NULL */
  ASSERT_NE (intranode, sc_MPI_COMM_NULL) << "intra node communicator not set.";
  ASSERT_NE (internode, sc_MPI_COMM_NULL) << "inter node communicator not set.";
#endif

  /* Compute ranks and size and print them */
  mpiret = sc_MPI_Comm_size (intranode, &intrasize);
  SC_CHECK_MPI (mpiret);
  mpiret = sc_MPI_Comm_rank (intranode, &intrarank);
  SC_CHECK_MPI (mpiret);
  mpiret = sc_MPI_Comm_size (internode, &intersize);
  SC_CHECK_MPI (mpiret);
  mpiret = sc_MPI_Comm_rank (internode, &interrank);
  SC_CHECK_MPI (mpiret);
  t8_debugf ("On intranode communicator I am %i of %i\n", intrarank, intrasize);
  t8_debugf ("On internode communicator I am rank %i of %i\n", interrank, intersize);

  /* finalize shared mem usage */
  t8_shmem_finalize (comm);

  /* Get intranode and internode comm */
  sc_mpi_comm_get_node_comms (comm, &intranode, &internode);
  /* Check that they are NULL */
  ASSERT_EQ (intranode, sc_MPI_COMM_NULL) << "intra node communicator not set.";
  ASSERT_EQ (internode, sc_MPI_COMM_NULL) << "inter node communicator not set.";
}

TEST_P (shmem, test_sc_shmem_alloc)
{
  int mpirank;
  int mpisize;
  int mpiret;
  sc_MPI_Comm intranode;
  sc_MPI_Comm internode;

  mpiret = sc_MPI_Comm_rank (comm, &mpirank);
  SC_CHECK_MPI (mpiret);
  mpiret = sc_MPI_Comm_size (comm, &mpisize);
  SC_CHECK_MPI (mpiret);

  /* Seed random number generator. */
  srand (0);

  /* Checking shared memory type */
  const sc_shmem_type_t shmem_type = (sc_shmem_type_t) shmem_type_int;
  int intrasize, intrarank;
  t8_debugf ("Checking shared memory type %s.\n", sc_shmem_type_to_string[shmem_type]);

  /* setup shared memory usage */
  const int intranode_size = t8_shmem_init (comm);
  ASSERT_GT (intranode_size, 0) << "Could not initialize shared memory.";

  t8_shmem_set_type (comm, shmem_type);

#if T8_ENABLE_MPI
  const sc_shmem_type_t control_shmem_type = sc_shmem_get_type (comm);
  ASSERT_EQ (shmem_type, control_shmem_type) << "Setting shmem type not successful.";
#endif

  sc_mpi_comm_get_node_comms (comm, &intranode, &internode);
#if T8_ENABLE_MPI
  ASSERT_NE (intranode, sc_MPI_COMM_NULL) << "intra node communicator not set.";
  ASSERT_NE (internode, sc_MPI_COMM_NULL) << "inter node communicator not set.";
#endif

  /* Get the size and rank on the shared memory region (intranode) */
  mpiret = sc_MPI_Comm_size (intranode, &intrasize);
  SC_CHECK_MPI (mpiret);
  mpiret = sc_MPI_Comm_rank (intranode, &intrarank);
  SC_CHECK_MPI (mpiret);

  /* Allocate one integer */
  int *shared_int = (int *) sc_shmem_malloc (t8_get_package_id (), sizeof (int), intrasize, comm);
  t8_debugf ("Allocated %i integers at pos %p.\n", intrasize, (void *) shared_int);

  /* Write something into the array */
  const int random_number = rand ();
  if (sc_shmem_write_start (shared_int, comm)) {
    if (
#if defined(SC_ENABLE_MPIWINSHARED)
      shmem_type == SC_SHMEM_WINDOW || shmem_type == SC_SHMEM_WINDOW_PRESCAN ||
#endif
#if defined(__bgq__)
      shmem_type == SC_SHMEM_BGQ || shmem_type == SC_SHMEM_BGQ_PRESCAN ||
#endif
      0 /* If neither SC_ENABLE_MPIWINSHARED or __bgq__ are defined, we need a condition in this if
                                * but do not want to execute any code. Hence, 0. */
    ) {
      ASSERT_EQ (intrarank, 0) << "Type is one of window, window_prescan, bgq or bgq_prescan and this process should "
                                  "not have write permissions.\n";
    }
    t8_debugf ("I have write permissions\n");
    shared_int[0] = random_number;
  }
  sc_shmem_write_end (shared_int, comm);

  ASSERT_EQ (shared_int[0], random_number) << "shared integer was not assigned correctly at position 0.";

  /* Free the memory */
  sc_shmem_free (t8_get_package_id (), shared_int, comm);
  /* Finalize shmem usage so that we can use a different shmem type
    * in the next loop iteration. */
  t8_shmem_finalize (comm);
}

TEST_P (shmem, test_shmem_array_allgatherv)
{
  const size_t element_size = sizeof (t8_gloidx_t);
  int mpirank;
  int mpisize;
  int mpiret;

  mpiret = sc_MPI_Comm_rank (comm, &mpirank);
  SC_CHECK_MPI (mpiret);
  mpiret = sc_MPI_Comm_size (comm, &mpisize);
  SC_CHECK_MPI (mpiret);

  /* Each process fills an array of size mpirank * 10, where the elements
   * in the arrays increase, such that in the shmem we have a contiguous increase. */
  const int base_size = 10;
  const t8_gloidx_t array_length = (mpirank + 1) * base_size;
  t8_gloidx_t *sendbuf = T8_TESTSUITE_ALLOC_ZERO (t8_gloidx_t, array_length);
  const t8_gloidx_t first_array_value = base_size * mpirank * (mpirank + 1) / 2;
  const int total_size = base_size * mpisize * (mpisize + 1) / 2;

  for (t8_gloidx_t i = 0; i < array_length; i++) {
    sendbuf[i] = first_array_value + i;
  }

  t8_debugf ("Checking shared memory type %s.\n", sc_shmem_type_to_string[shmem_type_int]);

  const sc_shmem_type_t shmem_type = (sc_shmem_type_t) shmem_type_int;

  /* setup shared memory usage */
  const int intranode_size = t8_shmem_init (comm);
  ASSERT_GT (intranode_size, 0) << "Could not initialize shared memory.";
  t8_shmem_set_type (comm, shmem_type);

#if T8_ENABLE_MPI
  const sc_shmem_type_t control_shmem_type = sc_shmem_get_type (comm);
  ASSERT_EQ (shmem_type, control_shmem_type) << "Setting shmem type not successful.";
#endif

  t8_shmem_array_t shmem_array;
  t8_shmem_array_init (&shmem_array, element_size, total_size, comm);

  sc_MPI_Comm check_comm = t8_shmem_array_get_comm (shmem_array);
  /* Check communicator of shared memory array. */
  ASSERT_EQ (comm, check_comm) << "Shared memory array has wrong communicator.";

  /* Check element count of shared memory array. */
  const int check_count = t8_shmem_array_get_elem_count (shmem_array);
  ASSERT_EQ (check_count, total_size) << "shared memory array has wrong element count.";

  /* Check element size of shared memory array. */
  const size_t check_size = t8_shmem_array_get_elem_size (shmem_array);
  ASSERT_EQ (check_size, element_size) << "shared memory has wrong element size.";

  t8_shmem_array_allgatherv ((void *) sendbuf, array_length, T8_MPI_GLOIDX, shmem_array, T8_MPI_GLOIDX, comm);

  for (int i = 0; i < total_size; ++i) {
    const int value = t8_shmem_array_get_gloidx (shmem_array, i);
    ASSERT_EQ (value, i) << "Value at position " << i << " not correct (expected " << i << " got " << value << ")";
  }

  t8_shmem_array_destroy (&shmem_array);
  t8_shmem_finalize (comm);
  T8_TESTSUITE_FREE (sendbuf);
}

TEST_P (shmem, test_shmem_array_prefix)
{
  int mpirank;
  int mpisize;
  int mpiret;

  mpiret = sc_MPI_Comm_rank (comm, &mpirank);
  SC_CHECK_MPI (mpiret);
  mpiret = sc_MPI_Comm_size (comm, &mpisize);
  SC_CHECK_MPI (mpiret);

  /* Checking shared memory type */
  const sc_shmem_type_t shmem_type = (sc_shmem_type_t) shmem_type_int;

  /* setup shared memory usage */
  const int intranode_size = t8_shmem_init (comm);
  ASSERT_GT (intranode_size, 0) << "Could not initialize shared memory.";
  t8_shmem_set_type (comm, shmem_type);

#if T8_ENABLE_MPI
  const sc_shmem_type_t control_shmem_type = sc_shmem_get_type (comm);
  ASSERT_EQ (shmem_type, control_shmem_type) << "Setting shmem type not successful.";
#endif

  /* Allocate one integer */
  t8_shmem_array_t shmem_array;
  t8_shmem_array_init (&shmem_array, sizeof (t8_gloidx_t), mpisize + 1, comm);

  t8_gloidx_t sendbuf = 1;
  t8_shmem_array_prefix ((void *) &sendbuf, shmem_array, 1, T8_MPI_GLOIDX, sc_MPI_SUM, comm);

  /* Check value at each position */
  for (int i = 0; i < mpisize; ++i) {
    t8_gloidx_t value = t8_shmem_array_get_gloidx (shmem_array, i);
    ASSERT_EQ (value, i) << "Value at position " << i << " not correct (expected " << i << " got " << value << ")";
  }

  t8_shmem_array_destroy (&shmem_array);

  t8_shmem_finalize (comm);
}

TEST_P (shmem, test_shmem_array)
{
  const int array_length = 100;
  const int element_size = sizeof (t8_gloidx_t);
  int mpirank, mpisize;
  int mpiret;

  mpiret = sc_MPI_Comm_rank (comm, &mpirank);
  SC_CHECK_MPI (mpiret);
  mpiret = sc_MPI_Comm_size (comm, &mpisize);
  SC_CHECK_MPI (mpiret);

  /* Seed random number generator. */
  srand (0);

  /* Checking shared memory type */
  const sc_shmem_type_t shmem_type = (sc_shmem_type_t) shmem_type_int;

  /* setup shared memory usage */
  const int intranode_size = t8_shmem_init (comm);
  ASSERT_GT (intranode_size, 0) << "Could not initialize shared memory.";
  t8_shmem_set_type (comm, shmem_type);

#if T8_ENABLE_MPI
  const sc_shmem_type_t control_shmem_type = sc_shmem_get_type (comm);
  ASSERT_EQ (shmem_type, control_shmem_type) << "Setting shmem type not successful.";
#endif

  /* Allocate one integer */
  t8_shmem_array_t shmem_array;
  t8_shmem_array_init (&shmem_array, element_size, array_length, comm);

  sc_MPI_Comm check_comm = t8_shmem_array_get_comm (shmem_array);
  /* Check communicator of shared memory array. */
  ASSERT_EQ (comm, check_comm) << "Shared memory array has wrong communicator.";

  /* Check element count of shared memory array. */
  const int check_count = t8_shmem_array_get_elem_count (shmem_array);
  ASSERT_EQ (check_count, array_length) << "shared memory array has wrong element count.";

  /* Check element size of shared memory array. */
  const int check_size = t8_shmem_array_get_elem_size (shmem_array);
  ASSERT_EQ (check_size, element_size) << "shared memory has wrong element size.";

  /* Write into array */
  /* In the first half we use the t8_shmem_array_set_gloidx function,
    * we then use a standard store to write 
    * and at the end we use t8_shmem_array_index_for_writing.
    */
  if (t8_shmem_array_start_writing (shmem_array)) {
    /* Double check that array_length is big enough so that each of the three cases
      * is covered. */
    ASSERT_TRUE (0 < array_length / 3 && array_length / 3 < (int) (2. / 3 * array_length)
                 && (int) (2. / 3 * array_length) < array_length)
      << "Please choose a larger value for array length.";

    t8_gloidx_t *array = t8_shmem_array_get_gloidx_array_for_writing (shmem_array);
    for (int i = 0; i < array_length / 3; ++i) {
      t8_shmem_array_set_gloidx (shmem_array, i, i);
    }
    for (int i = array_length / 3; i < 2. / 3 * array_length; ++i) {
      array[i] = i;
    }
    for (int i = 2. / 3 * array_length; i < array_length; ++i) {
      t8_gloidx_t *index = (t8_gloidx_t *) t8_shmem_array_index_for_writing (shmem_array, i);
      *index = i;
    }
  }
  t8_shmem_array_end_writing (shmem_array);

  /* Check value at each position */
  for (int i = 0; i < array_length; ++i) {
    t8_gloidx_t value = t8_shmem_array_get_gloidx (shmem_array, i);
    ASSERT_EQ (value, i) << "Value at position " << i << " not correct (expected " << i << " got " << value << ")";
  }

  /* Copy */
  t8_shmem_array_t copy_array;
  t8_shmem_array_init (&copy_array, element_size, array_length, comm);
  t8_shmem_array_copy (copy_array, shmem_array);
  /* Check equality of arrays after copying. */
  ASSERT_TRUE (t8_shmem_array_is_equal (copy_array, shmem_array)) << "Arrays are not equal after copy.";

  t8_shmem_array_destroy (&shmem_array);
  t8_shmem_array_destroy (&copy_array);

  t8_shmem_finalize (comm);
}

INSTANTIATE_TEST_SUITE_P (t8_gtest_shmem, shmem,
                          testing::Combine (testing::Range (0, T8_TEST_SHMEM_NUM_COMMS),
                                            testing::Values (sc_MPI_COMM_WORLD, sc_MPI_COMM_SELF),
                                            testing::Range ((int) SC_SHMEM_BASIC, (int) SC_SHMEM_NUM_TYPES)));
