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

#include <t8_data/t8_shmem.h>

static void
t8_test_shmem_init_finalize (sc_MPI_Comm mpic)
{
  t8_global_productionf ("Starting shmem init test.\n");
  sc_MPI_Comm         intranode, internode;
  int                 intrarank, interrank;
  int                 intrasize, intersize;
  int                 mpiret;

  /* setup shared memory usage */
  t8_shmem_init (mpic);

  /* Get intranode and internode comm */
  sc_mpi_comm_get_node_comms (mpic, &intranode, &internode);

#if T8_ENABLE_MPI
  /* Check that they are not NULL */
  SC_CHECK_ABORT (intranode != sc_MPI_COMM_NULL,
                  "intra node communicator not set.");
  SC_CHECK_ABORT (internode != sc_MPI_COMM_NULL,
                  "inter node communicator not set.");
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

  t8_productionf ("On intranode comm i am rank %i of %i\n", intrarank,
                  intrasize);
  t8_productionf ("On internode comm i am rank %i of %i\n", interrank,
                  intersize);

  /* finalize shared mem usage */
  t8_shmem_finalize (mpic);

  /* Get intranode and internode comm */
  sc_mpi_comm_get_node_comms (mpic, &intranode, &internode);
  /* Check that they are NULL */
  SC_CHECK_ABORT (intranode == sc_MPI_COMM_NULL,
                  "intra node communicator not set.");
  SC_CHECK_ABORT (internode == sc_MPI_COMM_NULL,
                  "inter node communicator not set.");

}

static void
t8_test_shmem_array (sc_MPI_Comm comm)
{
  t8_global_productionf ("Starting shmem array test.\n");
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

  for (int shmem_type_int = SC_SHMEM_BASIC;
       shmem_type_int < SC_SHMEM_NUM_TYPES; ++shmem_type_int) {
    const sc_shmem_type_t shmem_type = (sc_shmem_type_t) shmem_type_int;
    t8_productionf ("Checking shared memory type %s.\n",
                    sc_shmem_type_to_string[shmem_type]);

    /* setup shared memory usage */
    t8_shmem_init (comm);
    t8_shmem_set_type (comm, shmem_type);

#if T8_ENABLE_MPI
    const sc_shmem_type_t control_shmem_type = sc_shmem_get_type (comm);
    SC_CHECK_ABORT (shmem_type == control_shmem_type,
                    "Setting shmem type not succesful.\n");
#endif

    /* Allocate one integer */
    t8_shmem_array_t    shmem_array;
    t8_shmem_array_init (&shmem_array, element_size, array_length, comm);

    sc_MPI_Comm         check_comm = t8_shmem_array_get_comm (shmem_array);
    SC_CHECK_ABORT (comm == check_comm,
                    "shared memory array has wrong communicator.\n");

    const size_t        check_count =
      t8_shmem_array_get_elem_count (shmem_array);
    SC_CHECK_ABORT (check_count == array_length,
                    "shared memory array has wrong element count.\n");

    const size_t        check_size =
      t8_shmem_array_get_elem_size (shmem_array);
    SC_CHECK_ABORT (check_size == element_size,
                    "shared memory array has wrong element size.\n");

    /* Write into array */
    /* In the first half we use the t8_shmem_array_set_gloidx function,
     * we then use a standard store to write 
     * and at the end we use t8_shmem_array_index_for_writing.
     */
    if (t8_shmem_array_start_writing (shmem_array)) {
      /* Double check that array_length is big enough so that each of the three cases
       * is covered. */
      SC_CHECK_ABORT (0 < array_length / 3
                      && array_length / 3 < (int) (2. / 3 * array_length)
                      && (int) (2. / 3 * array_length) < array_length,
                      "Please choose a larger value for array_length.");
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

    for (int i = 0; i < array_length; ++i) {
      t8_gloidx_t         value = t8_shmem_array_get_gloidx (shmem_array, i);
      SC_CHECK_ABORTF (value == i,
                       "Value at position %i not correct (expected %i, got %li)\n",
                       i, i, value);
    }

    /* Copy */
    t8_shmem_array_t    copy_array;
    t8_shmem_array_init (&copy_array, element_size, array_length, comm);
    t8_debugf ("Copying array.\n");
    t8_shmem_array_copy (copy_array, shmem_array);
    SC_CHECK_ABORT (t8_shmem_array_is_equal (copy_array, shmem_array),
                    "Array are not equal after copy.");

    t8_shmem_array_destroy (&shmem_array);
    t8_shmem_array_destroy (&copy_array);

    t8_shmem_finalize (comm);
  }
}

static void
t8_test_sc_shmem_alloc (sc_MPI_Comm comm)
{
  t8_global_productionf ("Starting shmem alloc test.\n");
  int                 mpirank, mpisize;
  int                 mpiret;
  sc_MPI_Comm         intranode, internode;

  mpiret = sc_MPI_Comm_rank (comm, &mpirank);
  SC_CHECK_MPI (mpiret);
  mpiret = sc_MPI_Comm_size (comm, &mpisize);
  SC_CHECK_MPI (mpiret);

  /* Seed random number generator. */
  srand (0);

  for (int shmem_type_int = SC_SHMEM_BASIC;
       shmem_type_int < SC_SHMEM_NUM_TYPES; ++shmem_type_int) {
    const sc_shmem_type_t shmem_type = (sc_shmem_type_t) shmem_type_int;
    int                 intrasize, intrarank;
    t8_productionf ("Checking shared memory type %s.\n",
                    sc_shmem_type_to_string[shmem_type]);

    /* setup shared memory usage */
    t8_shmem_init (comm);
    t8_shmem_set_type (comm, shmem_type);

#if T8_ENABLE_MPI
    const sc_shmem_type_t control_shmem_type = sc_shmem_get_type (comm);
    SC_CHECK_ABORT (shmem_type == control_shmem_type,
                    "Setting shmem type not succesful.\n");
#endif

    sc_mpi_comm_get_node_comms (comm, &intranode, &internode);
#if T8_ENABLE_MPI
    SC_CHECK_ABORT (intranode != sc_MPI_COMM_NULL,
                    "intra node communicator not set.");
    SC_CHECK_ABORT (internode != sc_MPI_COMM_NULL,
                    "inter node communicator not set.");
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
    t8_productionf ("Allocated %i integers at pos %p.\n", intrasize,
                    (void *) shared_int);

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
        SC_CHECK_ABORT (intrarank == 0,
                        "Type is one of window, window_prescan, bgq or bgq_prescan and this process should not have write permissions.\n");
      }
      t8_productionf ("I have write permissions\n");
      shared_int[0] = random_number;
    }
    sc_shmem_write_end (shared_int, comm);

    SC_CHECK_ABORTF (shared_int[0] == random_number,
                     "shared integer was not assigned correctly at position 0.");

    /* Free the memory */
    sc_shmem_free (t8_get_package_id (), shared_int, comm);
    /* Finalize shmem usage so that we can use a different shmem type
     * in the next loop iteration. */
    t8_shmem_finalize (comm);
  }
}

int
main (int argc, char **argv)
{
  int                 mpiret;

  mpiret = sc_MPI_Init (&argc, &argv);
  SC_CHECK_MPI (mpiret);

  sc_init (sc_MPI_COMM_WORLD, 1, 1, NULL, SC_LP_PRODUCTION);
  p4est_init (NULL, SC_LP_ESSENTIAL);
  t8_init (SC_LP_DEFAULT);

#define T8_TEST_SHMEM_NUM_COMMS 2
  sc_MPI_Comm         test_comms[T8_TEST_SHMEM_NUM_COMMS] =
    { sc_MPI_COMM_WORLD, sc_MPI_COMM_SELF };
  const char         *test_comms_names[T8_TEST_SHMEM_NUM_COMMS] =
    { "comm world", "comm self" };

  for (int icomm = 0; icomm < T8_TEST_SHMEM_NUM_COMMS; ++icomm) {
    sc_MPI_Comm         comm = test_comms[icomm];
    t8_productionf ("Checking communicator %s\n\n", test_comms_names[icomm]);
    t8_test_shmem_init_finalize (comm);
    t8_test_sc_shmem_alloc (comm);
    t8_test_shmem_array (comm);
  }

  sc_finalize ();

  mpiret = sc_MPI_Finalize ();
  SC_CHECK_MPI (mpiret);

  return 0;
}
