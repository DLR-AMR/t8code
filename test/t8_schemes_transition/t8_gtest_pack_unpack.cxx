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
#include <t8_eclass.h>
#include <test/t8_gtest_custom_assertion.hxx>
#include <test/t8_gtest_macros.hxx>
#include <test/t8_schemes_transition/t8_gtest_dfs_base.hxx>
#include <t8_schemes/t8_transition/t8_transition_cxx.hxx>


/** Use DFS to check for all elements, if packing them, sending them to ourself and unpacking them results in the same element
 * Here, each element is sent individually.
 */
class class_test_pack: public TestDFS {
  virtual void

  /* pack the element and its children, send to ourself, unpack and check if it is the same element */
  check_element ()
  {
    size_t count = 1;
    int position = 0;

    /* Compute pack size and allocate send buffer */
    int pack_size;
    int num_children = ts->t8_element_num_children (element);
    ts->t8_element_MPI_Pack_size (count, comm, &pack_size);
    pack_size *= (num_children + 1);
    char *sendbuf = T8_ALLOC (char, pack_size);

    /* pack data */
    ts->t8_element_MPI_Pack (&element, count, sendbuf, pack_size, &position, comm);
    t8_element_t **children = T8_ALLOC (t8_element_t *, num_children);
    ts->t8_element_new (num_children, children);
    ts->t8_element_children (element, num_children, children);
    ts->t8_element_MPI_Pack (children, num_children, sendbuf, pack_size, &position, comm);

    int recvBufferSize = pack_size;
    char *recvbuf = T8_ALLOC (char, recvBufferSize);
#if T8_ENABLE_MPI
    /* Send data */
    sc_MPI_Request request;
    mpiret = sc_MPI_Isend (sendbuf, position, sc_MPI_PACKED, rank, T8_MPI_TEST_ELEMENT_PACK_TAG, comm, &request);
    SC_CHECK_MPI (mpiret);

    /* Probe size and allocate */
    sc_MPI_Status status;
    sc_MPI_Probe (rank, T8_MPI_TEST_ELEMENT_PACK_TAG, comm, &status);

    /* receive data */
    mpiret
      = sc_MPI_Recv (recvbuf, position, sc_MPI_PACKED, rank, T8_MPI_TEST_ELEMENT_PACK_TAG, comm, sc_MPI_STATUS_IGNORE);
    SC_CHECK_MPI (mpiret);

    /* Finalize non-blocking send communication */
    mpiret = sc_MPI_Wait (&request, sc_MPI_STATUS_IGNORE);
    SC_CHECK_MPI (mpiret);
#else
    /* just copy the data, if we did not compile with MPI*/
    mempcpy (recvbuf, sendbuf, pack_size);
#endif
    /* Unpack data */
    position = 0;
    ts->t8_element_MPI_Unpack (recvbuf, recvBufferSize, &position, &element_compare, count, comm);
    t8_element_t **children_compare = T8_ALLOC (t8_element_t *, num_children);
    ts->t8_element_new (num_children, children_compare);
    ts->t8_element_MPI_Unpack (recvbuf, recvBufferSize, &position, children_compare, num_children, comm);

    /* free buffers */
    T8_FREE (sendbuf);
    T8_FREE (recvbuf);

    /* Check that data was sent and received correctly */
    EXPECT_ELEM_EQ (ts, element, element_compare);
    for (int ichild = 0; ichild < num_children; ichild++) {
      EXPECT_ELEM_EQ (ts, children[ichild], children_compare[ichild]);
    }
    ts->t8_element_destroy (num_children, children);
    ts->t8_element_destroy (num_children, children_compare);
    T8_FREE (children);
    T8_FREE (children_compare);
  }

 protected:
  void
  SetUp () override
  {
    dfs_test_setup ();
    /* Get element and initialize it */
    ts->t8_element_new (1, &element_compare);

    comm = sc_MPI_COMM_WORLD;
    mpiret = sc_MPI_Comm_rank (comm, &rank);
    SC_CHECK_MPI (mpiret);
  }
  void
  TearDown () override
  {
    /* Destroy element */
    ts->t8_element_destroy (1, &element_compare);

    /* Destroy DFS test */
    dfs_test_teardown ();
  }
  t8_element_t *element_compare;
  sc_MPI_Comm comm;
  int rank;
  int mpiret;
};

TEST_P (class_test_pack, test_equal_dfs)
{
#ifdef T8_ENABLE_LESS_TESTS
  const int maxlvl = 4;
#else
  const int maxlvl = 6;
#endif
  check_recursive_dfs_to_max_lvl (maxlvl);
}

INSTANTIATE_TEST_SUITE_P (t8_gtest_test_all_imps, class_test_pack, testing::Values (T8_ECLASS_HEX));
