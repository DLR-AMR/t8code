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

// testing
#include <gtest/gtest.h>
#include <test/t8_gtest_schemes.hxx>

// t8code
#include <t8.h>
#include <t8_cmesh.h>
#include <t8_cmesh/t8_cmesh_examples.h>
#include <t8_forest/t8_forest_general.h>
#include <t8_forest/t8_forest_geometrical.h>
#include <t8_forest/t8_forest_io.h>
#include <t8_forest/t8_forest_types.h>
#include <t8_schemes/t8_default/t8_default.hxx>
#include <t8_types/t8_vec.h>

/** Data structure to define adaptation criterion: Refine or coarsen based on distance from a given midpoint. */
struct test_adapt_data
{
  double midpoint[3];               /* The midpoint of our sphere. */
  double refine_if_inside_radius;   /* if an element's center is smaller than this value, we refine the element. */
  double coarsen_if_outside_radius; /* if an element's center is larger this value, we coarsen its family. */
};

/** Example callback function used for mesh adaptation */
int
test_adapt_callback (t8_forest_t forest, t8_forest_t forest_from, t8_locidx_t which_tree,
                     [[maybe_unused]] t8_eclass_t tree_class, [[maybe_unused]] t8_locidx_t lelement_id,
                     [[maybe_unused]] const t8_scheme *scheme, const int is_family,
                     [[maybe_unused]] const int num_elements, t8_element_t *elements[])
{
  // Set adaptation criterion.
  const struct test_adapt_data *adapt_data = (const struct test_adapt_data *) t8_forest_get_user_data (forest);

  // Compute the element's centroid coordinates.
  double centroid[3];
  t8_forest_element_centroid (forest_from, which_tree, elements[0], centroid);

  // Compute the distance to our sphere midpoint and decide whether to adapt.
  double dist = t8_dist (centroid, adapt_data->midpoint);
  if (dist < adapt_data->refine_if_inside_radius) {
    return 1;
  }
  else if (is_family && dist > adapt_data->coarsen_if_outside_radius) {
    return -1;
  }
  return 0;
}

/** Test function that adapts the mesh around some hard-coded midpoint.*/
t8_forest_t
test_adapt_forest (t8_forest_t forest)
{
  // Set up an example adaptation criterion
  struct test_adapt_data adapt_data = { { 0.25, 0.5, 0.0 }, 0.125, 0.25 };

  // Adapt forest and return it.
  t8_forest_t forest_adapt = t8_forest_new_adapt (forest, test_adapt_callback, 0, 0, &adapt_data);
  return forest_adapt;
}

/** In this test we test the partition-for-coarsening functionality.
  * The test is based on the following main steps:
  *   - create and adapt an example forest and apply partition for coarsening
  *   - send the last 10 leaf elements of a processor to the next one
  *   - on each process, check whether a family is split across the border to the lower-rank process
  **/
class t8_test_partition_for_coarsening_test: public testing::TestWithParam<std::tuple<int, t8_eclass_t>> {

 public:
  /**
   * Obtain the last elements that of the neighboring process with the lower MPI rank.
   *
   * \details To make sure the forest is partitioned for coarsening, we have to make sure no family is split over
   *          process boundaries. To check this criterion, we need the elements of our neighboring process.
   *          By convention, each process checks its lower element boundary, so that we need the last few elements
   *          of the neighboring process with the lower MPI rank.
   *          This function therefore collects the last MYTODO elements of this process and stores them into the
   *          output array.
   *
   * \param[in]   forest      the considered forest
   * \param[out]  recvElems   array containing the elements
  */
  void
  find_elements_beyond_lower_process_bound (const t8_forest_t forest, t8_element_t **recvElems)
  {

    // (0) Allocation and initialization.
    // ----------------------------------
    //    - Forest has to be committed
    T8_ASSERT (t8_forest_is_committed (forest));

    //    - Determine parent of first local element.
    t8_locidx_t ltree_id;
    t8_element_t *element = t8_forest_get_leaf_element (forest, 0, &ltree_id);

    //    - Set scheme and eclass to the parameters currently tested
    const t8_scheme_c *newscheme = tested_scheme;
    t8_eclass_t eclass = tested_eclass;

    //    - Get parent of first element and determine number of children
    // TODO: Use function get_max_num_children() once it is merged into main
    int max_children = 10;

    //    - Throw error if partition gets too small.
    int num_loc_elems = t8_forest_get_local_num_leaf_elements (forest);
    ASSERT_LE (max_children, num_loc_elems)
      << "This test requires that each process has more elements than max_children";
    t8_debugf ("Packing the last %i elements into array.\n", max_children);

    //    - Allocate array to hold last elements of this process.
    t8_element_t **lastElems = T8_ALLOC (t8_element_t *, max_children);
    newscheme->element_new (eclass, max_children, lastElems);

    // (1) Loop over the last max_children elements of this process.
    // -------------------------------------------------------------
    for (int i = 0; i < max_children; i++) {
      int i_loc_elem = num_loc_elems - max_children + i;

      //  - Write element into array lastElems.
      element = t8_forest_get_leaf_element (forest, i_loc_elem, &ltree_id);
      ASSERT_TRUE (newscheme->element_is_valid (eclass, element));
      newscheme->element_copy (eclass, element, lastElems[i]);
    }

    // (2) Prepare sending and receiving elements via MPI.
    // ---------------------------------------------------
    t8_debugf ("Prepare them for MPI send.\n");

    //    - Determine size of MPI package to send.
    int pack_size;
    size_t count = 1;
    newscheme->element_MPI_Pack_size (eclass, count, forest->mpicomm, &pack_size);
    pack_size *= (max_children);

    //    - Allocate buffers
    int recvBufferSize = pack_size;
    char *sendbuf = T8_ALLOC (char, pack_size);
    char *recvbuf = T8_ALLOC (char, recvBufferSize);

    //    - Pack data
    int position = 0;
    newscheme->element_MPI_Pack (eclass, lastElems, max_children, sendbuf, pack_size, &position, forest->mpicomm);

    // (3) Send elements to neighboring process with higher MPI rank.
    // -------------------------------------------------------------
    //    - Determine process to send to.
    int sendToRank = forest->mpirank + 1;
    if (sendToRank == forest->mpisize)
      sendToRank = 0;
    t8_debugf ("Send data to rank %i...\n", sendToRank);

    //    - Non-blocking MPI send.
    sc_MPI_Request request;
    int mpiret = sc_MPI_Isend (sendbuf, position, sc_MPI_PACKED, sendToRank, T8_MPI_TEST_ELEMENT_PACK_TAG,
                               forest->mpicomm, &request);

    //    - Sanity check and MPI barrier.
    SC_CHECK_MPI (mpiret);
    ASSERT_EQ (pack_size, position);
    sc_MPI_Barrier (forest->mpicomm);

    // (4) Receive data from neighboring process with lower MPI rank.
    // --------------------------------------------------------------
    //    - Determine rank to receive data from.
    int recvFromRank = forest->mpirank - 1;
    if (recvFromRank < 0)
      recvFromRank = forest->mpisize - 1;
    t8_debugf ("Receive data from rank %i.\n", recvFromRank);

    //    - MPI receive
    mpiret = sc_MPI_Recv (recvbuf, pack_size, sc_MPI_PACKED, recvFromRank, T8_MPI_TEST_ELEMENT_PACK_TAG,
                          forest->mpicomm, sc_MPI_STATUS_IGNORE);
    SC_CHECK_MPI (mpiret);

    //    - Finalize MPI communication
    mpiret = sc_MPI_Wait (&request, sc_MPI_STATUS_IGNORE);
    SC_CHECK_MPI (mpiret);
    sc_MPI_Barrier (forest->mpicomm);

    // (5) Unpack received data into elements.
    // ---------------------------------------
    t8_debugf ("Unpack received data into elements.\n");
    position = 0;
    newscheme->element_MPI_Unpack (eclass, recvbuf, recvBufferSize, &position, recvElems, max_children,
                                   forest->mpicomm);
    ASSERT_EQ (position, pack_size);

    //    - Sanity checks: Make sure the received elements are valid
    for (int i = 0; i < max_children; i++) {
      ASSERT_TRUE (newscheme->element_is_valid (eclass, recvElems[i]));
    }

    // (6) Free dynamically allocated memory.
    // --------------------------------------
    t8_element_destroy (newscheme, eclass, max_children, lastElems);
    T8_FREE (lastElems);
    T8_FREE (sendbuf);
    T8_FREE (recvbuf);
  }

  /**
   * This function determines whether the partition of the given forest satisfies the PFC criterion.
   *
   * \details Specifically, it checks whether a family is split at the process boundary to the lower-rank neighbor.
   *
   * \param[in] forest  the forest
   *
   * \return zero if the check is passed by all processes, nonzero otherwise.
  */
  void
  check_if_pfc_partitioned (const t8_forest_t forest)
  {

    // (0) Allocation and initialization.
    // ----------------------------------
    // Get maximum number of children
    // TODO: Use function get_max_num_children() once it is merged into main
    int max_children = 10;

    //    - Set scheme and eclass to the current testing parameters.
    const t8_scheme_c *newscheme = tested_scheme;
    t8_eclass_t eclass = tested_eclass;

    //    - Determine elements beyond lower process bound
    t8_element_t **recvElems = T8_ALLOC (t8_element_t *, max_children);
    newscheme->element_new (eclass, max_children, recvElems);
    find_elements_beyond_lower_process_bound (forest, recvElems);

    //    - Allocate memory for the two parent elements required during the check.
    t8_element_t *start_parent;
    t8_element_new (newscheme, eclass, 1, &start_parent);
    t8_element_t *cur_parent;
    t8_element_new (newscheme, eclass, 1, &cur_parent);

    //    - Get the first local element and its parent
    t8_locidx_t ltree_id;
    t8_element_t *element = t8_forest_get_leaf_element (forest, 0, &ltree_id);
    newscheme->element_get_parent (eclass, element, start_parent);

    //    - Get the size of the potential family, i.e., the number of children of the start element's parent.
    int num_children = newscheme->element_get_num_children (eclass, start_parent);

    // (1) Determine how many siblings are in the beginning of this process' elements.
    // -------------------------------------------------------------------------------
    int sibling_counter = 1;

    //    - Loop over the first num_children elements
    for (int i_loc_elem = 1; i_loc_elem < num_children; i_loc_elem++) {

      //    - Get current element and its parent
      element = t8_forest_get_leaf_element (forest, i_loc_elem, &ltree_id);
      newscheme->element_get_parent (eclass, element, cur_parent);

      //  - Terminate the loop if the current element is no sibling of the first element,
      if (!newscheme->element_is_equal (eclass, cur_parent, start_parent))
        break;

      //  - Increase sibling counter.
      sibling_counter++;
    }

    // (2) Check whether any family is split over the boundary.
    // --------------------------------------------------------
    //    - If the partition does not start with a family, check whether the remaining siblings
    //      are held by the neighboring process.
    if (sibling_counter < num_children) {

      //  - Loop over elements received by neighboring process.
      for (int i_recv_elem = max_children - 1; i_recv_elem >= 0; i_recv_elem--) {

        //  - Get parent of current element
        newscheme->element_get_parent (eclass, recvElems[i_recv_elem], cur_parent);

        //  - Terminate the loop if the current element is no sibling of the first element,
        if (!newscheme->element_is_equal (eclass, cur_parent, start_parent))
          break;

        //  - Increase sibling counter.
        sibling_counter++;

        //  - If the sibling counter reaches num_children, this means a family is split across the process
        //    boundary, so the PFC criteria are not met and we throw an error.
        EXPECT_NE (sibling_counter, num_children)
          << "ERROR: Family is split between processes " << forest->mpirank - 1 << " and " << forest->mpirank;
      }
    }

    // (3) Finalize.
    // -------------
    //    - Deallocate memory
    t8_element_destroy (newscheme, eclass, 1, &cur_parent);
    t8_element_destroy (newscheme, eclass, 1, &start_parent);
    t8_element_destroy (newscheme, eclass, max_children, recvElems);
    T8_FREE (recvElems);
  }

  /** Create an example forest, adapt it, and apply partition for coarsening. */
  t8_forest_t
  create_and_partition_example_forest ()
  {

    // For debugging purpose: Set to true to write vtk output.
    bool debug_write_vtk = false;

    // (1) Create cmesh and initially uniform forest.
    // ----------------------------------------------
    const int level = 4;
    sc_MPI_Comm comm = sc_MPI_COMM_WORLD;
    t8_cmesh_t cmesh = t8_cmesh_new_from_class (tested_eclass, comm);
    t8_forest_t forest = t8_forest_new_uniform (cmesh, tested_scheme, level, 0, comm);
    if (debug_write_vtk)
      t8_forest_write_vtk (forest, "t8_test_pfc_uniform");

    // (2) Adapt and partition the forest - with PFC.
    // ----------------------------------------------
    //    - adapt
    forest = test_adapt_forest (forest);
    if (debug_write_vtk)
      t8_forest_write_vtk (forest, "t8_test_pfc_adapted");

    //    - partition for coarsening
    t8_forest_t pfc_forest;
    t8_forest_init (&pfc_forest);
    int pfc_flag = 1;
    t8_forest_set_partition (pfc_forest, forest, pfc_flag);
    t8_forest_commit (pfc_forest);
    if (debug_write_vtk)
      t8_forest_write_vtk (pfc_forest, "t8_test_pfc_pfc");

    // Barrier and return.
    sc_MPI_Barrier (pfc_forest->mpicomm);
    return pfc_forest;
  }

 protected:
  /** During SetUp, set the scheme and the eclass based on the current testing parameters.*/
  void
  SetUp () override
  {
    const int scheme_id = std::get<0> (GetParam ());
    tested_scheme = create_from_scheme_id (scheme_id);
    tested_eclass = std::get<1> (GetParam ());
  }

  /** During TearDown, decrease the scheme's reference pointer to trigger its destruction.*/
  void
  TearDown () override
  {
    tested_scheme->unref ();
  }

  // Member variables: The currently tested scheme and eclass.
  const t8_scheme *tested_scheme;
  t8_eclass_t tested_eclass;
};

// The test's main function.
TEST_P (t8_test_partition_for_coarsening_test, test_partition_for_coarsening)
{

  // Check eclass: This test case does not make sense for vertices.
  // TODO: Adjust for lines once PR adding get_max_num_children() is merged into main
  if (tested_eclass == T8_ECLASS_VERTEX || tested_eclass == T8_ECLASS_LINE) {
    t8_global_productionf ("Skipping test for this eclass = %i!\n", tested_eclass);
    return;
  }

  // Create an example forest and repartition it using partition for coarsening.
  t8_forest_t forest = create_and_partition_example_forest ();

  // Check whether the partitioning satisfies the PFC criteria.
  check_if_pfc_partitioned (forest);

  // Memory deallocation: Destroy the forest, but keep the scheme.
  tested_scheme->ref ();
  t8_forest_unref (&forest);
}

// Instantiate parameterized test to be run for all schemes.
INSTANTIATE_TEST_SUITE_P (t8_gtest_partititon_for_coarsening, t8_test_partition_for_coarsening_test, AllSchemes,
                          print_all_schemes);
