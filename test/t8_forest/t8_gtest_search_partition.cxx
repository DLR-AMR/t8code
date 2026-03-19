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
#include <test/t8_gtest_custom_assertion.hxx>
#include <t8_eclass.h>
#include <t8_cmesh.h>
#include <t8_cmesh/t8_cmesh_examples.h>
#include <t8_forest/t8_forest_general.h>
#include <t8_forest/t8_forest_search/t8_forest_search.hxx>
#include <t8_schemes/t8_default/t8_default.hxx>
#include <test/t8_gtest_macros.hxx>
#include <test/t8_gtest_schemes.hxx>

class forest_search_partition: public testing::TestWithParam<std::tuple<std::tuple<int, t8_eclass>, int>> {
 protected:
  void
  SetUp () override
  {
    const int scheme_id = std::get<0> (std::get<0> (GetParam ()));
    scheme = create_from_scheme_id (scheme_id);
    eclass = std::get<1> (std::get<0> (GetParam ()));
    level = std::get<1> (GetParam ());

    /* Construct a cube coarse mesh */
    cmesh = t8_cmesh_new_hypercube (eclass, sc_MPI_COMM_WORLD, 0, 0, 0);
    /* Build a uniform forest */
    forest = t8_forest_new_uniform (cmesh, scheme, level, 0, sc_MPI_COMM_WORLD);
  }
  void
  TearDown () override
  {
  }
  t8_eclass_t eclass;
  int level;
  t8_cmesh_t cmesh;
  t8_forest_t forest;
  const t8_scheme *scheme;
};

/* A search function that matches all elements.
 * This function verifies that ltreeid, pfirst and plast are sensible.
 */
bool
t8_test_search_partition_all_fn (const t8_forest_t forest, const t8_locidx_t ltreeid, const t8_element_t *element,
                                 const int pfirst, const int plast, void *user_data)
{
  T8_ASSERT (t8_forest_is_committed (forest));
  T8_ASSERT (user_data != NULL);
  int mpisize = *(int *) user_data;
  T8_ASSERT (0 <= pfirst && pfirst <= plast && plast < mpisize);
  T8_ASSERT (0 <= ltreeid && ltreeid < t8_forest_get_num_global_trees (forest));

  return true;
}

TEST_P (forest_search_partition, t8_test_search_partition_all_fn)
{
  sc_MPI_Comm comm = t8_forest_get_mpicomm (forest);
  int mpisize;
  int mpiret = sc_MPI_Comm_size (comm, &mpisize);
  SC_CHECK_MPI (mpiret);

  /* Call search. This search matches all elements. We check each call of the
   * callback for consistency. */
  t8_partition_search<int> search (t8_test_search_partition_all_fn);

  search.update_user_data (&mpisize);
  search.update_forest (forest);
  search.do_search ();

  t8_forest_unref (&forest);
}

#if T8_TEST_LEVEL_INT >= 2
const int maxlvl = 5;
#else
const int maxlvl = 6;
#endif

INSTANTIATE_TEST_SUITE_P (t8_gtest_search_partition, forest_search_partition,
                          testing::Combine (AllSchemes, testing::Range (0, maxlvl)));
