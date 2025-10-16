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

/* In this test we test the t8_forest_set/get_user_data and
 * t8_forest_set/get_user_function functions.
 */

#include <gtest/gtest.h>
#include <t8_schemes/t8_default/t8_default.hxx>
#include <t8_cmesh.h>
#include <t8_cmesh/t8_cmesh_examples.h>
#include <t8_forest/t8_forest_general.h>
#include <test/t8_gtest_schemes.hxx>

/* Test t8_forest_set/get_user_data.
 * We build a forest and set user data for it.
 * We then retrieve the data and check whether it is the same.
 */
class forest_user_data: public testing::TestWithParam<int> {
 protected:
  void
  SetUp () override
  {
    const int scheme_id = GetParam ();
    scheme = create_from_scheme_id (scheme_id);
  }
  const t8_scheme *scheme;
};

TEST_P (forest_user_data, test_user_data)
{
  /* Build a forest */
  t8_cmesh_t cmesh = t8_cmesh_new_hypercube (T8_ECLASS_TRIANGLE, sc_MPI_COMM_WORLD, 0, 0, 0, 0);
  t8_forest_t forest = t8_forest_new_uniform (cmesh, scheme, 1, 0, sc_MPI_COMM_WORLD);
  /* Define user data */
  double data = 42.42;
  double *get_data;
  int data2 = 42;
  int *get_data2;

  /* Set user data of forest */
  t8_forest_set_user_data (forest, (void *) &data);
  /* Get pointer to forest user data */
  get_data = (double *) t8_forest_get_user_data (forest);

  /* Check whether retrieved pointer points to data. */
  ASSERT_EQ (get_data, &data) << "Forest returned wrong user data pointer.";
  /* Check whether data was unchanged. */
  ASSERT_EQ (*get_data, 42.42) << "User data value has changed.";

  /* Add new data and repeat the check */

  /* Set user data of forest */
  t8_forest_set_user_data (forest, (void *) &data2);
  /* Get pointer to forest user data */
  get_data2 = (int *) t8_forest_get_user_data (forest);

  /* Check whether retrieved pointer points to data. */
  ASSERT_EQ (get_data2, &data2) << "Forest returned wrong user data pointer.";
  /* Check whether data was unchanged. */
  ASSERT_EQ (*get_data2, 42) << "User data value has changed.";

  /* Clean up */
  t8_forest_unref (&forest);
}

/* A test function that we can set for the 
 * user function pointer of a forest. */
static double
t8_test_function_42 (int i)
{
  return 42.42 + i;
}

/* A second test function that we can set for the 
 * user function pointer of a forest. */
static void
t8_test_function_second (void)
{
  /* do nothing */
}

/* Test t8_forest_set/get_user_function.
 * We build a forest and set a user function for it.
 * We then retrieve the function and check whether it is the same.
 */
TEST_P (forest_user_data, test_user_function)
{
  /* Build a forest */
  t8_cmesh_t cmesh = t8_cmesh_new_hypercube (T8_ECLASS_TRIANGLE, sc_MPI_COMM_WORLD, 0, 0, 0, 0);
  t8_forest_t forest = t8_forest_new_uniform (cmesh, scheme, 1, 0, sc_MPI_COMM_WORLD);

  double (*funpointer) (int);
  void (*funpointer_second) (void);

  /* Set the t8_test_function_42 as user function pointer. */
  t8_forest_set_user_function (forest, (void (*) (void)) & t8_test_function_42);
  /* Retrieve the function pointer from the forest. */
  funpointer = (double (*) (int)) t8_forest_get_user_function (forest);

  /* Check whether the function pointer is correct. */
  ASSERT_EQ (funpointer, &t8_test_function_42) << "Forest returned wrong user function pointer.";
  /* Additionally call the function and check correct return value. */
  ASSERT_EQ (funpointer (0), 42.42) << "Forest function pointer returned wrong result.";

  /* Overwrite the function user pointer with a second function. */
  t8_forest_set_user_function (forest, (void (*) (void)) & t8_test_function_second);
  /* Retrieve the function pointer from the forest. */
  funpointer_second = (void (*) (void)) t8_forest_get_user_function (forest);

  /* Check whether the function pointer is correct. */
  ASSERT_EQ (funpointer_second, &t8_test_function_second) << "Forest returned wrong user function pointer.";

  /* clean up */
  t8_forest_unref (&forest);
}

INSTANTIATE_TEST_SUITE_P (t8_gtest_ghost_user_data, forest_user_data, AllSchemeCollections);
