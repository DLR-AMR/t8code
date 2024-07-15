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

/* This test program performs checks with the t8_element_transform_face
 * routine.
 * The transformation need to satisfy certrain rules, for example 3 times
 * transforming with an orientation of 1 results in the identity.
 * We check whether such rules are fulfilled.
 */

#include <gtest/gtest.h>
#include <test/t8_gtest_custom_assertion.hxx>

#include <t8_eclass.h>
#include <t8_cmesh.h>
#include <t8_cmesh/t8_cmesh_examples.h>
#include <t8_forest/t8_forest_general.h>
#include <t8_forest/t8_forest_types.h>
#include <t8_schemes/t8_default/t8_default.hxx>
#include <t8_forest/t8_forest_partition.h>
#include <t8_forest/t8_forest_private.h>

class forest_transform: public testing::TestWithParam<std::tuple<t8_eclass, int>> {
 protected:
  void
  SetUp () override
  {
    eclass = std::get<0> (GetParam ());
    level = std::get<1> (GetParam ());

    t8_debugf ("\n\n\nTesting eclass %s with level %i", t8_eclass_to_string[eclass], level);
    default_scheme = t8_scheme_new_default_cxx ();
    /* Construct a coarse mesh of one tree */
    cmesh = t8_cmesh_new_from_class (eclass, sc_MPI_COMM_WORLD);

    /* Create a uniform forest */
    t8_forest_init (&forest);
    t8_forest_set_level (forest, level);
    t8_forest_set_cmesh (forest, cmesh, sc_MPI_COMM_WORLD);
    t8_forest_set_scheme (forest, default_scheme);
    t8_forest_commit (forest);
  }
  void
  TearDown () override
  {
    t8_forest_unref (&forest);
  }
  t8_eclass_t eclass;
  t8_cmesh_t cmesh;
  t8_scheme_cxx_t *default_scheme;
  t8_forest_t forest;
  int level;
};

static void
t8_test_transform_element (t8_eclass_scheme_c *ts, const t8_element_t *elem, t8_eclass_t eclass)
{
  t8_element_t *transform;

  ts->t8_element_new (1, &transform);

  ts->t8_element_transform_face (elem, transform, 0, 0, 0);
  EXPECT_ELEM_EQ (ts, elem, transform);
  ts->t8_element_transform_face (elem, transform, 0, 0, 1);
  EXPECT_ELEM_EQ (ts, elem, transform);
  if (eclass == T8_ECLASS_TRIANGLE) {
    /* For triangles we test:
     * 3 times ori = 1 sign = 0  == identity
     * ori = 1 sign = 0, then ori = 2 sign = 0  == identity
     * ori = 2 sign = 0, then ori = 1 sign = 0  == identity
     * ori = 1 sign = 1, then ori = 1 sign = 1  == identity
     * ori = 2 sign = 1, then ori = 2 sign = 1  == identity
     */
    ts->t8_element_copy (elem, transform);
    /* 3 time or = 1 sign = 0 */
    for (int itimes = 0; itimes < 3; itimes++) {
      ts->t8_element_transform_face (transform, transform, 1, 0, 0);
    }

    EXPECT_ELEM_EQ (ts, elem, transform);
    /* or = 1 sign = 0, then or = 2 sign = 0 */
    ts->t8_element_transform_face (transform, transform, 1, 0, 0);
    ts->t8_element_transform_face (transform, transform, 2, 0, 0);
    EXPECT_ELEM_EQ (ts, elem, transform);
    /* or = 2 sign = 0, then or = 1 sign = 0 */
    ts->t8_element_transform_face (transform, transform, 2, 0, 0);
    ts->t8_element_transform_face (transform, transform, 1, 0, 0);
    EXPECT_ELEM_EQ (ts, elem, transform);
    /* or = 1 sign = 1, then or = 1 sign = 1 */
    ts->t8_element_transform_face (transform, transform, 1, 1, 0);
    ts->t8_element_transform_face (transform, transform, 1, 1, 0);
    EXPECT_ELEM_EQ (ts, elem, transform);
    /* or = 2 sign = 1, then or = 2 sign = 1 */
    ts->t8_element_transform_face (transform, transform, 2, 1, 0);
    ts->t8_element_transform_face (transform, transform, 2, 1, 0);
    EXPECT_ELEM_EQ (ts, elem, transform);
  }
  else {
    T8_ASSERT (eclass == T8_ECLASS_QUAD);
    /* For quads we test:
     * 4 times ori = 1 sign = 0  == identity
     * ori = 1 sign = 0, then ori = 3 sign = 0, then ori = 1 sign = 0 == identity
     * ori = 2 sign = 0, then ori = 2 sign = 0  == identity
     *
     * ori = 1 sign = 1, then ori = 1 sign = 1  == identity
     * ori = 2 sign = 1, then ori = 1 sign = 1  == identity
     */

    ts->t8_element_copy (elem, transform);
    /* 4 times or = 1 sign = 0 */
    for (int itimes = 0; itimes < 4; itimes++) {
      ts->t8_element_transform_face (transform, transform, 1, 0, 1);
    }
    EXPECT_ELEM_EQ (ts, elem, transform);
    /* 4 times or = 1 sign = 0, if not smaller face */
    for (int itimes = 0; itimes < 4; itimes++) {
      ts->t8_element_transform_face (transform, transform, 1, 0, 0);
    }
    EXPECT_ELEM_EQ (ts, elem, transform);
    /* or = 1 sign = 0, then or = 3 sign = 0, then ori = 1 sign = 0 */
    ts->t8_element_transform_face (transform, transform, 1, 0, 1);
    ts->t8_element_transform_face (transform, transform, 3, 0, 1);
    ts->t8_element_transform_face (transform, transform, 1, 0, 1);
    EXPECT_ELEM_EQ (ts, elem, transform);
    /* or = 2 sign = 0, then or = 1 sign = 0 */
    ts->t8_element_transform_face (transform, transform, 2, 0, 1);
    ts->t8_element_transform_face (transform, transform, 1, 0, 1);
    EXPECT_ELEM_EQ (ts, elem, transform);
    /* TODO: Add tests */
  }

  /* Transforming back and forth must lead to the same element */
  for (int iorientation = 0; iorientation < t8_eclass_num_vertices[eclass]; iorientation++) {
    for (int sign = 0; sign < 2; sign++) {
      ts->t8_element_transform_face (elem, transform, iorientation, sign, 1);
      ts->t8_element_transform_face (transform, transform, iorientation, sign, 0);
      EXPECT_ELEM_EQ (ts, elem, transform) << "Orientation " << iorientation << " smaller sign " << sign;
      ts->t8_element_transform_face (elem, transform, iorientation, sign, 0);
      ts->t8_element_transform_face (transform, transform, iorientation, sign, 1);
      EXPECT_ELEM_EQ (ts, elem, transform) << "Orientation " << iorientation << " not smaller sign " << sign;
    }
  }

  ts->t8_element_destroy (1, &transform);
}

TEST_P (forest_transform, test_forest_transform_elements)
{
  for (int ielem = 0; ielem < t8_forest_get_local_num_elements (forest); ielem++) {
    /* Get a pointer to the element */
    t8_element_t *element = t8_forest_get_element (forest, ielem, NULL);
    /* perform the transform test */
    t8_test_transform_element (default_scheme->eclass_schemes[eclass], element, eclass);
  }
}

INSTANTIATE_TEST_SUITE_P (t8_gtest_forest_transform, forest_transform,
                          testing::Combine (testing::Range (T8_ECLASS_QUAD, T8_ECLASS_TRIANGLE),
                                            testing::Range (0, 6)));
