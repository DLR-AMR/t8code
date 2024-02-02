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
#include "t8_cmesh/t8_cmesh_trees.h"
#include "t8_cmesh/t8_cmesh_partition.h"
#include <t8_cmesh/t8_cmesh_examples.h>
#include <t8_eclass.h>
#include <test/t8_gtest_macros.hxx>
#include "test/t8_cmesh_generator/t8_gtest_cmesh_generator.hxx"
#include "test/t8_cmesh_generator/t8_gtest_cmesh_creator_base.hxx"
#include <iostream>

/* Test if a cmesh is committed properly and perform the face consistency check. */

class cmesh_copy_equality: public testing::TestWithParam<cmesh_sum_cart_prod> {
 protected:
  void
  SetUp () override
  {
    cmesh_gen = GetParam ();
    cmesh_original = cmesh_gen.get_cmesh ();
    /* Set up the cmesh copy */
    t8_cmesh_init (&cmesh_copy);
    /* We need the original cmesh later, so we ref it */
    t8_cmesh_ref (cmesh_original);
    t8_cmesh_set_derive (cmesh_copy, cmesh_original);
    t8_cmesh_commit (cmesh_copy, sc_MPI_COMM_WORLD);
  }
  void
  TearDown () override
  {
    t8_cmesh_unref (&cmesh_original);
    t8_cmesh_unref (&cmesh_copy);
  }

  t8_cmesh_t cmesh_original;
  t8_cmesh_t cmesh_copy;
  cmesh_sum_cart_prod cmesh_gen;
};

/* Test wheater the original cmaeh and its copy are committed and face consistent. Test will fail, if one of these is false. */
TEST_P (cmesh_copy_equality, check_cmeshes_and_their_trees)
{
  EXPECT_TRUE (t8_cmesh_is_committed (cmesh_original));
  EXPECT_TRUE (t8_cmesh_is_committed (cmesh_copy));
  EXPECT_TRUE (t8_cmesh_trees_is_face_consistent (cmesh_original, cmesh_original->trees));
  EXPECT_TRUE (t8_cmesh_trees_is_face_consistent (cmesh_copy, cmesh_copy->trees));
}

/* Test the equality of the original and copied cmeshs*/
TEST_P (cmesh_copy_equality, check_equality_of_copied_cmesh_with_original)
{

  EXPECT_TRUE (t8_cmesh_is_equal (cmesh_original, cmesh_copy));
}

namespace cmesh_list
{
std::vector<sc_MPI_Comm> my_comms = { sc_MPI_COMM_WORLD, sc_MPI_COMM_SELF };
std::vector<t8_eclass_t> eclasses = { T8_ECLASS_ZERO, T8_ECLASS_LINE, T8_ECLASS_TRIANGLE };

std::function<t8_cmesh_t (t8_eclass_t, sc_MPI_Comm)> new_from_class_wrapper = t8_cmesh_new_from_class;

//cmesh_args_cart_prod<t8_eclass_t, sc_MPI_Comm> cmesh_new_from_class( eclasses, my_comms,
//                                      new_from_class_wrapper);
cart_prod_base *new_form_class_prod
  = (cart_prod_base *) new cmesh_args_cart_prod<decltype (eclasses.begin ()), decltype (my_comms.begin ())> (
    std::make_pair (eclasses.begin (), eclasses.end ()), std::make_pair (my_comms.begin (), my_comms.end ()),
    new_from_class_wrapper);

std::vector<cart_prod_base *> cart_prod_vec = { new_form_class_prod };

cmesh_sum_cart_prod cmesh_sums (cart_prod_vec);
cmesh_sum_cart_prod begin = cmesh_sums.begin ();
cmesh_sum_cart_prod end = cmesh_sums.end ();
cmesh_sum_cart_prod step = cmesh_sums.step ();
}  // namespace cmesh_list

/* Test all cmeshes over all different inputs we get through their id */
INSTANTIATE_TEST_SUITE_P (t8_gtest_cmesh_copy, cmesh_copy_equality,
                          ::testing::Range (cmesh_list::begin, cmesh_list::end, cmesh_list::step));
