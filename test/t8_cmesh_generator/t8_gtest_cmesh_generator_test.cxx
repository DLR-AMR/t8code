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

#include "test/t8_cmesh_generator/t8_cmesh_example_sets.hxx"
#include "test/t8_cmesh_generator/t8_cmesh_parametrized_examples/t8_cmesh_new_bigmesh_param.hxx"
#include "test/t8_cmesh_generator/t8_gtest_cmesh_cartestian_product.hxx"

class t8_cmesh_iter: public testing::TestWithParam<cmesh_example_base *> {
 protected:
  void
  SetUp () override
  {
    cmesh_creator = GetParam ();
    cmesh = cmesh_creator->cmesh_create ();
    cmesh_creator->param_to_string (cmesh_param_string);
  }

  void
  TearDown () override
  {
    /* Unref both cmeshes */
    t8_cmesh_unref (&cmesh);
  }
  t8_cmesh_t cmesh;
  cmesh_example_base *cmesh_creator;
  std::string cmesh_param_string;
};

TEST_P (t8_cmesh_iter, cmesh_created)
{
  EXPECT_NE (cmesh, nullptr);
}

TEST_P (t8_cmesh_iter, print_string)
{
  EXPECT_FALSE (cmesh_param_string.empty ());
}

TEST (param_generator, check_num_parameters_combination)
{
  const size_t num_combs_generated = new_bigmesh::cmesh_example->example_all_combination.size ();
  const size_t num_expected
    = cmesh_params::large_mesh.size () * cmesh_params::my_comms.size () * cmesh_params::eclasses.size ();
  EXPECT_EQ (num_combs_generated, num_expected);
}

INSTANTIATE_TEST_SUITE_P (t8_gtest_create_cmeshes, t8_cmesh_iter, AllCmeshsParam, pretty_print_base_example);
