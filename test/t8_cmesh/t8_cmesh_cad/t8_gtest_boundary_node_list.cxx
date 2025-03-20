/*
  This file is part of t8code.
  t8code is a C library to manage a collection (a forest) of multiple
  connected adaptive space-trees of general element classes in parallel.

  Copyright (C) 2025 the developers

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

#include <t8_cmesh.h>
#include <gtest/gtest.h>
#include <test/t8_gtest_macros.hxx>
#include <t8_cmesh/t8_cmesh_types.h>
#include <t8_cmesh/t8_cmesh_trees.h>
#include <t8_cmesh.hxx>
#include <t8_cmesh/t8_cmesh_examples.h>
#include <t8_cmesh/t8_cmesh_cad/t8_cmesh_boundary_node_list.hxx>
#include <test/t8_cmesh_generator/t8_cmesh_example_sets.hxx>
#include <unordered_set>

class t8_test_boundary_node_list: public testing::Test {
 protected:
  void
  SetUp () override
  {
    cmesh = t8_cmesh_hyperquad ();
  }
  void
  TearDown () override
  {
    t8_cmesh_unref (&cmesh);
  }

  t8_cmesh_t cmesh;
};

TEST_F (t8_test_boundary_node_list, some_random_ass_name)
{
  std::unordered_set<t8_gloidx_t> boundary_list = cmesh->boundary_node_list->get_boundary_node_list ();
  const int num_boundary_nodes = cmesh->boundary_node_list->get_boundary_node_list ().size ();

  EXPECT_TRUE (cmesh->compute_boundary_node_list);
  EXPECT_TRUE (num_boundary_nodes == 10);
  EXPECT_TRUE (t8_cmesh_is_committed (cmesh));
}
