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

#include <t8.h>
#include <gtest/gtest.h>

#include "test/t8_cmesh_generator/t8_cmesh_example_sets.hxx"
#include "test/t8_cmesh_generator/t8_gtest_cmesh_sum_of_sets.hxx"

TEST (t8_gtest_cmesh_iterator, begin_end_not_equal)
{
  EXPECT_TRUE (cmesh_list::cmesh_sums.begin () != cmesh_list::cmesh_sums.end ());
}

TEST (t8_gtest_cmesh_iterator, test_iteration)
{
  for (cmesh_sum_of_sets::Iterator iter = cmesh_list::cmesh_sums.begin (); iter != cmesh_list::cmesh_sums.end ();
       iter++) {
    std::string out;
    (*iter)->param_to_string (out);
    EXPECT_FALSE (out.empty ());
  }
}

TEST (t8_gtest_cmesh_iterator, test_iteration_with_stl)
{
  std::for_each (cmesh_list::cmesh_sums.begin (), cmesh_list::cmesh_sums.end (), [] (base_example *elem) {
    std::string out;
    elem->param_to_string (out);
    EXPECT_FALSE (out.empty ());
  });
}