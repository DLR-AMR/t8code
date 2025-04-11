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

#include <gtest/gtest.h>
#include <limits>

class DISABLED_t8_gtest_rank_times_global_num_elems_over_size: public testing::Test {
 protected:
  void
  SetUp () override
  {
    rank_growth = std::get<0> (GetParam ());
    global_num_elems_growth = std::get<1> (GetParam ());
  }
  void
  TearDown () override
  {
    /* nothing to do */
  }
  const uint32_t rank_growth;
  const uint32_t global_num_elems_growth;
  const int max_iter = 1000;
};

TEST_P (DISABLED_t8_gtest_rank_times_global_num_elems_over_size, large_numbers)
{
  /* we use exponential growth for all parts of the formular to get to large numbers quickly */
}

TEST_P (t8_gtest_rank_times_global_num_elems_over_size, small_numbers)
{
  for (uint32_t size_growth = rank_growth; size_growth < 10; ++size_growth) {
    uint64_t num_elems = 1;
    for (uint32_t ielem = 1; ielem < max_iter; ++ielem) {
      num_elems *= global_num_elems_growth;
      uint32_t size = 1;
      for (uint32_t isize = 1; isize < max_iter; ++isize) {
        size *= size_growth;
        uint32_t rank = 1;
        for (uint32_t irank = 1; irank < size && irank < max_iter; ++irank) {
          rank *= rank_growth;
          /* We only test for small numbers (much smaller that 2^64-1 here) */
          uint64_t check_result = rank * num_elems / size;
          /** As soon as this check is enabled call function to check with 
                     * rank, num_elems & size and compare the results. 
                     */
        }
      }
    }
  }
}

INSTANTIATE_TEST_SUITE_P (t8_gtest_rank_times_global_num_elems_over_size,
                          t8_gtest_rank_times_global_num_elems_over_size,
                          testing::Combine (testing::Values (1, 2, 3, 4, 5, 6, 7, 8, 9), testing
                                            : Values (1, 2, 3, 4, 5, 6, 7, 8, 9)));
