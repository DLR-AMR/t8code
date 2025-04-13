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
#include <cmath>
#include <t8_cmesh.hxx>

class t8_gtest_rank_times_global_num_elems_over_size: public testing::TestWithParam<std::tuple<int, int, int>> {
 protected:
  void
  SetUp () override
  {
    rank_growth = std::get<0> (GetParam ());
    elem_growth = std::get<1> (GetParam ());
    size_growth = std::get<2> (GetParam ());
    rank_iter = log (std::numeric_limits<uint32_t>::max ()) / log (rank_growth);
    elem_iter = log (std::numeric_limits<uint64_t>::max ()) / log (elem_growth);
    size_iter = log (std::numeric_limits<uint32_t>::max ()) / log (size_growth);
  }
  void
  TearDown () override
  {
    /* nothing to do */
  }
  uint32_t rank_growth;
  uint32_t elem_growth;
  uint32_t size_growth;
#if T8CODE_TEST_LEVEL == 0
  const uint32_t max_iter = 100;
#elif T8CODE_TEST_LEVEL == 1
  const uint32_t max_iter = 50;
#else
  const uint32_t max_iter = 10;
#endif
  uint32_t rank_iter;
  uint32_t elem_iter;
  uint32_t size_iter;
};

TEST_P (t8_gtest_rank_times_global_num_elems_over_size, large_numbers)
{
  /** 
   * We test the formula rank * num_elems / size for large numbers.
   * The formula is used in the t8code library to compute the number of elements
   * that are owned by a rank.
   * 
   * We use a recursive approach to compute the result.
   * The result in the innermost loop is computed by:
   * check_result(n, m) = num_elems^n * rank^m / size
   * 
   * To prevent overflow we use the following approach:
   * floor (A^n*B/C) = A * floor (A^(n-1)*B/C) + floor(A/C)*(A^(n-1)*B % C) 
   * and + (A % C)*((A^(n-1)B) % C) / C is used to compute the remainder of the next step. 
   * for the inner and outer loop. 
   * We use integer division, therefore we store the remainder of each update to 
   * prevent rounding errors.
  */
  uint64_t size = 1;
  for (uint32_t isize = 1; isize < size_iter; ++isize) {
    /* The very first result is 1 * 1 / size */
    uint64_t check_result_elem = 1 / size;
    uint64_t check_result_elem_remain = 1;

    uint64_t num_elems = 1;
    /* Initialize factors */
    for (uint32_t ielem = 1; ielem < elem_iter; ++ielem) {
      uint32_t rank = 1;

      /** Used to compute elem^n * rank^m / size, where n is fixed. */
      uint64_t check_result = check_result_elem;
      uint64_t rank_remainder = check_result_elem_remain;
      for (uint32_t irank = 1; irank < rank_iter && rank <= size; ++irank) {
        const uint64_t computed_result = t8_cmesh_get_first_element_of_process (rank, size, num_elems);

        check_result = (rank == size) ? num_elems : check_result;
        ASSERT_EQ (computed_result, check_result)
          << "rank: " << rank << " num_elems: " << num_elems << " size: " << size;

        /* Update the result with respect to the updated rank */
        check_result *= rank_growth;
        check_result += rank_growth * rank_remainder / size;
        rank_remainder = (rank_growth * rank_remainder) % size;

        rank *= rank_growth;
      }
      /* Update the result with respect to the updated number of elements. */
      check_result_elem *= elem_growth;
      check_result_elem += elem_growth * check_result_elem_remain / size;
      check_result_elem_remain = (elem_growth * check_result_elem_remain) % size;

      num_elems *= elem_growth;
    }
    size *= size_growth;
  }
}

TEST_P (t8_gtest_rank_times_global_num_elems_over_size, small_numbers)
{
  uint64_t num_elems = 1;
  for (uint32_t ielem = 1; ielem < max_iter; ++ielem) {
    num_elems += elem_growth;
    uint32_t size = 1;
    for (uint32_t isize = 1; isize < max_iter; ++isize) {
      size += size_growth;
      uint32_t rank = 1;
      for (uint32_t irank = 1; irank * rank_growth < size && irank < max_iter; ++irank) {
        rank += rank_growth;
        /* We only test for small numbers (much smaller that 2^64-1 here) */
        const uint64_t check_result = rank * num_elems / size;
        const uint64_t computed_result = t8_cmesh_get_first_element_of_process (rank, size, num_elems);
        EXPECT_EQ (check_result, computed_result)
          << "rank: " << rank << " num_elems: " << num_elems << " size: " << size;
      }
    }
  }
}

INSTANTIATE_TEST_SUITE_P (t8_gtest_rank_times_global_num_elems_over_size,
                          t8_gtest_rank_times_global_num_elems_over_size,
                          testing::Combine (testing::Range (1, 10), testing::Range (1, 10), testing::Range (1, 10)));
