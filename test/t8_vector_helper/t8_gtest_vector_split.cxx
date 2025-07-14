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
#include <t8_vector_helper/t8_vector_algorithms.hxx>
#include <numeric>
#include <t8.h>

// Integer division function, used as example for the category_func argument of vector_split.
constexpr unsigned int
split (const int value, const int div)
{
  T8_ASSERT (div != 0);
  T8_ASSERT ((value >= 0 && div > 0) || (value < 0 && div < 0));
  return value / div;
}

class test_vector_split: public testing::TestWithParam<int> {
 public:
  void
  SetUp () override
  {
    div = GetParam ();
    num_types = num_entries / div + 1;
    values.resize (num_entries);
    offsets.resize (num_types + 1);
    std::iota (values.begin (), values.end (), 0);
  }

  void
  TearDown () override
  {
  }

#if T8_TEST_LEVEL_INT == 2
  const size_t num_entries = 10;
#elif T8_TEST_LEVEL_INT == 1
  const size_t num_entries = 100;
#else
  const size_t num_entries = 100000;
#endif
  int div;
  size_t num_types;             // divisor for split(...) function, set via testing parameters
  std::vector<size_t> offsets;  // vector of category offsets
  std::vector<int> values;      // example vector of integers used to test vector_split
};

TEST_P (test_vector_split, test_split)
{
  vector_split (values.begin (), values.end (), offsets, std::function<size_t (int, int)> (split), div);

  EXPECT_EQ (offsets[0], 0);
  EXPECT_EQ (offsets.size (), num_types + 1);
  EXPECT_EQ (offsets[num_types], num_entries);

  for (size_t i = 0; i < num_types; ++i) {
    EXPECT_LE (offsets[i], offsets[i + 1]) << " at " << i;
  }
  for (size_t i = 0; i < num_types; ++i) {
    for (size_t j = offsets[i]; j < offsets[i + 1]; ++j) {
      EXPECT_EQ (values[j] / div, i) << " at " << j << " with value " << values[j] << ", num_types " << num_types
                                     << " num_entries " << num_entries;
    }
  }
}

INSTANTIATE_TEST_SUITE_P (test_vector_split, test_vector_split, testing::Values (2, 3, 4, 5, 6, 7, 8, 9, 10));
