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

size_t
split (const int &value, const size_t &num_types)
{
  return (size_t) value / (num_types - 1);
}

class test_vector_split: public testing::TestWithParam<int> {
 public:
  void
  SetUp () override
  {
    num_types = GetParam ();
    values.resize (num_entries);
    std::iota (values.begin (), values.end (), 0);
  }

  void
  TearDown () override
  {
  }

#if T8CODE_TEST_LEVEL == 2
  const size_t num_entries = 10;
#elif T8CODE_TEST_LEVEL == 1
  const size_t num_entries = 100;
#else
  const size_t num_entries = 1000;
#endif
  size_t num_types;
  std::vector<size_t> offsets;
  std::vector<int> values;
};

TEST_P (test_vector_split, test_split)
{
  vector_split<int, const size_t &> (values, offsets, num_types,
                                     std::function<size_t (const int &, const size_t &)> (split), num_types);
  EXPECT_EQ (offsets[0], 0);
  EXPECT_EQ (offsets.size (), num_types + 1);
  for (size_t i = 0; i < num_types; ++i) {
    EXPECT_LE (offsets[i], offsets[i + 1]);
  }
  for (size_t i = 0; i < num_types; ++i) {
    for (size_t j = offsets[i]; j < offsets[i + 1]; ++j) {
      EXPECT_EQ (values[j] / (num_types - 1), i);
    }
  }
}

INSTANTIATE_TEST_SUITE_P (test_vector_split, test_vector_split, testing::Values (2, 3, 4, 5, 6, 7, 8, 9, 10));
