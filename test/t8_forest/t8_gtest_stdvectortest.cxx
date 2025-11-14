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
#include <vector>
#include <typeinfo>
#include <t8_data/t8_stdvector_conversion.hxx>  // header for sc_array conversion
#include "sc_containers.h"                      // header for sc_array_t and related functions

// Define sc_array_count function
size_t
sc_array_count (const sc_array_t* sc_arr)
{
  return sc_arr->elem_count;  // Access the count member
}

template <typename T>
class VectorTest: public ::testing::Test {
 protected:
  void
  SetUp () override
  {
    // Initialize the vector with some test values
    vec = { T (1), T (2), T (3), T (4), T (5) };
    arr = vec.data ();
    // Convert std::vector to sc_array
    sc_arr = t8_create_sc_array_view_from_vector (vec);
  }

  void
  TearDown () override
  {
    // Clean up the sc_array to prevent memory leaks
    sc_array_destroy (sc_arr);
  }

  std::vector<T> vec;  // Standard vector
  const T* arr;        // Raw pointer to the vector data
  sc_array_t* sc_arr;  // Converted sc_array
};

TYPED_TEST_SUITE_P (VectorTest);

TYPED_TEST_P (VectorTest, LengthTest)
{
  // Test that the vector length is as expected
  ASSERT_EQ (this->vec.size (), 5u) << "Vector length should be 5 for testing purposes";
  // Test that the sc_array length matches the vector length
  ASSERT_EQ (sc_array_count (this->sc_arr), this->vec.size ()) << "sc_array length should match vector length";
}

TYPED_TEST_P (VectorTest, ElementComparisonTest)
{
  for (size_t i = 0; i < this->vec.size (); ++i) {
    // Compare elements in the vector
    EXPECT_EQ (this->vec[i], this->arr[i]) << "Element mismatch at index " << i;
    // Compare elements between the vector and sc_array
    EXPECT_EQ (this->vec[i], *(reinterpret_cast<const TypeParam*> (sc_array_index (this->sc_arr, i))))
      << "Element mismatch between vector and sc_array at index " << i;
  }
}

REGISTER_TYPED_TEST_SUITE_P (VectorTest, LengthTest, ElementComparisonTest);

using MyTypes = ::testing::Types<int, double, float, char, short, long, long long, unsigned int, unsigned char,
                                 unsigned short, unsigned long, unsigned long long>;

INSTANTIATE_TYPED_TEST_SUITE_P (MyVectorTests, VectorTest, MyTypes, );