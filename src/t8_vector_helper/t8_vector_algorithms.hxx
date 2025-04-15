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

/**
 * \file t8_vector_algorithms.hxx
 * 
 * Contains specialized algorithms on std::vector
 * 
 */

#ifndef T8_VECTOR_ALGORITHMS
#define T8_VECTOR_ALGORITHMS

#include <vector>
#include <algorithm>
#include <functional>
#include <t8.h>

/**
 * Compute the offsets of a categories of elements in a sorted vector.
 * T should be a type that can be compared with <.
 * 
 * /tparam T 
 * /param[in] vector            A vector holding elements of type T
 * /param[in, out] offsets      A vector. Will be resized to num_categories + 1. Will hold indices
 *                              j of \a vector that contain objects of category k, such that offsets[k] <0 j < offset[k+1]
 *                              If there are no elements of category k then offsets[k] = offsets[k +1]
 * /param[in] num_categories    The number of categories
 * /param[in] category_func     A function that takes an element of type T and returns the category of the element.
 * /param[in] data              A pointer to data that is passed to the category_func.
 */
template <typename T, typename... Args>
void
vector_split (const std::vector<T> &vector, std::vector<size_t> &offsets, const size_t num_categories,
              std::function<size_t (const T &, Args...)> &&category_func, Args... args)
{
  T8_ASSERT (std::is_sorted (vector.begin (), vector.end ()));
  const size_t count = vector.size ();
  /* Initialize everything with count, except for the first value. */
  offsets.resize (num_categories + 1);
  std::fill (offsets.begin (), offsets.end (), count);
  /* The first offset is set to zero */
  offsets[0] = 0;

  if (count == 0 || num_categories <= 1) {
    return;
  }

  /* We search between low and high for the next category */
  size_t low = 0;
  size_t high = count;
  size_t step = 1;
  while (step < num_categories) {
    size_t guess = low + (high - low) / 2;
    const size_t category = category_func (vector[guess], args...);
    T8_ASSERT (category < num_categories);
    if (category < step) {
      low = guess + 1;
    }
    else {
      std::fill (offsets.begin () + step, offsets.begin () + category + 1, guess);
      high = guess;
    }
    while (low == high) {
      ++step;
      high = offsets[step];
      if (step == num_categories) {
        return;
      }
    }
  }
}

#endif /* T8_VECTOR_ALGORITHMS */
