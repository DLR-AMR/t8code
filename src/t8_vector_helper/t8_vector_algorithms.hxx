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

#include <algorithm>
#include <functional>
#include <numeric>
#include <iterator>
#include <t8.h>

/**
 * Compute the offsets of a categories of elements in a sorted iterable range given by \a begin and \a end.
 * The value type of the iterator should be comparable with <.
 * This is a re-implementation of sc_array_split
 * 
 * /tparam TIterator            An input iterator type
 * /tparam TSentinel            A sentinel type for the iterator
 * /tparam TContainer           A container type that holds the offsets. It should be a contiguous container
 * /tparam TCategory            The type of the category. It should be an integral and totally ordered type
 * /tparam Args                 The type of the arguments passed to the category_func
 *
 * /param[in] begin             An iterator pointing to the first element of the range
 * /param[in] end               An iterator pointing to the last element of the range
 * /param[in, out] offsets      A Container holding num_categories + 1 elements. Will hold indices
 *                              j of the range \a begin and \a end that contain objects of category k, such that offsets[k] <0 j < offset[k+1]
 *                              If there are no elements of category k then offsets[k] = offsets[k +1]
 * /param[in] category_func     A function that takes an element of the value type of the iterators \a begin / \a end and
 *                              returns the category of the element. A category should be in [0, num_categories) (However TCategory does not have to be unsigned)
 * /param[in] args              A parameter pack of arguments passed to the category_func
 */
template <std::input_iterator TIterator, std::sentinel_for<TIterator> TSentinel,
          std::ranges::contiguous_range TContainer, std::unsigned_integral TCategory, typename... Args>
constexpr void
vector_split (
  const TIterator begin, const TSentinel end, TContainer &offsets,
  std::function<TCategory (const typename std::iterator_traits<TIterator>::value_type, Args...)> &&category_func,
  Args... args)
{
  T8_ASSERT (std::is_sorted (begin, end));
  T8_ASSERT (begin != end);
  T8_ASSERT (!offsets.empty ());
  const TCategory count = std::distance (begin, end);
  const TCategory num_categories = offsets.size () - 1;
  /* Initialize everything with count, except for the first value. */
  std::fill (std::begin (offsets), std::end (offsets), count);
  /* The first offset is set to zero */
  offsets[0] = 0;

  if (count == 0 || num_categories <= 1) {
    return;
  }

  /* We search between low and high for the next category */
  TCategory low = 0;
  TCategory high = count;
  TCategory step = 1;
  while (step < num_categories) {
    // Using binary search to find the next category boundary
    const TCategory guess = std::midpoint (low, high);
    const TCategory category = category_func (*(begin + guess), args...);
    T8_ASSERT (0 <= category && category < num_categories);

    if (category < step) {
      // If the category is smaller than the current step, adjust low
      low = guess + 1;
    }
    else {
      // Fill offsets for all categories between step and category
      std::fill (std::begin (offsets) + step, std::begin (offsets) + category + 1, guess);
      // Minimize the high value to narrow the search space
      high = guess;
    }

    // Advance step and update high when low equals high
    while (low == high) {
      ++step;
      if (step == num_categories) {
        return;
      }
      high = offsets[step];
    }
  }
}

#endif /* T8_VECTOR_ALGORITHMS */
