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
 * \file t8_unrolled_for.hxx
 *
 * Contains a compile time unrolled for implementation.
 *
 */

#ifndef T8_UNROLLED_FOR_HXX
#define T8_UNROLLED_FOR_HXX

#include <utility>

namespace t8_hidden_functions
{

/**
 * Index for an unrolled_for
 * \tparam TIndex The index in the for loop
 */
template <std::size_t TIndex>
struct unrolled_index
{
  static constexpr std::size_t value = TIndex;

  /**
   * Implicit conversion to size_t
   * \return The index
   */
  constexpr
  operator std::size_t () const noexcept
  {
    return TIndex;
  }
};

/** Implementation of a compile-time unrolled for. Has 0 runtime overhead compared to normal loops.
 *  Index can be used in compile-time contexts like array length or template arguments.
 * \tparam first        Starting variable, has to be an integer
 * \tparam last         Ending variable, has to be an integer (last iteration is \a iter < \a last)
 * \tparam TFunction    Body of the for loop as a lambda function.
 * \param [in] body     The body of the loop as a lambda function.
 *
 *  Usage:
 *  int x = 1
 *  unrolled_for<0, 3>([&](auto loop_count) {
 *      do_something<loop_count>(x);
 *  };
 *  expands to:
 *  do_something<0>(x);
 *  do_something<1>(x);
 *  do_something<2>(x);
 */
template <std::size_t first, std::size_t last, typename TFunction>
constexpr void
unrolled_for_impl (TFunction&& body)
{
  [&]<std::size_t... indexes> (std::index_sequence<indexes...>) {
    (..., body (t8_hidden_functions::unrolled_index<indexes + first> {}));
  }(std::make_index_sequence<last - first> {});
}

} /* namespace t8_hidden_functions */

/** Implementation of a compile-time unrolled for. Has 0 runtime overhead compared to normal loops.
 *  Index can be used in compile-time contexts like array length or template arguments.
 * \param[in]       FIRST   Starting variable, has to be an integer.
 * \param[in]       LAST    Ending variable, has to be an integer (last iteration is \a ITER < \a LAST).
 * \param[in, out]  ITER    Loop index.
 * \param[in]       BODY    The loop body.
 *
 *  Usage:
 *  int x = 1
 *  unrolled_for (0, 3, loop_index, {
        do_something<loop_index>(x);
    });

    expands to:

    do_something<0>(x);
    do_something<1>(x);
    do_something<2>(x);
 */
#define unrolled_for(FIRST, LAST, ITER, BODY) t8_hidden_functions::unrolled_for_impl<FIRST, LAST> ([&](auto ITER) BODY)

#endif /* T8_UNROLLED_FOR_HXX */
