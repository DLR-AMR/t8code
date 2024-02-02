/*
This file is part of t8code.
t8code is a C library to manage a collection (a forest) of multiple
connected adaptive space-trees of general element classes in parallel.

Copyright (C) 2023 the developers

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

#ifndef T8_GTEST_CMESH_CREATOR_BASE_HXX
#define T8_GTEST_CMESH_CREATOR_BASE_HXX

#include <t8_cmesh.h>
#include <tuple>
#include <functional>
#include <iterator>
#include <vector>

class cart_prod_base {
 public:
  virtual bool
  next ()
    = 0;

  virtual t8_cmesh_t
  gen_cmesh ()
    = 0;

  virtual void
  set_to_first ()
    = 0;

  virtual void
  set_to_end ()
    = 0;

  virtual bool
  operator< (const cart_prod_base& other)
    = 0;

  size_t index = 0;
};

template <typename Args>
auto
vector_to_iter_pair (std::vector<Args> vec)
{
  return std::make_pair (vec.begin (), vec.end ());
}

template <typename T, typename B>
bool
increment (const B& begins, std::pair<T, T>& r)
{
  ++r.first;
  if (r.first == r.second)
    return true;
  return false;
}
template <typename T, typename... TT, typename B>
bool
increment (const B& begins, std::pair<T, T>& r, std::pair<TT, TT>&... rr)
{
  ++r.first;
  if (r.first == r.second) {
    r.first = std::get<std::tuple_size<B>::value - sizeof...(rr) - 1> (begins);
    return increment (begins, rr...);
  }
  return false;
}

template <typename OutputIterator, typename... Iter>
void
cartesian_product (OutputIterator out, std::pair<Iter, Iter>... ranges)
{
  const auto begins = std::make_tuple (ranges.first...);
  while (increment (begins, ranges...)) {
    out = { *ranges.first... };
  }
}

/**
 * A base class for cmesh_creators. 
 * 
 */
template <class... Iter>
class cmesh_args_cart_prod: cart_prod_base {
 public:
  cmesh_args_cart_prod (std::pair<Iter, Iter>... ranges,
                        std::function<t8_cmesh_t (typename Iter::value_type...)> cmesh_function)
    : cmesh_example (cmesh_function)
  {
    cartesian_product (std::back_inserter (cart_prod), ranges...);
  }

  t8_cmesh_t
  gen_cmesh ()
  {
    return std::apply (cmesh_example, cart_prod[index]);
  }

  virtual void
  set_to_first ()
  {
    index = 0;
  }

  virtual void
  set_to_end ()
  {
    index = cart_prod.size ();
  }

  virtual bool
  next ()
  {
    ++index;
    if (index >= cart_prod.size ()) {
      return false;
    }
    else {
      return true;
    }
  }

  virtual bool
  operator< (const cart_prod_base& other)
  {
    return index < other.index;
  }

  std::function<t8_cmesh_t (typename Iter::value_type...)> cmesh_example;
  std::vector<std::tuple<typename Iter::value_type...>> cart_prod;
};

#endif /* T8_GTEST_CMESH_CREATOR_BASE_HXX */