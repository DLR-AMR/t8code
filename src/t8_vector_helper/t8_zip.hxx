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
 * \file t8_zip.hxx
 *
 * Implements std zip for C++20 and injects it into the std namespace.
 * This can just be deleted when upgrading to C++23.
 *
 * This file is a copied and modified version from https://github.com/alemuntoni/zip-views which is licensed under MIT:
 *
 * Copyright 2020 CommitThis Ltd
 * Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"),
 * to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense,
 * and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:
 * The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM,
 * DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE
 * OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
 */

#ifndef ZIP_VIEW_HPP
#define ZIP_VIEW_HPP

#include <cassert>
#include <functional>
#include <iomanip>
#include <iostream>
#include <list>
#include <string>
#include <vector>
#include <typeinfo>
#include <ranges>

namespace c9::detail
{

template <typename... Args, std::size_t... Index>
auto
any_match_impl (std::tuple<Args...> const &lhs, std::tuple<Args...> const &rhs, std::index_sequence<Index...>) -> bool
{
  auto result = false;
  result = (... || (std::get<Index> (lhs) == std::get<Index> (rhs)));
  return result;
}

template <typename... Args>
auto
any_match (std::tuple<Args...> const &lhs, std::tuple<Args...> const &rhs) -> bool
{
  return any_match_impl (lhs, rhs, std::index_sequence_for<Args...> {});
}

template <std::ranges::viewable_range... Rng>
class zip_iterator {
 public:
  using value_type = std::tuple<std::ranges::range_reference_t<Rng>...>;

  zip_iterator () = delete;

  zip_iterator (std::ranges::iterator_t<Rng> &&...iters)
    : m_iters { std::forward<std::ranges::iterator_t<Rng>> (iters)... }
  {
  }

  auto
  operator++ () -> zip_iterator &
  {
    std::apply ([] (auto &&...args) { ((++args), ...); }, m_iters);
    return *this;
  }

  auto
  operator++ (int) -> zip_iterator
  {
    auto tmp = *this;
    ++*this;
    return tmp;
  }

  auto
  operator!= (zip_iterator const &other) const
  {
    return !(*this == other);
  }

  auto
  operator== (zip_iterator const &other) const
  {
    auto result = false;
    return any_match (m_iters, other.m_iters);
  }

  auto
  operator* () -> value_type
  {
    return std::apply ([] (auto &&...args) { return value_type (*args...); }, m_iters);
  }

 private:
  std::tuple<std::ranges::iterator_t<Rng>...> m_iters;
};

template <std::ranges::viewable_range... T>
class zipper {
 public:
  using zip_type = zip_iterator<T...>;

  template <typename... Args>
  zipper (Args &&...args): m_args { std::forward<Args> (args)... }
  {
  }

  auto
  begin () -> zip_type
  {
    return std::apply ([] (auto &&...args) { return zip_type (std::ranges::begin (args)...); }, m_args);
  }
  auto
  end () -> zip_type
  {
    return std::apply ([] (auto &&...args) { return zip_type (std::ranges::end (args)...); }, m_args);
  }

 private:
  std::tuple<T...> m_args;
};

}  // namespace c9::detail

namespace std::ranges::views
{

template <std::ranges::viewable_range... T>
auto
zip (T &&...t)
{
  return c9::detail::zipper<T...> { std::forward<T> (t)... };
}

}  // namespace std::ranges::views

#endif /* T8_ZIP */
