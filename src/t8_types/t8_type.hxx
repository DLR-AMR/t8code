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

/**
 * \file This files gives a template for strong types in t8code.
 */

#ifndef T8_TYPE_HXX
#define T8_TYPE_HXX

#include <t8.h>

/**
 * \brief An implementation of strong type with additional competences.
 *
 * This class template allows the creation of a type that can be extended with
 * multiple competences. Each competence is a template class that takes the
 * main type as a template parameter.
 * 
 * This is heavily inspired by (and taken from) https://www.fluentcpp.com/2016/12/08/strong-types-for-strong-interfaces/
 *
 * \tparam T The type of the value to be stored.
 * \tparam Parameter An additional parameter for the type.
 * \tparam competence Variadic template parameter for the competences.
 */
template <typename T, typename Parameter, template <typename> class... competence>
class T8Type: public competence<T8Type<T, Parameter, competence...>>... {
 public:
  using value_type = T;

  explicit constexpr T8Type () = default;

  explicit constexpr T8Type (const T& value): value_ (value)
  {
  }

  /**
   * \brief Construct a new T8Type object
   * 
   * \tparam T_ref 
   * \param value 
   * 
   * \note This constructor is only enabled if T is not a reference.
   */
  template <typename T_ref = T>
  explicit constexpr T8Type (T&& value,
                             typename std::enable_if<!std::is_reference<T_ref> {}, std::nullptr_t>::type = nullptr)
    : value_ (std::move (value))
  {
  }

  constexpr T8Type&
  operator= (const T& value)
  {
    value_ = value;
    return *this;
  }

  constexpr T&
  get () noexcept
  {
    return value_;
  }

  constexpr T const&
  get () const noexcept
  {
    return std::move (value_);
  }

 private:
  T value_;
};

namespace std
{
/**
 * \brief Functor for hashing T8Type objects.
 *
 * This struct defines a functor that computes the hash value of a T8Type object.
 * It uses the std::hash function to generate the hash value based on the underlying
 * type T of the T8Type object.
 *
 * \tparam T The underlying type of the T8Type object.
 * \tparam Parameter Additional parameters for the T8Type object.
 * \tparam competence Variadic template parameters for additional competencies.
 *
 * \note This functor is enabled only if the T8Type object is hashable, as determined
 *       by the is_hashable member of the T8TypeImpl type.
 */
template <typename T, typename Parameter, template <typename> class... competence>
struct hash<T8Type<T, Parameter, competence...>>
{
  using T8TypeImpl = T8Type<T, Parameter, competence...>;
  using checkIfHashable = typename std::enable_if<T8TypeImpl::is_hashable, void>::type;

  size_t
  operator() (T8Type<T, Parameter, competence...> const& x) const
  {
    return std::hash<T> {}(x.get ());
  }
};
}  // namespace std

#endif /* T8_TYPE_HXX */
