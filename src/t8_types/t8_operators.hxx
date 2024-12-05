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

#ifndef T8_OPERATORS_HXX
#define T8_OPERATORS_HXX

#include <iostream>
#include <t8_types/t8_type.hxx>

/**
 * \file This file provides the CRTP pattern for operators.
 * The operators can be used by a \a T8Type to extend the functionality of the type.
 */

/**
 * \brief The CRTP pattern for operators.
 * 
 * \tparam T 
 * \tparam crtpType 
 */
template <typename T, template <typename> class crtpType>
struct crtp
{
  T&
  underlying ()
  {
    return static_cast<T&> (*this);
  }

  T const&
  underlying () const
  {
    return static_cast<T const&> (*this);
  }
};

/*
 * The following is a list of competences that can be added to a type.
 * Each competence provides access to an operator of the underlying type. That way instead of
 * typing `my_int.get() + my_other_int.get()` you can type `my_int + my_other_int`. 
 */

/**
 * \brief A template for addable types. Provides the + operator.
 * 
 * \tparam T 
 */
template <typename T>
struct Addable: crtp<T, Addable>
{

  T
  operator+ (T const& other)
  {
    return T (this->underlying ().get () + other.get ());
  }
};

/**
 * \brief A template for subtractable types. Provides the - operator.
 * 
 * \tparam T 
 */
template <typename T>
struct Subtractable: crtp<T, Subtractable>
{
  T
  operator- (T const& other)
  {
    return T (this->underlying ().get () - other.get ());
  }
};

/**
 * \brief A template for multipliable types. Provides the * operator.
 * 
 * \tparam T 
 */
template <typename T>
struct Multipliable: crtp<T, Multipliable>
{
  T
  operator* (T const& other)
  {
    return T (this->underlying ().get () * other.get ());
  }
};

/**
 * \brief A template for dividable types. Provides the / operator.
 * 
 * \tparam T 
 */
template <typename T>
struct Dividable: crtp<T, Dividable>
{
  T
  operator/ (T const& other)
  {
    return T (this->underlying ().get () / other.get ());
  }
};

/**
 * \brief A template for add-assignable types. Provides the += operator.
 * 
 * \tparam T 
 */
template <typename T>
struct AddAssignable: crtp<T, AddAssignable>
{
  T&
  operator+= (T const& other)
  {
    this->underlying ().get () += other.get ();
    return this->underlying ();
  }
};

/**
 * \brief A template for incrementable types. Provides the ++ operator.
 * 
 * \tparam T 
 * 
 * \note The operator is a prefix operator.
 */
template <typename T>
struct PrefixIncrementable: crtp<T, PrefixIncrementable>
{
  T&
  operator++ ()
  {
    this->underlying ().get ()++;
    return this->underlying ();
  }
};

/**
 * \brief A template for decrementable types. Provides the -- operator.
 * 
 * \tparam T 
 * 
 * \note The operator is a prefix operator.
 */
template <typename T>
struct PrefixDecrementable: crtp<T, PrefixDecrementable>
{
  T&
  operator-- ()
  {
    this->underlying ().get ()--;
    return this->underlying ();
  }
};

template <typename T>
struct Printable: crtp<T, Printable>
{
  void
  print (std::ostream& os) const
  {
    os << this->underlying ().get ();
  }
};

/**
 * \brief A template for swapping types. Used to make a type swappable.
 * 
 * \tparam T 
 */
template <typename T>
struct Swapable: crtp<T, Swapable>
{
  void
  swap (T& other)
  {
    std::swap (this->underlying ().get (), other.get ());
  }
};

/**
 * \brief A template for equality comparable types. Provides the == operator.
 * 
 * \tparam T 
 */
template <typename T>
struct EqualityComparable: crtp<T, EqualityComparable>
{
  bool
  operator== (T const& other) const
  {
    return this->underlying ().get () == other.get ();
  }

  bool
  operator!= (T const& other) const
  {
    return !(*this == other);
  }
};

template <typename T>
struct Hashable
{
  static constexpr bool is_hashable = true;
};

/**
 * \brief A template for random accessible types. Provides the [] operator.
 * 
 * \tparam T 
 */
template <typename T>
struct RandomAccessible: crtp<T, RandomAccessible>
{
  auto
  operator[] (std::size_t index) -> decltype (auto)
  {
    return this->underlying ().get ()[index];
  }

  auto
  operator[] (std::size_t index) const -> decltype (auto)
  {
    return this->underlying ().get ()[index];
  }

  auto
  begin () -> decltype (auto)
  {
    return this->underlying ().get ().begin ();
  }

  auto
  begin () const -> decltype (auto)
  {
    return this->underlying ().get ().begin ();
  }

  auto
  end () -> decltype (auto)
  {
    return this->underlying ().get ().end ();
  }

  auto
  end () const -> decltype (auto)
  {
    return this->underlying ().get ().end ();
  }

  auto
  data () -> decltype (auto)
  {
    return this->underlying ().get ().data ();
  }

  auto
  data () const -> decltype (auto)
  {
    return this->underlying ().get ().data ();
  }
};

#endif  // T8_OPERATORS_HXX
