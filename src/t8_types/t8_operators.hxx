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
 * \tparam TUnderlying
 * \tparam crtpType 
 */
template <typename TUnderlying, template <typename> class crtpType>
struct t8_crtp_operator
{
  constexpr TUnderlying&
  underlying () noexcept
  {
    return static_cast<TUnderlying&> (*this);
  }

  constexpr const TUnderlying&
  underlying () const noexcept
  {
    return static_cast<const TUnderlying&> (*this);
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
 * \tparam TUnderlying
 */
template <typename TUnderlying>
struct Addable: t8_crtp_operator<TUnderlying, Addable>
{

  constexpr TUnderlying
  operator+ (const TUnderlying& other) const noexcept
  {
    return TUnderlying (this->underlying ().get () + other.get ());
  }
};

/**
 * \brief A template for subtractable types. Provides the - operator.
 * 
 * \tparam TUnderlying 
 */
template <typename TUnderlying>
struct Subtractable: t8_crtp_operator<TUnderlying, Subtractable>
{
  constexpr TUnderlying
  operator- (const TUnderlying& other) const noexcept
  {
    return TUnderlying (this->underlying ().get () - other.get ());
  }
};

/**
 * \brief A template for multipliable types. Provides the * operator.
 * 
 * \tparam TUnderlying 
 */
template <typename TUnderlying>
struct Multipliable: t8_crtp_operator<TUnderlying, Multipliable>
{
  constexpr TUnderlying
  operator* (const TUnderlying& other) const noexcept
  {
    return TUnderlying (this->underlying ().get () * other.get ());
  }
};

/**
 * \brief A template for dividable types. Provides the / operator.
 * 
 * \tparam TUnderlying 
 */
template <typename TUnderlying>
struct Dividable: t8_crtp_operator<TUnderlying, Dividable>
{
  constexpr TUnderlying
  operator/ (const TUnderlying& other) const noexcept
  {
    return TUnderlying (this->underlying ().get () / other.get ());
  }
};

/**
 * \brief A template for add-assignable types. Provides the += operator.
 * 
 * \tparam TUnderlying 
 */
template <typename TUnderlying>
struct AddAssignable: t8_crtp_operator<TUnderlying, AddAssignable>
{
  constexpr TUnderlying&
  operator+= (const TUnderlying& other) noexcept
  {
    this->underlying ().get () += other.get ();
    return this->underlying ();
  }
};

/**
 * \brief A template for incrementable types. Provides the ++ operator.
 * 
 * \tparam TUnderlying
 * 
 * \note The operator is a prefix operator.
 */
template <typename TUnderlying>
struct PrefixIncrementable: t8_crtp_operator<TUnderlying, PrefixIncrementable>
{
  TUnderlying&
  operator++ () noexcept
  {
    ++this->underlying ().get ();
    return this->underlying ();
  }
};

/**
 * \brief A template for decrementable types. Provides the -- operator.
 * 
 * \tparam TUnderlying 
 * 
 * \note The operator is a prefix operator.
 */
template <typename TUnderlying>
struct PrefixDecrementable: t8_crtp_operator<TUnderlying, PrefixDecrementable>
{
  TUnderlying&
  operator-- () noexcept
  {
    --this->underlying ().get ();
    return this->underlying ();
  }
};

template <typename TUnderlying>
struct Printable: t8_crtp_operator<TUnderlying, Printable>
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
 * \tparam TUnderlying 
 */
template <typename TUnderlying>
struct Swapable: t8_crtp_operator<TUnderlying, Swapable>
{
  constexpr void
  swap (TUnderlying& lhs, TUnderlying& other) noexcept
  {
    std::swap (lhs.get (), other.get ());
  }
};

/**
 * \brief A template for equality comparable types. Provides the == operator.
 * 
 * \tparam TUnderlying 
 */
template <typename TUnderlying>
struct EqualityComparable: t8_crtp_operator<TUnderlying, EqualityComparable>
{
  constexpr bool
  operator== (TUnderlying const& other) const noexcept
  {
    return this->underlying ().get () == other.get ();
  }
};

/**
 * \brief A template for << types. Provides the << operator.
 * 
 * \tparam T 
 */
template <typename TUnderlying, typename Parameter>
std::ostream&
operator<< (std::ostream& os, T8Type<TUnderlying, Parameter> const& p)
{
  p.print (os);
  return os;
}

/**
 * \brief A template for hashable types. Used to make a type hashable.
 * 
 * \tparam T 
 */
template <typename TUnderlying>
struct Hashable
{
  static constexpr bool is_hashable = true;
};

#endif /* T8_OPERATORS_HXX */
