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
 * \file t8_operators.hxx This file provides the CRTP pattern for operators.
 * The operators can be used by a \a T8Type to extend the functionality of the type.
 * There also is a more basic CRTP pattern in \ref t8_crtp.hxx.
 */

#ifndef T8_OPERATORS_HXX
#define T8_OPERATORS_HXX

#include <iostream>
#include <t8_types/t8_type.hxx>

/**
 * The CRTP pattern for operators.
 *
 * \tparam TUnderlying The Underlying type.
 * \tparam crtpType    The type to add via CRTP.
 *
 * \note There is also the basic CRTP pattern \ref t8_crtp without nested templates in \ref t8_crtp.hxx.
 *       Use this pattern only if the nested templates are needed.
 */
template <typename TUnderlying, template <typename> class crtpType>
struct t8_crtp_operator
{
  /**
   * Get the underlying type.
   */
  constexpr TUnderlying&
  underlying () noexcept
  {
    return static_cast<TUnderlying&> (*this);
  }

  /**
   * Get the underlying type.
   */
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
 *  A template for addable types. Provides the + operator.
 *
 * \tparam TUnderlying
 */
template <typename TUnderlying>
struct Addable: t8_crtp_operator<TUnderlying, Addable>
{
  /**
   * Add the value of \a other to the underlying type.
   *
   * \param [in] other The value to add to the underlying type.
   * \return The underlying type after the addition.
   */
  constexpr TUnderlying
  operator+ (const TUnderlying& other) const noexcept
  {
    return TUnderlying (this->underlying ().get () + other.get ());
  }
};

/**
 *  A template for subtractable types. Provides the - operator.
 *
 * \tparam TUnderlying
 */
template <typename TUnderlying>
struct Subtractable: t8_crtp_operator<TUnderlying, Subtractable>
{
  /**
   * Subtract the value of \a other from the underlying type.
   *
   * \param [in] other The value to subtract from the underlying type.
   * \return The underlying type after the subtraction.
   */
  constexpr TUnderlying
  operator- (const TUnderlying& other) const noexcept
  {
    return TUnderlying (this->underlying ().get () - other.get ());
  }
};

/**
 *  A template for multipliable types. Provides the * operator.
 *
 * \tparam TUnderlying
 */
template <typename TUnderlying>
struct Multipliable: t8_crtp_operator<TUnderlying, Multipliable>
{
  /**
   * Multiply the underlying type with \a other.
   *
   * \param [in] other The value to multiply the underlying type with.
   * \return The underlying type after the multiplication.
   */
  constexpr TUnderlying
  operator* (const TUnderlying& other) const noexcept
  {
    return TUnderlying (this->underlying ().get () * other.get ());
  }
};

/**
 *  A template for dividable types. Provides the / operator.
 *
 * \tparam TUnderlying
 */
template <typename TUnderlying>
struct Dividable: t8_crtp_operator<TUnderlying, Dividable>
{
  /**
   * Divide the underlying type by \a other.
   *
   * \param [in] other The value to divide the underlying type by.
   * \return The underlying type after the division.
   */
  constexpr TUnderlying
  operator/ (const TUnderlying& other) const noexcept
  {
    return TUnderlying (this->underlying ().get () / other.get ());
  }
};

/**
 *  A template for add-assignable types. Provides the += operator.
 *
 * \tparam TUnderlying
 */
template <typename TUnderlying>
struct AddAssignable: t8_crtp_operator<TUnderlying, AddAssignable>
{
  /**
   * Add-assign the value of \a other to the underlying type.
   *
   * \param [in] other The value to add to the underlying type.
   * \return The underlying type after the addition.
   */
  constexpr TUnderlying&
  operator+= (const TUnderlying& other) noexcept
  {
    this->underlying ().get () += other.get ();
    return this->underlying ();
  }
};

/**
 *  A template for incrementable types. Provides the ++ operator.
 *
 * \tparam TUnderlying
 *
 * \note The operator is a prefix operator.
 */
template <typename TUnderlying>
struct PrefixIncrementable: t8_crtp_operator<TUnderlying, PrefixIncrementable>
{
  /**
   * Increment the underlying type.
   *
   * \return The underlying type after the increment.
   */
  TUnderlying&
  operator++ () noexcept
  {
    ++this->underlying ().get ();
    return this->underlying ();
  }
};

/**
 *  A template for decrementable types. Provides the -- operator.
 *
 * \tparam TUnderlying
 *
 * \note The operator is a prefix operator.
 */
template <typename TUnderlying>
struct PrefixDecrementable: t8_crtp_operator<TUnderlying, PrefixDecrementable>
{
  /**
   * Decrement the underlying type.
   *
   * \return The underlying type after the decrement.
   */
  TUnderlying&
  operator-- () noexcept
  {
    --this->underlying ().get ();
    return this->underlying ();
  }
};

/**
 *  A template for printable types. Provides the << operator.
 *
 * \tparam TUnderlying
 */
template <typename TUnderlying>
struct Printable: t8_crtp_operator<TUnderlying, Printable>
{
  /**
   * Print the underlying type to the output stream.
   *
   * \param [in] os The output stream to print to.
   * \param [in] obj The object to print.
   * \return The output stream after printing.
   */
  friend std::ostream&
  operator<< (std::ostream& os, const TUnderlying& obj)
  {
    os << obj.get ();
    return os;
  }
};

/**
 *  A template for swapping types. Used to make a type swappable.
 *
 * \tparam TUnderlying
 */
template <typename TUnderlying>
struct Swapable: t8_crtp_operator<TUnderlying, Swapable>
{
  /**
   * Swap the underlying type with another underlying type.
   *
   * \param [in,out] lhs The left-hand side of the swap.
   * \param [in,out] other The right-hand side of the swap.
   */
  constexpr void
  swap (TUnderlying& lhs, TUnderlying& other) noexcept
  {
    std::swap (lhs.get (), other.get ());
  }
};

/**
 *  A template for equality comparable types. Provides the == operator.
 *
 * \tparam TUnderlying
 */
template <typename TUnderlying>
struct EqualityComparable: t8_crtp_operator<TUnderlying, EqualityComparable>
{
  /**
   * Check if the underlying types are equal.
   *
   * \param [in] lhs The left-hand side of the equality check.
   * \param [in] rhs The right-hand side of the equality check.
   * \return True if the underlying types are equal, false otherwise.
   */
  friend constexpr bool
  operator== (const TUnderlying& lhs, const TUnderlying& rhs) noexcept
  {
    return lhs.get () == rhs.get ();
  }

  /**
   * Check if the underlying type is not equal to another underlying type.
   *
   * \param [in] other The other underlying type to compare with.
   * \return True if the underlying types are not equal, false otherwise.
   */
  constexpr bool
  operator!= (TUnderlying const& other) const
  {
    return this->underlying ().get () != other.get ();
  }
};

/**
 *  A template for hashable types. Used to make a type hashable.
 *
 * \tparam TUnderlying
 */
template <typename TUnderlying>
struct Hashable
{
  /** Set if the underlying type is hashable. */
  static constexpr bool is_hashable = true;
};

/**
 *  A template for random accessible types. Provides the [] operator.
 *
 * \tparam TUnderlying
 */
template <typename TUnderlying>
struct RandomAccessible: t8_crtp_operator<TUnderlying, RandomAccessible>
{
  /**
   * Get the element at the given index.
   *
   * \param [in] index The index of the element to get.
   * \return The element at the given index.
   */
  auto
  operator[] (std::size_t index) -> decltype (auto)
  {
    return this->underlying ().get ()[index];
  }

  /**
   * Get the element at the given index.
   *
   * \param [in] index The index of the element to get.
   * \return The element at the given index.
   */
  auto
  operator[] (std::size_t index) const -> decltype (auto)
  {
    return this->underlying ().get ()[index];
  }

  /**
   * Get an iterator to the beginning of the underlying type.
   *
   * \return An iterator to the beginning of the underlying type.
   */
  auto
  begin () -> decltype (auto)
  {
    return this->underlying ().get ().begin ();
  }

  /**
   * Get an iterator to the begin of the underlying type.
   *
   * \return An iterator to the begin of the underlying type.
   */
  auto
  begin () const -> decltype (auto)
  {
    return this->underlying ().get ().begin ();
  }

  /**
   * Get an iterator to the end of the underlying type.
   *
   * \return An iterator to the end of the underlying type.
   */
  auto
  end () -> decltype (auto)
  {
    return this->underlying ().get ().end ();
  }

  /**
   * Get an iterator to the end of the underlying type.
   *
   * \return An iterator to the end of the underlying type.
   */
  auto
  end () const -> decltype (auto)
  {
    return this->underlying ().get ().end ();
  }

  /**
   * Get a pointer to the data of the underlying type.
   *
   * \return A pointer to the data of the underlying type.
   */
  auto
  data () -> decltype (auto)
  {
    return this->underlying ().get ().data ();
  }

  /**
   * Get a pointer to the data of the underlying type.
   *
   * \return A pointer to the data of the underlying type.
   */
  auto
  data () const -> decltype (auto)
  {
    return this->underlying ().get ().data ();
  }
};

#endif  // T8_OPERATORS_HXX
