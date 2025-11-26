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
 * \file t8_crtp.hxx This file provides the basic CRTP pattern. There is also a more complex
 * operator pattern in \a t8_operators.hxx.
 * The pattern can be used for compile-time polymorphism of types.
 */

#ifndef T8_CRTP_HXX
#define T8_CRTP_HXX

/**
  * The basic CRTP pattern.
  * \tparam TUnderlying The Underlying type.
  *
  * \note There is also a CRTP pattern for operators \ref t8_crtp_operator in \ref t8_operators.hxx.
  */
template <typename TUnderlying>
struct t8_crtp_basic
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

#endif /* T8_CRTP_HXX */
