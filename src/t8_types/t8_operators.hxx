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

template <typename T>
struct Addable: crtp<T, Addable>
{

  T
  operator+ (T const& other)
  {
    return T (this->underlying ().get () + other.get ());
  }
};

template <typename T>
struct Subtractable: crtp<T, Subtractable>
{
  T
  operator- (T const& other)
  {
    return T (this->underlying ().get () - other.get ());
  }
};

template <typename T>
struct Multipliable: crtp<T, Multipliable>
{
  T
  operator* (T const& other)
  {
    return T (this->underlying ().get () * other.get ());
  }
};

template <typename T>
struct Dividable: crtp<T, Dividable>
{
  T
  operator/ (T const& other)
  {
    return T (this->underlying ().get () / other.get ());
  }
};

template <typename T>
struct Incrementable: crtp<T, Incrementable>
{
  T&
  operator++ ()
  {
    this->underlying ().get ()++;
    return this->underlying ();
  }
};

template <typename T>
struct Decrementable: crtp<T, Decrementable>
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

template <typename T>
struct Swapable: crtp<T, Swapable>
{
  void
  swap (T& other)
  {
    std::swap (this->underlying ().get (), other.get ());
  }
};

template <typename T>
struct EqualityComparable: crtp<T, EqualityComparable>
{
  bool
  operator== (T const& other) const
  {
    return this->underlying ().get () == other.get ();
  }
};

template <typename T, typename Parameter>
std::ostream&
operator<< (std::ostream& os, T8Type<T, Parameter> const& p)
{
  p.print (os);
  return os;
}

#endif /* T8_OPERATORS_HXX */
