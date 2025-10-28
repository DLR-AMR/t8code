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

/** \file t8_scheme_helpers.hxx
 * Provides helper functions for the easier implementation of schemes.
 */

#ifndef T8_SCHEME_HELPERS_HXX
#define T8_SCHEME_HELPERS_HXX

#include <t8_element.h>
#include <t8_types/t8_crtp.hxx>
#include <t8_eclass.h>
#include <t8_types/t8_vec.hxx>

/**
  * Class which provides helper functions and default implementations for different schemes.
  * Functions defined here are overridden by functions implemented in the schemes themselves.
  * \tparam TUnderlyingEclassScheme The scheme this helper class is adding functionality to.
  */
template <t8_eclass_t TEclass, class TUnderlyingEclassScheme>
class t8_scheme_helpers: public t8_crtp_basic<TUnderlyingEclassScheme> {
 protected:
  /**
   * Default constructor which is only accessible by derived classes.
   * This class cannot be constructed on its own.
   */
  t8_scheme_helpers () noexcept {};

 public:
  /** Return the tree dimension of this scheme.
   * \return The tree dimension of this scheme.
   */
  static constexpr size_t
  get_dimension (void) noexcept
  {
    return t8_eclass_to_dimension[TEclass];
  }

  /** Return the tree class of this scheme.
   * \return The tree class of this scheme.
   */
  static constexpr t8_eclass_t
  get_eclass (void) noexcept
  {
    return TEclass;
  }

 protected:
  /** Checks if the contents of a container have the same or a higher dimension than the element.
   * \tparam TCoordinateContainer The container type.
  */
  template <typename TCoordContainerType>
  static constexpr void
  assert_coord_container_dimensionality () noexcept
  {
    using coord_type = typename TCoordContainerType::value_type;
    assert_coord_dimensionality<coord_type> ();
  }

  /** Checks if a tuple has the same or a higher dimension than the element.
   * \tparam TCoord The tuple type.
  */
  template <typename TCoordType>
  static constexpr void
  assert_coord_dimensionality () noexcept
  {
    constexpr std::size_t coord_type_dim = t8_dim_v<TCoordType>;
    constexpr std::size_t elem_dim = t8_eclass_to_dimension[TEclass];
    static_assert (coord_type_dim >= elem_dim, "Input dimension for element_get_reference_coords is too small.");
  }
};

#endif /* T8_SCHEME_HELPERS_HXX */
