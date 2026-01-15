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

  /** Given a face of an element and a level coarser than (or equal to)
   * the element's level, return the face number
   * of the ancestor of the element that matches the element's face. Or return -1 if
   * no face of the ancestor matches the face.
   * \param [in]  element    The element.
   * \param [in]  ancestor_level A refinement level smaller than (or equal to) \a element's level.
   * \param [in]  face    Then number of a face of \a element.
   * \return              If \a face of \a element is a subface of a face of \a element's ancestor at level \a ancestor_level,
   *                      the face number of this face. Otherwise -1.
   * \note For the root element this function always returns \a face.
   */
  inline int
  element_face_get_ancestor_face (const t8_element_t *element, const int ancestor_level, const int face) const
  {
    auto underlying_impl = this->underlying ();  // Reference to the underlying scheme implementation

    const int element_level = underlying_impl.element_get_level (element);
    T8_ASSERT (element_level >= ancestor_level);
    if (element_level == ancestor_level) {
      // On the same level, the return value is the face itself
      return face;
    }
    // Allocate memory for a temporary element.
    t8_element_t *parent;
    underlying_impl.element_new (1, &parent);
    // Pointer to a temoporary element, that will move up the refinement hierarchy
    const t8_element_t *temp_element = element;
    int temp_face = face;
    for (int ilevel = element_level; ilevel > ancestor_level; --ilevel) {
      // Go one level up in the refinement hierarchy with the face
      temp_face = underlying_impl.element_face_get_parent_face (temp_element, temp_face);
      if (temp_face == -1) {
        // This face is not a subface of an ancestor face.
        underlying_impl.element_destroy (1, &parent);
        return -1;
      }
      underlying_impl.element_get_parent (temp_element, parent);
      temp_element = parent;
    }
    underlying_impl.element_destroy (1, &parent);
    return temp_face;
  }
};

#endif /* T8_SCHEME_HELPERS_HXX */
