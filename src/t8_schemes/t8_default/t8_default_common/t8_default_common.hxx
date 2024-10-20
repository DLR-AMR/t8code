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

/** \file t8_default_common.hxx
 * We provide some functions that are useful across element classes.
 */

#pragma once

#include <t8_element.hxx>
#include <t8_schemes/t8_crtp.hxx>

/* Macro to check whether a pointer (VAR) to a base class, comes from an
 * implementation of a child class (TYPE). */
#define T8_COMMON_IS_TYPE(VAR, TYPE) ((dynamic_cast<TYPE> (VAR)) != NULL)

class t8_default_scheme_common_c {
 public:
  /** Destructor for all default schemes */
  ~t8_default_scheme_common_c ();

  /** Compute the number of corners of a given element. */
  int
  element_get_num_corners (const t8_element_t *elem) const;

  /** Allocate space for a bunch of elements. */
  void
  element_new (int length, t8_element_t **elem) const;

  /** Deallocate space for a bunch of elements. */
  void
  element_destroy (int length, t8_element_t **elem) const;

  /** Return the shape of an element */
  t8_element_shape_t
  element_get_shape (const t8_element_t *elem) const;

  /** Count how many leaf descendants of a given uniform level an element would produce.
   * \param [in] t     The element to be checked.
   * \param [in] level A refinement level.
   * \return Suppose \a t is uniformly refined up to level \a level. The return value
   * is the resulting number of elements (of the given level).
   * Each default element (except pyramids) refines into 2^{dim * (level - level(t))}
   * children.
   */
  t8_gloidx_t
  element_count_leaves (const t8_element_t *t, int level) const;

  /** Compute the number of siblings of an element. That is the number of 
   * Children of its parent.
   * \param [in] elem The element.
   * \return          The number of siblings of \a element.
   * Note that this number is >= 1, since we count the element itself as a sibling.
   */
  int
  element_get_num_siblings (const t8_element_t *elem) const;

  /** Count how many leaf descendants of a given uniform level the root element will produce.
   * \param [in] level A refinement level.
   * \return The value of \ref t8_element_count_leaves if the input element
   *      is the root (level 0) element.
   */
  t8_gloidx_t
  count_leaves_from_root (int level) const;

#if T8_ENABLE_DEBUG
  void
  element_debug_print (const t8_element_t *elem) const;
#endif
};
