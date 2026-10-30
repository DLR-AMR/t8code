/*
  This file is part of t8code.
  t8code is a C library to manage a collection (a forest) of multiple
  connected adaptive space-trees of general element classes in parallel.

  Copyright (C) 2026 the developers

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

/** \file t8_element_array.h
 * This file defines the C interface of the t8_element_array_view. For more information
 * refer to the C++ interface in \ref t8_element_array.hxx.
 */

#ifndef T8_ELEMENT_ARRAY_H
#define T8_ELEMENT_ARRAY_H

#include <t8.h>
#include <t8_element/t8_element.h>
#include <t8_schemes/t8_scheme.h>

/** A non-owning, read-only view of contiguous elements of one eclass scheme.
 *  Opaque pointer for C interface.
 *  Detailed documentation at \ref t8_element_array_view.
 */
typedef struct t8_element_array_view t8_element_array_view_c;

T8_EXTERN_C_BEGIN ();

/** Return the number of elements in a view.
 * \param [in] view The view.
 * \return          The number of elements in \a view.
 */
size_t
t8_element_array_view_get_count (const t8_element_array_view_c *view);

/** Return an element of a view.
 * \param [in] view  The view.
 * \param [in] index The index of the element. Must be smaller than the number of elements in \a view.
 * \return           A pointer to the element at position \a index in \a view.
 */
const t8_element_t *
t8_element_array_view_get_element (const t8_element_array_view_c *view, const size_t index);

/** Return the scheme of the elements in a view.
 * \param [in] view The view.
 * \return          The scheme of the elements in \a view.
 */
const t8_scheme_c *
t8_element_array_view_get_scheme (const t8_element_array_view_c *view);

/** Return the tree class of the elements in a view.
 * \param [in] view The view.
 * \return          The tree class of the elements in \a view.
 */
t8_eclass_t
t8_element_array_view_get_tree_class (const t8_element_array_view_c *view);

T8_EXTERN_C_END ();

#endif /* !T8_ELEMENT_ARRAY_H */
