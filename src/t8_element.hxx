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

/** \file t8_element.hxx
 * This file defines basic operations on an element in a refinement tree.
 *
 * All operations work for all element classes by providing a inline function table.
 * For each element class, one implementation of the type and inline table is required.
 */

#ifndef T8_ELEMENT_HXX
#define T8_ELEMENT_HXX

#include <sc_refcount.h>
#include <t8_eclass.h>
#include <t8_element.h>
#include <t8_schemes/t8_crtp.hxx>

/** This class holds functions for a particular element class. */
template <class TUnderlyingEclassScheme>
class t8_eclass_scheme: public t8_crtp<TUnderlyingEclassScheme> {
 private:
  /** Private constructor which can only be used by derived schemes.
   * \param [in] tree_class The tree class of this element scheme.
   * \param [in] elem_size  The size of the elements this scheme holds.
  */
  t8_eclass_scheme (const t8_eclass_t tree_class, const size_t elem_size, void *ts_context)
    : element_size (elem_size), eclass (tree_class), ts_context (ts_context) {};
  friend TUnderlyingEclassScheme;

 protected:
  const size_t element_size; /**< The size in bytes of an element of class \a eclass */
  void *ts_context;          /**< Anonymous implementation context. */

 public:
  const t8_eclass_t eclass; /**< The tree class */

  /** Return the size of any element of a given class.
   * \return                      The size of an element of class \b ts.
   * We provide a default implementation of this routine that should suffice
   * for most use cases.
   */
  inline size_t
  get_element_size (void) const
  {
    return element_size;
  }
};

#endif /* !T8_ELEMENT_HXX */
