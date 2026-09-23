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

/** \file t8_multilevel.hxx
 * Helper functions for the multilevel scheme.
 */

#ifndef T8_MULTILEVEL_HXX
#define T8_MULTILEVEL_HXX

#include <t8_schemes/t8_scheme.hxx>
#include <t8_schemes/t8_multilevel/t8_multilevel_implementation.hxx>

/** Check whether a given eclass_scheme is one of the standalone schemes.
 * \tparam TUnderlyingEclassScheme  The underlying eclass scheme of the multilevel.
 * \tparam TUnderlyingElementType   The underlying element type of the multilevel scheme.
 * \param [in] scheme   A (pointer to a) scheme
 * \param [in] eclass   The eclass to check
 * \return              True if scheme is a multilevel scheme,
 *                      false otherwise.
 */
template <class TUnderlyingEclassScheme, typename TUnderlyingElementType>
bool
t8_eclass_scheme_is_multilevel (const t8_scheme *scheme, const t8_eclass_t eclass)
{
  return scheme->check_eclass_scheme_type<t8_multilevel_scheme<TUnderlyingEclassScheme>, TUnderlyingElementType> (
    eclass);
}

#endif /* !T8_MULTILEVEL_HXX */
