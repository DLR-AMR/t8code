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

#ifndef T8_STANDALONE_HXX
#define T8_STANDALONE_HXX

#include <t8_schemes/t8_scheme.hxx>
#include <t8_schemes/t8_standalone/t8_standalone_implementation.hxx>

T8_EXTERN_C_BEGIN ();

/** Return the standalone element implementation of t8code. */
const t8_scheme *
t8_scheme_new_standalone (void);

/** Check whether a given eclass_scheme is one of the standalone schemes.
 * \param [in] scheme   A (pointer to a) scheme
 * \param [in] eclass   The eclass to check
 * \return              True (non-zero) if \a scheme is one of the default schemes,
 *                      false (zero) otherwise.
 */
int
t8_eclass_scheme_is_standalone (const t8_scheme *scheme, const t8_eclass_t eclass);

T8_EXTERN_C_END ();

#endif /* !T8_STANDALONE_HXX */
