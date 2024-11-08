/*
  This file is part of t8code.
  t8code is a C library to manage a collection (a forest) of multiple
  connected adaptive space-trees of general element classes in parallel.

  Copyright (C) 2015 the developers

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

/** \file t8_default.hxx
 * This file is the point of entry for our default element implementation.
 *
 * This scheme points to a consistent implementation of all element classes.
 */

#pragma once

#include <t8_element.hxx>
#include <t8_schemes/t8_scheme.hxx>
#include <t8_schemes/t8_default/t8_default_vertex/t8_default_vertex.hxx>
#include <t8_schemes/t8_default/t8_default_line/t8_default_line.hxx>
#include <t8_schemes/t8_default/t8_default_quad/t8_default_quad.hxx>
#include <t8_schemes/t8_default/t8_default_hex/t8_default_hex.hxx>
#include <t8_schemes/t8_default/t8_default_tri/t8_default_tri.hxx>
#include <t8_schemes/t8_default/t8_default_tet/t8_default_tet.hxx>
#include <t8_schemes/t8_default/t8_default_prism/t8_default_prism.hxx>
#include <t8_schemes/t8_default/t8_default_pyramid/t8_default_pyramid.hxx>

T8_EXTERN_C_BEGIN ();

/** Return the default element implementation of t8code. */
t8_scheme *
t8_scheme_new_default (void);

/** Check whether a given eclass_scheme is one of the default schemes.
 * \param [in] scheme   A (pointer to a) scheme
 * \param [in] eclass   The eclass to check
 * \return              True (non-zero) if \a ts is one of the default schemes,
 *                      false (zero) otherwise.
 */
int
t8_eclass_scheme_is_default (const t8_scheme *scheme, const t8_eclass_t eclass);

T8_EXTERN_C_END ();
