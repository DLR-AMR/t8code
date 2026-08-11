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

/** \file t8_subelement_traits.hxx
 * Traits for subelement schemes.
 * Each trait defines the underlying scheme and subelement type associated with the subelement scheme.
 * This allows the common subelement implementation in \ref t8_subelement_scheme.hxx to remain generic 
 * while using the types required by each specialization.
 */

#pragma once

#include <t8_schemes/t8_standalone/t8_standalone_elements.hxx>
#include <t8_schemes/t8_standalone/t8_standalone_implementation.hxx>
#include <t8_schemes/t8_default/t8_default_tri/t8_dtri.h>
#include <t8_eclass/t8_eclass.h>
#include <t8_schemes/t8_subelement/t8_subelement_type.hxx>

/** Forward declaration of the quadrilateral subelement scheme. */
struct t8_subelementquad_scheme;

/** Forward declaration of the triangular subelement scheme. */
struct t8_subelementtri_scheme;

/** Traits associating a subelement scheme with its underlying scheme and subelement type.
 * \tparam TScheme The subelement scheme for which the traits are defined.
 */
template <typename TScheme>
struct t8_subelement_traits;

/** Traits specialization for quadrilateral subelements. */
template <>
struct t8_subelement_traits<t8_subelementquad_scheme>
{
  using UnderlyingScheme = t8_standalone_scheme<T8_ECLASS_QUAD>;
  using SubelementType = t8_subelement_element<t8_standalone_element<T8_ECLASS_QUAD>>;
};

/** Traits specialization for triangular subelements. */
template <>
struct t8_subelement_traits<t8_subelementtri_scheme>
{
  using UnderlyingScheme = t8_default_scheme_tri;
  using SubelementType = t8_subelement_element<t8_dtri>;
};
