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

/** \file t8_forest_standard_adapt.hxx
 * The adaptor for forests without holes and without subelements.
 */

#ifndef T8_FOREST_STANDARD_ADAPT_HXX
#define T8_FOREST_STANDARD_ADAPT_HXX

#include <t8_forest/t8_forest_adapt/t8_forest_adapt.hxx>

/** Adaptor for a source forest without holes and without subelements. It uses all hooks of \ref t8_adaptor.
 * \tparam TRecursive If true, adapt recursively, see \ref t8_adaptor.
 */
template <bool TRecursive>
struct t8_standard_adaptor: t8_adaptor<t8_standard_adaptor<TRecursive>, TRecursive>
{
  using t8_adaptor<t8_standard_adaptor, TRecursive>::t8_adaptor;
};

#endif /* T8_FOREST_STANDARD_ADAPT_HXX */
