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
/** \file t8_forest_adapt.cxx
 * Implements functions declared in \ref t8_forest_adapt.h.
 */

#include <t8_forest/t8_forest_adapt.h>
#include <t8_forest/t8_forest_types.h>
#include <t8_forest/t8_forest_general.h>
#include <t8_forest/t8_forest_subelement.hxx>
#include <t8_schemes/t8_scheme.hxx>
#include <t8_schemes/t8_subelement/t8_subelement.hxx>
#include <t8_forest/t8_forest_adapt/t8_forest_standard_adapt.hxx>
#include <t8_forest/t8_forest_adapt/t8_forest_subelement_adapt.hxx>

/** Adapt \a forest with the serial or the recursive variant of \a TAdaptor, as set in \a forest.
 * \tparam TAdaptor The adaptor, templated on whether it adapts recursively.
 * \param [in, out] forest The forest under construction.
 */
template <template <bool> typename TAdaptor>
static void
t8_forest_adapt_with (t8_forest_t forest)
{
  if (forest->set_adapt_recursive) {
    TAdaptor<true> adaptor (forest);
    adaptor.adapt ();
  }
  else {
    TAdaptor<false> adaptor (forest);
    adaptor.adapt ();
  }
}

/* We want to export the whole implementation to be callable from "C" */
T8_EXTERN_C_BEGIN ();

void
t8_forest_adapt (t8_forest_t forest)
{
  T8_ASSERT (forest != NULL);
  T8_ASSERT (forest->set_from != NULL);
  T8_ASSERT (forest->set_adapt_recursive != -1);
  const t8_forest_t forest_from = forest->set_from;

  if (t8_scheme_has_subelement_scheme (t8_forest_get_scheme (forest_from))) {
    SC_CHECK_ABORT (!forest->set_adapt_recursive || !t8_forest_has_subelements (forest_from),
                    "Recursive adaptation is currently not implemented for subelement schemes.");
    SC_CHECK_ABORT (!forest_from->incomplete_trees,
                    "Adapting forests with holes is currently not implemented for subelement schemes.");
    t8_forest_adapt_with<t8_subelement_adaptor> (forest);
  }
  else if (forest_from->incomplete_trees) {
    t8_forest_adapt_with<t8_incomplete_adaptor> (forest);
  }
  else {
    t8_forest_adapt_with<t8_standard_adaptor> (forest);
  }
}

T8_EXTERN_C_END ();
