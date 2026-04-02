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

/** \file t8_standalone.cxx
 *  Implements functions declared in \ref t8_standalone.hxx.
 */

#include <t8_schemes/t8_standalone/t8_standalone.hxx>
#include <t8_schemes/t8_scheme_builder.hxx>

/* We want to export the whole implementation to be callable from "C" */
T8_EXTERN_C_BEGIN ();

const t8_scheme *
t8_scheme_new_standalone (void)
{
  t8_scheme_builder builder;

  builder.add_eclass_scheme<t8_standalone_scheme<T8_ECLASS_VERTEX>> ();
  builder.add_eclass_scheme<t8_standalone_scheme<T8_ECLASS_LINE>> ();
  builder.add_eclass_scheme<t8_standalone_scheme<T8_ECLASS_QUAD>> ();
  builder.add_eclass_scheme<t8_standalone_scheme<T8_ECLASS_TRIANGLE>> ();
  builder.add_eclass_scheme<t8_standalone_scheme<T8_ECLASS_HEX>> ();
  builder.add_eclass_scheme<t8_standalone_scheme<T8_ECLASS_TET>> ();
  builder.add_eclass_scheme<t8_standalone_scheme<T8_ECLASS_PRISM>> ();
  builder.add_eclass_scheme<t8_standalone_scheme<T8_ECLASS_PYRAMID>> ();
  // builder.add_eclass_scheme<t8_standalone_scheme<T8_ECLASS_PYRAMID>> (); //ADD if PYRAMID in "element_get_shape"
  return builder.build_scheme ();
}

int
t8_eclass_scheme_is_standalone (const t8_scheme *scheme, const t8_eclass_t eclass)
{
  switch (eclass) {
  case T8_ECLASS_VERTEX:
    return scheme->check_eclass_scheme_type<t8_standalone_scheme<T8_ECLASS_VERTEX>> (T8_ECLASS_VERTEX);
  case T8_ECLASS_LINE:
    return scheme->check_eclass_scheme_type<t8_standalone_scheme<T8_ECLASS_LINE>> (T8_ECLASS_LINE);
  case T8_ECLASS_QUAD:
    return scheme->check_eclass_scheme_type<t8_standalone_scheme<T8_ECLASS_QUAD>> (T8_ECLASS_QUAD);
  // case T8_ECLASS_TRIANGLE:
  //   return scheme->check_eclass_scheme_type<t8_standalone_scheme<T8_ECLASS_TRIANGLE>> (T8_ECLASS_TRIANGLE);
  case T8_ECLASS_HEX:
    return scheme->check_eclass_scheme_type<t8_standalone_scheme<T8_ECLASS_HEX>> (T8_ECLASS_HEX);
  // case T8_ECLASS_TET:
  //   return scheme->check_eclass_scheme_type<t8_standalone_scheme<T8_ECLASS_PRISM>> (T8_ECLASS_TET);
  // case T8_ECLASS_PRISM:
  //   return scheme->check_eclass_scheme_type<t8_standalone_scheme<T8_ECLASS_PRISM>> (T8_ECLASS_PRISM);
  // case T8_ECLASS_PYRAMID:
  //   return scheme->check_eclass_scheme_type<t8_standalone_scheme<T8_ECLASS_PYRAMID>> (T8_ECLASS_PYRAMID);
  default:
    SC_ABORT_NOT_REACHED ();
  }
  return 0; /* Default return value false */
}

T8_EXTERN_C_END ();
