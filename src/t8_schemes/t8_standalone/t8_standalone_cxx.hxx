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

#ifndef T8_STANDALONE_CXX_HXX
#define T8_STANDALONE_CXX_HXX

#include <new>
#include <t8_refcount.h>
#include <t8_eclass.h>
#include <t8_schemes/t8_standalone/t8_standalone_element_cxx.hxx>

/** Return the standalone element implementation of t8code. */
t8_scheme_cxx_t    *
t8_scheme_new_standalone_cxx (void)
{
  t8_scheme_cxx_t    *s;

  s = T8_ALLOC_ZERO (t8_scheme_cxx_t, 1);
  t8_refcount_init (&s->rc);

  s->eclass_schemes[T8_ECLASS_VERTEX]   = new t8_standalone_scheme_c <T8_ECLASS_VERTEX> ();
  s->eclass_schemes[T8_ECLASS_LINE]     = new t8_standalone_scheme_c <T8_ECLASS_LINE> ();
  s->eclass_schemes[T8_ECLASS_QUAD]     = new t8_standalone_scheme_c <T8_ECLASS_QUAD> ();
  s->eclass_schemes[T8_ECLASS_HEX]      = new t8_standalone_scheme_c <T8_ECLASS_HEX> ();
  s->eclass_schemes[T8_ECLASS_TRIANGLE] = new t8_standalone_scheme_c <T8_ECLASS_TRIANGLE> ();
  s->eclass_schemes[T8_ECLASS_TET]      = new t8_standalone_scheme_c <T8_ECLASS_TET> ();
  s->eclass_schemes[T8_ECLASS_PRISM]    = new t8_standalone_scheme_c <T8_ECLASS_PRISM> ();
  s->eclass_schemes[T8_ECLASS_PYRAMID]  = new t8_standalone_scheme_c <T8_ECLASS_PYRAMID> ();

  return s;
}

#endif /* !T8_STANDALONE_CXX_HXX */
