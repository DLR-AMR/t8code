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

#include <new>
#include <t8_schemes/t8_default/t8_default.hxx>
#include <t8_refcount.h>

/* We want to export the whole implementation to be callable from "C" */
T8_EXTERN_C_BEGIN ();

t8_scheme *
t8_scheme_new_default_cxx (void)
{
  t8_scheme *s;

  s = T8_ALLOC_ZERO (t8_scheme, 1);
  t8_refcount_init (&s->rc);

  s->eclass_schemes[T8_ECLASS_VERTEX] = new t8_default_scheme_vertex ();
  s->eclass_schemes[T8_ECLASS_LINE] = new t8_default_scheme_line ();
  s->eclass_schemes[T8_ECLASS_QUAD] = new t8_default_scheme_quad ();
  s->eclass_schemes[T8_ECLASS_HEX] = new t8_default_scheme_hex ();
  s->eclass_schemes[T8_ECLASS_TRIANGLE] = new t8_default_scheme_tri ();
  s->eclass_schemes[T8_ECLASS_TET] = new t8_default_scheme_tet ();
  s->eclass_schemes[T8_ECLASS_PRISM] = new t8_default_scheme_prism ();
  s->eclass_schemes[T8_ECLASS_PYRAMID] = new t8_default_scheme_pyramid ();

  T8_ASSERT (s->eclass_schemes[T8_ECLASS_VERTEX]->t8_element_maxlevel ()
             >= s->eclass_schemes[T8_ECLASS_LINE]->t8_element_maxlevel ());
  T8_ASSERT (s->eclass_schemes[T8_ECLASS_LINE]->t8_element_maxlevel ()
             >= s->eclass_schemes[T8_ECLASS_QUAD]->t8_element_maxlevel ());
  T8_ASSERT (s->eclass_schemes[T8_ECLASS_LINE]->t8_element_maxlevel ()
             >= s->eclass_schemes[T8_ECLASS_TRIANGLE]->t8_element_maxlevel ());
  T8_ASSERT (s->eclass_schemes[T8_ECLASS_TRIANGLE]->t8_element_maxlevel ()
             >= s->eclass_schemes[T8_ECLASS_TET]->t8_element_maxlevel ());
  T8_ASSERT (s->eclass_schemes[T8_ECLASS_TRIANGLE]->t8_element_maxlevel ()
             >= s->eclass_schemes[T8_ECLASS_PRISM]->t8_element_maxlevel ());
  T8_ASSERT (s->eclass_schemes[T8_ECLASS_TRIANGLE]->t8_element_maxlevel ()
             >= s->eclass_schemes[T8_ECLASS_PYRAMID]->t8_element_maxlevel ());
  T8_ASSERT (s->eclass_schemes[T8_ECLASS_QUAD]->t8_element_maxlevel ()
             >= s->eclass_schemes[T8_ECLASS_HEX]->t8_element_maxlevel ());
  T8_ASSERT (s->eclass_schemes[T8_ECLASS_QUAD]->t8_element_maxlevel ()
             >= s->eclass_schemes[T8_ECLASS_PRISM]->t8_element_maxlevel ());
  T8_ASSERT (s->eclass_schemes[T8_ECLASS_QUAD]->t8_element_maxlevel ()
             >= s->eclass_schemes[T8_ECLASS_PYRAMID]->t8_element_maxlevel ());

  return s;
}

int
t8_eclass_scheme_is_default (t8_scheme *ts)
{
  switch (ts->eclass) {
  case T8_ECLASS_VERTEX:
    return T8_COMMON_IS_TYPE (ts, t8_default_scheme_vertex *);
  case T8_ECLASS_LINE:
    return T8_COMMON_IS_TYPE (ts, t8_default_scheme_line *);
  case T8_ECLASS_QUAD:
    return T8_COMMON_IS_TYPE (ts, t8_default_scheme_quad *);
  case T8_ECLASS_TRIANGLE:
    return T8_COMMON_IS_TYPE (ts, t8_default_scheme_tri *);
  case T8_ECLASS_HEX:
    return T8_COMMON_IS_TYPE (ts, t8_default_scheme_hex *);
  case T8_ECLASS_TET:
    return T8_COMMON_IS_TYPE (ts, t8_default_scheme_tet *);
  case T8_ECLASS_PRISM:
    return T8_COMMON_IS_TYPE (ts, t8_default_scheme_prism *);
  case T8_ECLASS_PYRAMID:
    return T8_COMMON_IS_TYPE (ts, t8_default_scheme_pyramid *);
  default:
    SC_ABORT_NOT_REACHED ();
  }
  return 0; /* Default return value false */
}

T8_EXTERN_C_END ();
