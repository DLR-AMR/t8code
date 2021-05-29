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
#include <t8_schemes/t8_new_feature/t8_subelements_cxx.hxx>
#include <t8_refcount.h>

#include "t8_subelements_quad_cxx.hxx"

/* We want to export the whole implementation to be callable from "C" */
T8_EXTERN_C_BEGIN ();

t8_scheme_cxx_t    *
t8_scheme_new_subelement_cxx (void)
{
  t8_scheme_cxx_t    *s;

  s = T8_ALLOC_ZERO (t8_scheme_cxx_t, 1);
  t8_refcount_init (&s->rc);

  s->eclass_schemes[T8_ECLASS_VERTEX] = NULL;
  s->eclass_schemes[T8_ECLASS_LINE] = NULL;
  s->eclass_schemes[T8_ECLASS_QUAD] = new t8_default_scheme_sub_c ();
  s->eclass_schemes[T8_ECLASS_HEX] = NULL;
  s->eclass_schemes[T8_ECLASS_TRIANGLE] = NULL;
  s->eclass_schemes[T8_ECLASS_TET] = NULL;
  s->eclass_schemes[T8_ECLASS_PRISM] = NULL;
  /* s->eclass_schemes[T8_ECLASS_SUB] = new t8_default_scheme_sub_c (); */
  /* should be implemented in the future -> change t8_eclass.c ? */

  return s;
}

int
t8_eclass_scheme_is_default (t8_eclass_scheme_c * ts)
{
  switch (ts->eclass) {
  case T8_ECLASS_QUAD:
    return T8_COMMON_IS_TYPE (ts, t8_default_scheme_sub_c *);
  default:
    SC_ABORT_NOT_REACHED ();
    /* TODO: Add pyramid as soon as pyramid scheme is implemented */
    /* TODO: Add a test for this function */
  }
  return 0;                     /* Default return value false */
}

T8_EXTERN_C_END ();
