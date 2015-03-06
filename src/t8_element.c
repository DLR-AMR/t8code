/*
  This file is part of t8code.
  t8code is a C library to manage a collection (a forest) of multiple
  connected adaptive space-trees of general element types in parallel.

  Copyright (C) 2010 The University of Texas System
  Written by Carsten Burstedde, Lucas C. Wilcox, and Tobin Isaac

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

#include <t8_element_quad.h>

struct t8_element
{
  int                 t8_element_dummy;
};

t8_scheme_t        *
t8_type_scheme_new_default (void)
{
  t8_scheme_t        *s;

  s = T8_ALLOC_ZERO (t8_scheme_t, 1);
  s->type_schemes[T8_TYPE_QUAD] = t8_type_scheme_new_quad ();

  return s;
}

void
t8_escheme_destroy_default (t8_scheme_t * s)
{
  int                 t;

  T8_ASSERT (s != NULL);

  for (t = T8_TYPE_FIRST; t < T8_TYPE_LAST; ++t) {
    if (s->type_schemes[t] != NULL) {
      t8_type_scheme_destroy (s->type_schemes[t]);
    }
  }
  T8_FREE (s);
}

void
t8_type_scheme_destroy (t8_type_scheme_t * ts)
{
  T8_ASSERT (ts != NULL);

  if (ts->ts_destroy != NULL) {
    ts->ts_destroy (ts);
  }
  T8_FREE (ts);
}

void
t8_element_parent (t8_type_scheme_t * ts,
                   const t8_element_t * elem, t8_element_t * parent)
{
  T8_ASSERT (ts != NULL && ts->elem_parent != NULL);
  ts->elem_parent (elem, parent);
}

void
t8_element_sibling (t8_type_scheme_t * ts,
                    const t8_element_t * elem, int sibid,
                    t8_element_t * sibling)
{
  T8_ASSERT (ts != NULL && ts->elem_sibling != NULL);
  ts->elem_sibling (elem, sibid, sibling);
}

void
t8_type_scheme_mempool_destroy (t8_type_scheme_t * ts)
{
  T8_ASSERT (ts->ts_context != NULL);
  sc_mempool_destroy ((sc_mempool_t *) ts->ts_context);
}

void
t8_element_mempool_new (void *ts_context, int length, t8_element_t ** elem)
{
  int                 i;

  T8_ASSERT (ts_context != NULL);
  T8_ASSERT (0 <= length);
  T8_ASSERT (elem != NULL);

  for (i = 0; i < length; ++i) {
    elem[i] = (t8_element_t *) sc_mempool_alloc ((sc_mempool_t *) ts_context);
  }
}

void
t8_element_mempool_destroy (void *ts_context, int length,
                            t8_element_t ** elem)
{
  int                 i;

  T8_ASSERT (ts_context != NULL);
  T8_ASSERT (0 <= length);
  T8_ASSERT (elem != NULL);

  for (i = 0; i < length; ++i) {
    sc_mempool_free ((sc_mempool_t *) ts_context, elem[i]);
  }
}

void
t8_element_new (t8_type_scheme_t * ts, int length, t8_element_t ** elems)
{
  T8_ASSERT (ts != NULL && ts->elem_new != NULL);
  ts->elem_new (ts->ts_context, length, elems);
}

void
t8_element_destroy (t8_type_scheme_t * ts, int length, t8_element_t ** elems)
{
  T8_ASSERT (ts != NULL && ts->elem_destroy != NULL);
  ts->elem_destroy (ts->ts_context, length, elems);
}
