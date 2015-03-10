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
#include <t8_element_hex.h>
#include <t8_element_tet.h>

struct t8_element
{
  int                 t8_element_dummy;
};

/* *INDENT-OFF* */
const int t8_type_to_dimension[T8_TYPE_LAST] =
  { 0, 1, 2, 2, 3, 3, 3, 3 };

static const int t8_type_boundary_count[T8_TYPE_LAST][T8_TYPE_LAST] =
  {{ 0,  0, 0, 0, 0, 0, 0, 0 },
   { 2,  0, 0, 0, 0, 0, 0, 0 },
   { 4,  4, 0, 0, 0, 0, 0, 0 },
   { 3,  3, 0, 0, 0, 0, 0, 0 },
   { 8, 12, 6, 0, 0, 0, 0, 0 },
   { 4,  6, 0, 4, 0, 0, 0, 0 },
   { 6,  9, 3, 2, 0, 0, 0, 0 },
   { 5,  8, 1, 4, 0, 0, 0, 0 }};
/* *INDENT-ON* */

t8_scheme_t        *
t8_scheme_new_default (void)
{
  t8_scheme_t        *s;

  s = T8_ALLOC_ZERO (t8_scheme_t, 1);
  s->type_schemes[T8_TYPE_QUAD] = t8_type_scheme_new_quad ();
  s->type_schemes[T8_TYPE_HEX] = t8_type_scheme_new_hex ();
  s->type_schemes[T8_TYPE_TET] = t8_type_scheme_new_tet ();

  return s;
}

void
t8_scheme_destroy (t8_scheme_t * s)
{
  int                 t;

  T8_ASSERT (s != NULL);

  for (t = 0; t < T8_TYPE_LAST; ++t) {
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

int
t8_type_count_boundary (t8_type_t thetype, int min_dim, int *per_type)
{
  int                 t;
  int                 sum;

  sum = 0;
  for (t = 0; t < T8_TYPE_LAST; ++t) {
    if (t8_type_to_dimension[t] >= min_dim) {
      sum += (per_type[t] = t8_type_boundary_count[thetype][t]);
    }
    else {
      per_type[t] = 0;
    }
  }

  return sum;
}

void
t8_type_boundary_alloc (t8_scheme_t * scheme, t8_type_t thetype,
                        int min_dim, int length, t8_element_t ** boundary)
{
  int                 t, offset, per;
#ifdef T8_ENABLE_DEBUG
  int                 per_type[T8_TYPE_LAST];
#endif

  T8_ASSERT (length == t8_type_count_boundary (thetype, min_dim, per_type));

  for (offset = t = 0; t < T8_TYPE_LAST; ++t) {
    if (t8_type_to_dimension[t] >= min_dim) {
      per = t8_type_boundary_count[thetype][t];
      if (per > 0) {
        t8_element_new (scheme->type_schemes[t], per, boundary + offset);
        offset += per;
      }
    }
  }
  T8_ASSERT (offset == length);
}

void
t8_type_boundary_destroy (t8_scheme_t * scheme, t8_type_t thetype,
                          int min_dim, int length, t8_element_t ** boundary)
{
  int                 t, offset, per;
#ifdef T8_ENABLE_DEBUG
  int                 per_type[T8_TYPE_LAST];
#endif

  T8_ASSERT (length == t8_type_count_boundary (thetype, min_dim, per_type));

  for (offset = t = 0; t < T8_TYPE_LAST; ++t) {
    if (t8_type_to_dimension[t] >= min_dim) {
      per = t8_type_boundary_count[thetype][t];
      if (per > 0) {
        t8_element_destroy (scheme->type_schemes[t], per, boundary + offset);
        offset += per;
      }
    }
  }
  T8_ASSERT (offset == length);
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
t8_element_child (t8_type_scheme_t * ts, const t8_element_t * elem,
                  int childid, t8_element_t * child)
{
  T8_ASSERT (ts != NULL && ts->elem_child != NULL);
  ts->elem_child (elem, childid, child);
}

void
t8_element_nca (t8_type_scheme_t * ts, const t8_element_t * elem1,
                const t8_element_t * elem2, t8_element_t * nca)
{
  T8_ASSERT (ts != NULL && ts->elem_nca != NULL);
  ts->elem_nca (elem1, elem2, nca);
}

void
t8_element_boundary (t8_type_scheme_t * ts,
                     const t8_element_t * elem,
                     int min_dim, int length, t8_element_t ** boundary)
{
  T8_ASSERT (ts != NULL && ts->elem_boundary != NULL);
  ts->elem_boundary (elem, min_dim, length, boundary);
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
