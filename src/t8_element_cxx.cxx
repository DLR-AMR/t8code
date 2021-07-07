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

#include <t8_element_cxx.hxx>

/* We want to export the whole implementation to be callable from "C" */
T8_EXTERN_C_BEGIN ();

#if 0
t8_eclass_scheme_c::~t8_eclass_scheme_c ()
{
}
#endif

/* This belongs here since it uses c++ function,
 * see t8_element.c/.h */
void
t8_scheme_cxx_destroy (t8_scheme_cxx_t * s)
{
  int                 t;

  T8_ASSERT (s != NULL);
  T8_ASSERT (s->rc.refcount == 0);

  for (t = 0; t < T8_ECLASS_COUNT; ++t) {
    if (s->eclass_schemes[t] != NULL) {
      delete              s->eclass_schemes[t];
    }
  }
  T8_FREE (s);
}

/* *INDENT-OFF* */
/* Default implementation for the element size */
size_t
t8_eclass_scheme::t8_element_size ()
{
  return element_size;
}
/* *INDENT-ON* */

/* Default implementation for array_index */
t8_element_t       *
t8_eclass_scheme::t8_element_array_index (sc_array_t * array, size_t it)
{
  T8_ASSERT (it < array->elem_count);
  T8_ASSERT (element_size == array->elem_size);
  return (t8_element_t *) sc_array_index (array, it);
}

T8_EXTERN_C_END ();

#if 0
void
t8_eclass_scheme_destroy (t8_eclass_scheme_t * ts)
{
  T8_ASSERT (ts != NULL);

  if (ts->ts_destroy != NULL) {
    ts->ts_destroy (ts);
  }
  T8_FREE (ts);
}

void
t8_eclass_boundary_alloc (t8_scheme_t * scheme, t8_eclass_t theclass,
                          int min_dim, int length, t8_element_t ** boundary)
{
  int                 t, offset, per;
#ifdef T8_ENABLE_DEBUG
  int                 per_eclass[T8_ECLASS_COUNT];
#endif

  T8_ASSERT (length ==
             t8_eclass_count_boundary (theclass, min_dim, per_eclass));

  for (offset = t = 0; t < T8_ECLASS_COUNT; ++t) {
    if (t8_eclass_to_dimension[t] >= min_dim) {
      per = t8_eclass_boundary_count[theclass][t];
      if (per > 0) {
        t8_element_new (scheme->eclass_schemes[t], per, boundary + offset);
        offset += per;
      }
    }
  }
  T8_ASSERT (offset == length);
}

void
t8_eclass_boundary_destroy (t8_scheme_t * scheme, t8_eclass_t theclass,
                            int min_dim, int length, t8_element_t ** boundary)
{
  int                 t, offset, per;
#ifdef T8_ENABLE_DEBUG
  int                 per_eclass[T8_ECLASS_COUNT];
#endif

  T8_ASSERT (length ==
             t8_eclass_count_boundary (theclass, min_dim, per_eclass));

  for (offset = t = 0; t < T8_ECLASS_COUNT; ++t) {
    if (t8_eclass_to_dimension[t] >= min_dim) {
      per = t8_eclass_boundary_count[theclass][t];
      if (per > 0) {
        t8_element_destroy (scheme->eclass_schemes[t], per,
                            boundary + offset);
        offset += per;
      }
    }
  }
  T8_ASSERT (offset == length);
}

size_t
t8_element_size (t8_eclass_scheme_t * ts)
{
  T8_ASSERT (ts != NULL && ts->elem_size != NULL);
  return ts->elem_size ();
}

#endif /* if 0 */
