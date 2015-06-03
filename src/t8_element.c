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

#include <t8_element.h>

static void
t8_scheme_destroy (t8_scheme_t * s)
{
  int                 t;

  T8_ASSERT (s != NULL);
  T8_ASSERT (s->rc.refcount == 0);

  for (t = 0; t < T8_ECLASS_LAST; ++t) {
    if (s->eclass_schemes[t] != NULL) {
      t8_eclass_scheme_destroy (s->eclass_schemes[t]);
    }
  }
  T8_FREE (s);
}

void
t8_scheme_ref (t8_scheme_t * scheme)
{
  T8_ASSERT (scheme != NULL);

  sc_refcount_ref (&scheme->rc);
}

void
t8_scheme_unref (t8_scheme_t ** pscheme)
{
  t8_scheme_t        *scheme;

  T8_ASSERT (pscheme != NULL);
  scheme = *pscheme;
  T8_ASSERT (scheme != NULL);

  if (sc_refcount_unref (&scheme->rc)) {
    t8_scheme_destroy (scheme);
    *pscheme = NULL;
  }
}

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
  int                 per_eclass[T8_ECLASS_LAST];
#endif

  T8_ASSERT (length ==
             t8_eclass_count_boundary (theclass, min_dim, per_eclass));

  for (offset = t = 0; t < T8_ECLASS_LAST; ++t) {
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
  int                 per_eclass[T8_ECLASS_LAST];
#endif

  T8_ASSERT (length ==
             t8_eclass_count_boundary (theclass, min_dim, per_eclass));

  for (offset = t = 0; t < T8_ECLASS_LAST; ++t) {
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

int
t8_element_maxlevel (t8_eclass_scheme_t * ts)
{
  T8_ASSERT (ts != NULL && ts->elem_maxlevel != NULL);
  return ts->elem_maxlevel ();
}

t8_eclass_t
t8_element_child_eclass (t8_eclass_scheme_t * ts, int childid)
{
  T8_ASSERT (ts != NULL && ts->elem_child_eclass != NULL);
  T8_ASSERT (0 <= childid && childid < t8_eclass_num_children[ts->eclass]);
  return ts->elem_child_eclass (childid);
}

int
t8_element_level (t8_eclass_scheme_t * ts, const t8_element_t * elem)
{
  T8_ASSERT (ts != NULL && ts->elem_level != NULL);
  return ts->elem_level (elem);
}

void
t8_element_parent (t8_eclass_scheme_t * ts,
                   const t8_element_t * elem, t8_element_t * parent)
{
  T8_ASSERT (ts != NULL && ts->elem_parent != NULL);
  T8_ASSERT (t8_element_level (ts, elem) > 0);
  ts->elem_parent (elem, parent);
}

void
t8_element_sibling (t8_eclass_scheme_t * ts,
                    const t8_element_t * elem, int sibid,
                    t8_element_t * sibling)
{
  T8_ASSERT (ts != NULL && ts->elem_sibling != NULL);
  T8_ASSERT (0 <= sibid && sibid < t8_eclass_num_children[ts->eclass]);
  ts->elem_sibling (elem, sibid, sibling);
}

void
t8_element_child (t8_eclass_scheme_t * ts, const t8_element_t * elem,
                  int childid, t8_element_t * child)
{
  T8_ASSERT (ts != NULL && ts->elem_child != NULL);
  T8_ASSERT (t8_element_level (ts, elem) < t8_element_maxlevel (ts));
  T8_ASSERT (0 <= childid && childid < t8_eclass_num_children[ts->eclass]);
  ts->elem_child (elem, childid, child);
}

void
t8_element_children (t8_eclass_scheme_t * ts, const t8_element_t * elem,
                     int length, t8_element_t * c[])
{
  T8_ASSERT (ts != NULL && ts->elem_children != NULL);
  T8_ASSERT (t8_element_level (ts, elem) < t8_element_maxlevel (ts));
  T8_ASSERT (length == t8_eclass_num_children[ts->eclass]);
  ts->elem_children (elem, length, c);
}

void
t8_element_nca (t8_eclass_scheme_t * ts, const t8_element_t * elem1,
                const t8_element_t * elem2, t8_element_t * nca)
{
  T8_ASSERT (ts != NULL && ts->elem_nca != NULL);
  ts->elem_nca (elem1, elem2, nca);
}

void
t8_element_boundary (t8_eclass_scheme_t * ts,
                     const t8_element_t * elem,
                     int min_dim, int length, t8_element_t ** boundary)
{
  T8_ASSERT (ts != NULL && ts->elem_boundary != NULL);
  ts->elem_boundary (elem, min_dim, length, boundary);
}

void
t8_element_set_linear_id (t8_eclass_scheme_t * ts,
                          t8_element_t * elem, int level, uint64_t id)
{
  T8_ASSERT (ts != NULL && ts->elem_set_linear_id != NULL);

  ts->elem_set_linear_id (elem, level, id);
}

void
t8_element_new (t8_eclass_scheme_t * ts, int length, t8_element_t ** elems)
{
  T8_ASSERT (ts != NULL && ts->elem_new != NULL);
  ts->elem_new (ts->ts_context, length, elems);
}

void
t8_element_destroy (t8_eclass_scheme_t * ts, int length,
                    t8_element_t ** elems)
{
  T8_ASSERT (ts != NULL && ts->elem_destroy != NULL);
  ts->elem_destroy (ts->ts_context, length, elems);
}
