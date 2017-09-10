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
t8_element_copy (t8_eclass_scheme_t * ts, const t8_element_t * source,
                 t8_element_t * dest)
{
  T8_ASSERT (ts != NULL && ts->elem_copy != NULL);
  ts->elem_copy (source, dest);
}

int
t8_element_compare (t8_eclass_scheme_t * ts, const t8_element_t * elem1,
                    const t8_element_t * elem2)
{
  T8_ASSERT (ts != NULL && ts->elem_compare != NULL);

  return ts->elem_compare (elem1, elem2);
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

int
t8_element_child_id (t8_eclass_scheme_t * ts, const t8_element_t * elem)
{
  T8_ASSERT (ts != NULL && ts->elem_child_id != NULL);
  return ts->elem_child_id (elem);
}

int
t8_element_is_family (t8_eclass_scheme_t * ts, t8_element_t ** fam)
{
  T8_ASSERT (ts != NULL && ts->elem_is_family != NULL);

  return ts->elem_is_family (fam);
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
                          t8_element_t * elem, int level, t8_linearidx_t id)
{
  T8_ASSERT (ts != NULL && ts->elem_set_linear_id != NULL);

  ts->elem_set_linear_id (elem, level, id);
}

t8_linearidx_t
t8_element_get_linear_id (t8_eclass_scheme_t * ts,
                          const t8_element_t * elem, int level)
{
  T8_ASSERT (ts != NULL && ts->elem_get_linear_id != NULL);

  return ts->elem_get_linear_id (elem, level);
}

void
t8_element_first_descendant (t8_eclass_scheme_t * ts,
                             const t8_element_t * elem, t8_element_t * desc)
{
  T8_ASSERT (ts != NULL && ts->elem_first_desc != NULL);

  ts->elem_first_desc (elem, desc);
}

void
t8_element_last_descendant (t8_eclass_scheme_t * ts,
                            const t8_element_t * elem, t8_element_t * desc)
{
  T8_ASSERT (ts != NULL && ts->elem_last_desc != NULL);

  ts->elem_last_desc (elem, desc);
}

void
t8_element_successor (t8_eclass_scheme_t * ts, const t8_element_t * elem1,
                      t8_element_t * elem2, int level)
{
  T8_ASSERT (ts != NULL && ts->elem_successor != NULL);

  ts->elem_successor (elem1, elem2, level);
}

void
t8_element_anchor (t8_eclass_scheme_t * ts, const t8_element_t * elem,
                   int anchor[3])
{
  T8_ASSERT (ts != NULL && ts->elem_anchor != NULL);

  ts->elem_anchor (elem, anchor);
}

int
t8_element_root_len (t8_eclass_scheme_t * ts, const t8_element_t * elem)
{
  T8_ASSERT (ts != NULL && ts->elem_root_len != NULL);
  return ts->elem_root_len (elem);
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

t8_element_t       *
t8_element_array_index (t8_eclass_scheme_t * ts, sc_array_t * array,
                        size_t it)
{
  T8_ASSERT (ts != NULL && ts->elem_size != NULL);
  T8_ASSERT (array->elem_size == ts->elem_size ());
  T8_ASSERT (it < array->elem_count);

  return (t8_element_t *) (array->array + array->elem_size * it);
}
#endif /* if 0 */
