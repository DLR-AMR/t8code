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

#include <p4est_bits.h>
#include <t8_element_quad.h>

static void
t8_element_quad_sibling (const t8_element_t * elem,
                         int sibid, t8_element_t * sibling)
{
  p4est_quadrant_sibling ((const p4est_quadrant_t *) elem,
                          (p4est_quadrant_t *) sibling, sibid);
}

static void
t8_element_quad_child (const t8_element_t * elem,
                       int childid, t8_element_t * child)
{
  const p4est_quadrant_t *q = (const p4est_quadrant_t *) elem;
  const p4est_qcoord_t shift = P4EST_QUADRANT_LEN (q->level + 1);
  p4est_quadrant_t   *r = (p4est_quadrant_t *) child;

  P4EST_ASSERT (p4est_quadrant_is_extended (q));
  P4EST_ASSERT (q->level < P4EST_QMAXLEVEL);
  P4EST_ASSERT (childid >= 0 && childid < P4EST_CHILDREN);

  r->x = childid & 0x01 ? (q->x | shift) : q->x;
  r->y = childid & 0x02 ? (q->y | shift) : q->y;
  r->level = q->level + 1;
  P4EST_ASSERT (p4est_quadrant_is_parent (q, r));
}

static void
t8_element_quad_new (void *ts_context, int length, t8_element_t ** elem)
{
  int                 i;

  T8_ASSERT (ts_context != NULL);
  T8_ASSERT (0 <= length);
  T8_ASSERT (elem != NULL);

  for (i = 0; i < length; ++i) {
    elem[i] = (t8_element_t *) sc_mempool_alloc ((sc_mempool_t *) ts_context);
  }
}

static void
t8_element_quad_destroy (void *ts_context, int length, t8_element_t ** elem)
{
  int                 i;

  T8_ASSERT (ts_context != NULL);
  T8_ASSERT (0 <= length);
  T8_ASSERT (elem != NULL);

  for (i = 0; i < length; ++i) {
    sc_mempool_free ((sc_mempool_t *) ts_context, elem[i]);
  }
}

static void
t8_type_scheme_quad_destroy (t8_type_scheme_t * ts)
{
  T8_ASSERT (ts->ts_context != NULL);

  sc_mempool_destroy ((sc_mempool_t *) ts->ts_context);
}

t8_type_scheme_t   *
t8_type_scheme_new_quad (void)
{
  t8_type_scheme_t   *ts;

  ts = T8_ALLOC (t8_type_scheme_t, 1);

  ts->elem_parent = (t8_element_parent_t) p4est_quadrant_parent;
  ts->elem_sibling = t8_element_quad_sibling;
  ts->elem_child = t8_element_quad_child;
  ts->elem_nca = (t8_element_nca_t) p4est_nearest_common_ancestor;

  ts->elem_new = t8_element_quad_new;
  ts->elem_destroy = t8_element_quad_destroy;

  ts->ts_destroy = t8_type_scheme_quad_destroy;
  ts->ts_context = sc_mempool_new (sizeof (p4est_quadrant_t));

  return ts;
}
