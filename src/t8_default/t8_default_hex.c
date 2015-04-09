/*
  This file is part of t8code.
  t8code is a C library to manage a collection (a forest) of multiple
  connected adaptive space-trees of general element classes in parallel.

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

#include <p8est_bits.h>
#include "t8_default_common.h"
#include "t8_default_hex.h"

static              size_t
t8_default_hex_size (void)
{
  return sizeof (t8_phex_t);
}

static void
t8_default_hex_sibling (const t8_element_t * elem,
                        int sibid, t8_element_t * sibling)
{
  p8est_quadrant_sibling ((const p8est_quadrant_t *) elem,
                          (p8est_quadrant_t *) sibling, sibid);
}

static void
t8_default_hex_child (const t8_element_t * elem,
                      int childid, t8_element_t * child)
{
  const p8est_quadrant_t *q = (const p8est_quadrant_t *) elem;
  const p4est_qcoord_t shift = P8EST_QUADRANT_LEN (q->level + 1);
  p8est_quadrant_t   *r = (p8est_quadrant_t *) child;

  P4EST_ASSERT (p8est_quadrant_is_extended (q));
  P4EST_ASSERT (q->level < P8EST_QMAXLEVEL);
  P4EST_ASSERT (childid >= 0 && childid < P8EST_CHILDREN);

  r->x = childid & 0x01 ? (q->x | shift) : q->x;
  r->y = childid & 0x02 ? (q->y | shift) : q->y;
  r->z = childid & 0x04 ? (q->y | shift) : q->z;
  r->level = q->level + 1;
  P4EST_ASSERT (p8est_quadrant_is_parent (q, r));
}

t8_eclass_scheme_t *
t8_default_scheme_new_hex (void)
{
  t8_eclass_scheme_t *ts;

  ts = T8_ALLOC (t8_eclass_scheme_t, 1);

  ts->elem_size = t8_default_hex_size;

  ts->elem_parent = (t8_element_parent_t) p8est_quadrant_parent;
  ts->elem_sibling = t8_default_hex_sibling;
  ts->elem_child = t8_default_hex_child;
  ts->elem_nca = (t8_element_nca_t) p8est_nearest_common_ancestor;

  ts->elem_new = t8_default_mempool_alloc;
  ts->elem_destroy = t8_default_mempool_free;
  ts->ts_destroy = t8_default_scheme_mempool_destroy;
  ts->ts_context = sc_mempool_new (sizeof (t8_phex_t));

  return ts;
}
