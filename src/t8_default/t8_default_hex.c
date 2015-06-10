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

#include <p8est_bits.h>
#include "t8_default_common.h"
#include "t8_default_hex.h"

static              size_t
t8_default_hex_size (void)
{
  return sizeof (t8_phex_t);
}

static int
t8_default_hex_maxlevel (void)
{
  return P8EST_QMAXLEVEL;
}

static              t8_eclass_t
t8_default_hex_child_eclass (int childid)
{
  T8_ASSERT (0 <= childid && childid < P8EST_CHILDREN);

  return T8_ECLASS_HEX;
}

static int
t8_default_hex_level (const t8_element_t * elem)
{
  return (int) ((const p8est_quadrant_t *) elem)->level;
}

static void
t8_default_hex_copy (const t8_element_t * source, t8_element_t * dest)
{
  *(p8est_quadrant_t *) dest = *(const p8est_quadrant_t *) source;
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

  T8_ASSERT (p8est_quadrant_is_extended (q));
  T8_ASSERT (q->level < P8EST_QMAXLEVEL);
  T8_ASSERT (0 <= childid && childid < P8EST_CHILDREN);

  r->x = childid & 0x01 ? (q->x | shift) : q->x;
  r->y = childid & 0x02 ? (q->y | shift) : q->y;
  r->z = childid & 0x04 ? (q->z | shift) : q->z;
  r->level = q->level + 1;
  T8_ASSERT (p8est_quadrant_is_parent (q, r));
}

static void
t8_default_hex_children (const t8_element_t * elem,
                         int length, t8_element_t * c[])
{
  T8_ASSERT (length == P8EST_CHILDREN);

  p8est_quadrant_childrenpv ((const p8est_quadrant_t *) elem,
                             (p8est_quadrant_t **) c);
}

static int
t8_default_hex_child_id (const t8_element_t * elem)
{
  return p8est_quadrant_child_id ((p8est_quadrant_t *) elem);
}

static int
t8_default_hex_is_family (const t8_element_t ** fam)
{
  return p8est_quadrant_is_familypv ((p8est_quadrant_t **) fam);
}

static void
t8_default_hex_set_linear_id (t8_element_t * elem, int level, uint64_t id)
{
  T8_ASSERT (0 <= level && level <= P8EST_QMAXLEVEL);
  T8_ASSERT (0 <= id && id < (uint64_t) 1 << P8EST_DIM * level);

  p8est_quadrant_set_morton ((p8est_quadrant_t *) elem, level, id);
}

static void
t8_default_hex_successor (const t8_element_t * elem1,
                          t8_element_t * elem2, int level)
{
  uint64_t            id;
  T8_ASSERT (0 <= level && level <= P8EST_QMAXLEVEL);

  id = p8est_quadrant_linear_id ((const p8est_quadrant_t *) elem1, level);
  T8_ASSERT (id + 1 < (1 << P8EST_DIM * level));
  p8est_quadrant_set_morton ((p8est_quadrant_t *) elem2, level, id + 1);
}

t8_eclass_scheme_t *
t8_default_scheme_new_hex (void)
{
  t8_eclass_scheme_t *ts;

  ts = T8_ALLOC (t8_eclass_scheme_t, 1);
  ts->eclass = T8_ECLASS_HEX;

  ts->elem_size = t8_default_hex_size;
  ts->elem_maxlevel = t8_default_hex_maxlevel;
  ts->elem_child_eclass = t8_default_hex_child_eclass;

  ts->elem_level = t8_default_hex_level;
  ts->elem_copy = t8_default_hex_copy;
  ts->elem_parent = (t8_element_parent_t) p8est_quadrant_parent;
  ts->elem_sibling = t8_default_hex_sibling;
  ts->elem_child = t8_default_hex_child;
  ts->elem_children = t8_default_hex_children;
  ts->elem_child_id = t8_default_hex_child_id;
  ts->elem_is_family = t8_default_hex_is_family;
  ts->elem_nca = (t8_element_nca_t) p8est_nearest_common_ancestor;
  ts->elem_boundary = NULL;
  ts->elem_set_linear_id = t8_default_hex_set_linear_id;
  ts->elem_successor = t8_default_hex_successor;

  ts->elem_new = t8_default_mempool_alloc;
  ts->elem_destroy = t8_default_mempool_free;

  ts->ts_destroy = t8_default_scheme_mempool_destroy;
  ts->ts_context = sc_mempool_new (sizeof (t8_phex_t));

  return ts;
}
