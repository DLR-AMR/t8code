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
#include "t8_default_common_cxx.hxx"
#include "t8_default_hex_cxx.hxx"


/* *INDENT-OFF* */
size_t
t8_default_scheme_hex_c::t8_element_size (void)
{
  return sizeof (t8_phex_t);
}
/* *INDENT-ON* */

int
t8_default_scheme_hex_c::t8_element_maxlevel (void)
{
  return P8EST_QMAXLEVEL;
}

/* *INDENT-OFF* */
t8_eclass_t
t8_default_scheme_hex_c::t8_element_child_eclass (int childid)
{
  T8_ASSERT (0 <= childid && childid < P8EST_CHILDREN);

  return T8_ECLASS_HEX;
}
/* *INDENT-ON* */

int
t8_default_scheme_hex_c::t8_element_level (const t8_element_t * elem)
{
  return (int) ((const p8est_quadrant_t *) elem)->level;
}

void
t8_default_scheme_hex_c::t8_element_copy (const t8_element_t * source,
                                          t8_element_t * dest)
{
  *(p8est_quadrant_t *) dest = *(const p8est_quadrant_t *) source;
}

int
t8_default_scheme_hex_c::t8_element_compare (const t8_element_t * elem1,
                                             const t8_element_t * elem2)
{
  int                 maxlvl;
  u_int64_t           id1, id2;

  /* Compute the bigger level of the two */
  maxlvl = SC_MAX (t8_default_scheme_hex_c::t8_element_level (elem1),
                   t8_default_scheme_hex_c::t8_element_level (elem2));
  /* Compute the linear ids of the elements */
  id1 = t8_default_scheme_hex_c::t8_element_get_linear_id (elem1, maxlvl);
  id2 = t8_default_scheme_hex_c::t8_element_get_linear_id (elem2, maxlvl);
  /* return negativ if id1 < id2, zero if id1 = id2, positive if id1 > id2 */
  return id1 < id2 ? -1 : id1 != id2;
}

void
t8_default_scheme_hex_c::t8_element_parent (const t8_element_t * elem,
                                            t8_element_t * parent)
{
  p8est_quadrant_parent ((const p8est_quadrant_t *) elem,
                         (p8est_quadrant_t *) parent);
}

void
t8_default_scheme_hex_c::t8_element_sibling (const t8_element_t * elem,
                                             int sibid,
                                             t8_element_t * sibling)
{
  p8est_quadrant_sibling ((const p8est_quadrant_t *) elem,
                          (p8est_quadrant_t *) sibling, sibid);
}

void
t8_default_scheme_hex_c::t8_element_child (const t8_element_t * elem,
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

void
t8_default_scheme_hex_c::t8_element_children (const t8_element_t * elem,
                                              int length, t8_element_t * c[])
{
  T8_ASSERT (length == P8EST_CHILDREN);

  p8est_quadrant_childrenpv ((const p8est_quadrant_t *) elem,
                             (p8est_quadrant_t **) c);
}

int
t8_default_scheme_hex_c::t8_element_child_id (const t8_element_t * elem)
{
  return p8est_quadrant_child_id ((p8est_quadrant_t *) elem);
}

int
t8_default_scheme_hex_c::t8_element_is_family (t8_element_t ** fam)
{
  return p8est_quadrant_is_familypv ((p8est_quadrant_t **) fam);
}

void
t8_default_scheme_hex_c::t8_element_nca (const t8_element_t * elem1,
                                         const t8_element_t * elem2,
                                         t8_element_t * nca)
{
  p8est_nearest_common_ancestor ((const p8est_quadrant_t *) elem1,
                                 (const p8est_quadrant_t *) elem2,
                                 (p8est_quadrant_t *) nca);
}

void
t8_default_scheme_hex_c::t8_element_set_linear_id (t8_element_t * elem,
                                                   int level, uint64_t id)
{
  T8_ASSERT (0 <= level && level <= P8EST_QMAXLEVEL);
  T8_ASSERT (0 <= id && id < ((uint64_t) 1) << P8EST_DIM * level);

  p8est_quadrant_set_morton ((p8est_quadrant_t *) elem, level, id);
}

uint64_t
  t8_default_scheme_hex_c::t8_element_get_linear_id (const t8_element_t *
                                                     elem, int level)
{
  T8_ASSERT (0 <= level && level <= P8EST_QMAXLEVEL);

  return p8est_quadrant_linear_id ((p8est_quadrant_t *) elem, level);
}

void
t8_default_scheme_hex_c::t8_element_first_descendant (const t8_element_t *
                                                      elem,
                                                      t8_element_t * desc)
{
  p8est_quadrant_first_descendant ((p8est_quadrant_t *) elem,
                                   (p8est_quadrant_t *) desc,
                                   P8EST_QMAXLEVEL);
}

void
t8_default_scheme_hex_c::t8_element_last_descendant (const t8_element_t *
                                                     elem,
                                                     t8_element_t * desc)
{
  p8est_quadrant_last_descendant ((p8est_quadrant_t *) elem,
                                  (p8est_quadrant_t *) desc, P8EST_QMAXLEVEL);
}

void
t8_default_scheme_hex_c::t8_element_successor (const t8_element_t * elem1,
                                               t8_element_t * elem2,
                                               int level)
{
  uint64_t            id;
  T8_ASSERT (0 <= level && level <= P8EST_QMAXLEVEL);

  id = p8est_quadrant_linear_id ((const p8est_quadrant_t *) elem1, level);
  T8_ASSERT (id + 1 < ((uint64_t) 1) << P8EST_DIM * level);
  p8est_quadrant_set_morton ((p8est_quadrant_t *) elem2, level, id + 1);
}

void
t8_default_scheme_hex_c::t8_element_anchor (const t8_element_t * elem,
                                            int coord[3])
{
  p8est_quadrant_t   *q;

  q = (p8est_quadrant_t *) elem;
  coord[0] = q->x;
  coord[1] = q->y;
  coord[2] = q->z;
}

int
t8_default_scheme_hex_c::t8_element_root_len (const t8_element_t * elem)
{
  return P8EST_ROOT_LEN;
}

/* Constructor */
t8_default_scheme_hex_c::t8_default_scheme_hex_c (void)
{
  eclass = T8_ECLASS_HEX;
  ts_context = sc_mempool_new (sizeof (t8_phex_t));
}

t8_default_scheme_hex_c::~t8_default_scheme_hex_c ()
{
  /* This destructor is empty since the destructor of the
   * default_common scheme is called automatically and it
   * suffices to destroy the quad_scheme.
   * However we need to provide an implementation of the destructor
   * and hence this empty function. */
}
