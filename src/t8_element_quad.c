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
t8_escheme_quad_destroy (void * context)
{
  sc_mempool_destroy ((sc_mempool_t *) context);
}

t8_escheme_t *
t8_escheme_new_quad (void)
{
  t8_escheme_t * es;

  es = T8_ALLOC (t8_escheme_t, 1);

  es->ttype = T8_TYPE_QUAD;
  es->elem_parent = (t8_element_parent_t) p4est_quadrant_parent;
  es->elem_sibling = t8_element_quad_sibling;
  es->elem_nca = (t8_element_nca_t) p4est_nearest_common_ancestor;

  es->escheme_destroy = t8_escheme_quad_destroy;
  es->context = sc_mempool_new (sizeof (p4est_quadrant_t));

  return es;
}
