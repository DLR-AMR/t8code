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

#include "t8_default_tet_bits.h"
#include "t8_default_common.h"
#include "t8_default_tet.h"

static              size_t
t8_default_tet_size (void)
{
  return sizeof (t8_tet_t);
}

t8_default_tet_type_t
t8_default_tet_get_type (const t8_tet_t * t)
{
  return t->type;
}

void
t8_default_tet_set_type (t8_tet_t * t, t8_default_tet_type_t type)
{
  T8_ASSERT (0 <= type && type <= 5);

  t->type = type;
}

t8_tcoord_t
t8_default_tet_get_coordinate (const t8_tet_t * t, int i)
{
  T8_ASSERT (0 <= i && i < 3);
  return t->anchor_coordinates[i];
}

void
t8_default_tet_set_coordinate (t8_tet_t * t, int i, t8_tcoord_t value)
{
  T8_ASSERT (0 <= i && i < 3);
  t->anchor_coordinates[i] = value;
}



t8_eclass_scheme_t *
t8_default_scheme_new_tet (void)
{
  t8_eclass_scheme_t *ts;

  ts = T8_ALLOC (t8_eclass_scheme_t, 1);

  ts->elem_size = t8_default_tet_size;
  ts->elem_parent = t8_default_tet_parent;
  ts->elem_sibling = t8_default_tet_sibling;
  ts->elem_child = t8_default_tet_child;

  ts->elem_new = t8_default_mempool_alloc;
  ts->elem_destroy = t8_default_mempool_free;
  ts->ts_destroy = t8_default_scheme_mempool_destroy;
  ts->ts_context = sc_mempool_new (sizeof (t8_tet_t));

  return ts;
}
