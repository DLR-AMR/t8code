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

#include "t8_default_common.h"
#include "t8_default_line.h"
#include "t8_dline_bits.h"
#include "t8_dline.h"

typedef t8_dline_t  t8_default_line_t;

static size_t
t8_default_line_size (void)
{
  return sizeof (t8_default_line_t);
}

static int
t8_default_line_maxlevel (void)
{
  return T8_DLINE_MAXLEVEL;
}

static int
t8_default_line_level (const t8_element_t * elem)
{
  return t8_dline_get_level ((t8_dline_t *) elem);
}

static void
t8_default_tline_copy (const t8_element_t * source, t8_element_t * dest)
{
  t8_dline_copy ((const t8_dline_t *) source, (t8_dline_t *) dest);
}

static void
t8_default_line_parent (const t8_element_t * elem, t8_element_t * parent)
{
  const t8_default_line_t *l = (const t8_default_line_t *) elem;
  t8_default_line_t  *p = (t8_default_line_t *) parent;

  t8_dline_parent (l, p);
}

t8_eclass_scheme_t *
t8_default_scheme_new_line (void)
{
  t8_eclass_scheme_t *ts;

  ts = T8_ALLOC_ZERO (t8_eclass_scheme_t, 1);
  ts->eclass = T8_ECLASS_LINE;

  ts->elem_size = t8_default_line_size;
  ts->elem_maxlevel = t8_default_line_maxlevel;

  ts->elem_level = t8_default_line_level;
  ts->elem_copy = t8_default_tline_copy;
  ts->elem_compare = NULL;
  ts->elem_parent = t8_default_line_parent;
  ts->elem_sibling = NULL;
  ts->elem_child = NULL;
  ts->elem_children = NULL;
  ts->elem_is_family = NULL;
  ts->elem_child_id = NULL;
  ts->elem_nca = NULL;
  ts->elem_set_linear_id = NULL;
  ts->elem_get_linear_id = NULL;
  ts->elem_first_desc = NULL;
  ts->elem_last_desc = NULL;
  ts->elem_successor = NULL;
  ts->elem_anchor = NULL;
  ts->elem_root_len = NULL;

  ts->elem_new = t8_default_mempool_alloc;
  ts->elem_destroy = t8_default_mempool_free;

  ts->ts_destroy = t8_default_scheme_mempool_destroy;
  ts->ts_context = sc_mempool_new (sizeof (t8_default_line_t));

  return ts;
}
