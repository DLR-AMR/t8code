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

#include "t8_default_common_cxx.hxx"
#include "t8_default_line_cxx.hxx"
#include "t8_dline_bits.h"
#include "t8_dline.h"

typedef t8_dline_t  t8_default_line_t;

int
t8_default_scheme_line_c::t8_element_maxlevel (void)
{
  return T8_DLINE_MAXLEVEL;
}

int
t8_default_scheme_line_c::t8_element_level (const t8_element_t * elem)
{
  return t8_dline_get_level ((const t8_dline_t *) elem);
}

void
t8_default_scheme_line_c::t8_element_copy (const t8_element_t * source,
                                           t8_element_t * dest)
{
  t8_dline_copy ((const t8_dline_t *) source, (t8_dline_t *) dest);
}

int
t8_default_scheme_line_c::t8_element_compare (const t8_element_t * elem1,
                                              const t8_element_t * elem2)
{
  return t8_dline_compare ((const t8_dline_t *) elem1,
                           (const t8_dline_t *) elem2);
}

void
t8_default_scheme_line_c::t8_element_parent (const t8_element_t * elem,
                                             t8_element_t * parent)
{
  const t8_default_line_t *l = (const t8_default_line_t *) elem;
  t8_default_line_t  *p = (t8_default_line_t *) parent;

  t8_dline_parent (l, p);
}

void
t8_default_scheme_line_c::t8_element_child (const t8_element_t * elem,
                                            int childid, t8_element_t * child)
{
  const t8_default_line_t *l = (const t8_default_line_t *) elem;
  t8_default_line_t  *c = (t8_default_line_t *) child;

  t8_dline_child (l, childid, c);
}

t8_eclass_t
  t8_default_scheme_line_c::t8_element_face_class (const t8_element_t * elem,
                                                   int face)
{
  return T8_ECLASS_VERTEX;
}

void
t8_default_scheme_line_c::t8_element_transform_face (const t8_element_t *
                                                     elem1,
                                                     t8_element_t * elem2,
                                                     int orientation,
                                                     int is_smaller_face)
{
  T8_ASSERT (orientation == 0 || orientation == 1);

  /* We can ignore is_smaller_face, since for lines the orientation is independent
   * of the face. */
  t8_dline_transform_face ((const t8_dline_t *) elem1, (t8_dline_t *) elem2,
                           orientation);
}

void
t8_default_scheme_line_c::t8_element_set_linear_id (t8_element_t * elem,
                                                    int level, uint64_t id)
{
  T8_ASSERT (0 <= level && level <= T8_DLINE_MAXLEVEL);
  T8_ASSERT (0 <= id && id < ((u_int64_t) 1) << level);

  t8_dline_init_linear_id ((t8_default_line_t *) elem, level, id);
}

void
t8_default_scheme_line_c::t8_element_successor (const t8_element_t * elem1,
                                                t8_element_t * elem2,
                                                int level)
{
  T8_ASSERT (1 <= level && level <= T8_DLINE_MAXLEVEL);

  t8_dline_successor ((const t8_default_line_t *) elem1,
                      (t8_default_line_t *) elem2, level);
}

void
t8_default_scheme_line_c::t8_element_first_descendant (const t8_element_t *
                                                       elem,
                                                       t8_element_t * desc)
{
  t8_dline_first_descendant ((const t8_dline_t *) elem, (t8_dline_t *) desc,
                             T8_DLINE_MAXLEVEL);
}

void
t8_default_scheme_line_c::t8_element_last_descendant (const t8_element_t *
                                                      elem,
                                                      t8_element_t * desc)
{
  t8_dline_last_descendant ((const t8_dline_t *) elem, (t8_dline_t *) desc,
                            T8_DLINE_MAXLEVEL);
}

void
t8_default_scheme_line_c::t8_element_vertex_coords (const t8_element_t * t,
                                                    int vertex, int coords[])
{
  t8_dline_vertex_coords ((const t8_dline_t *) t, vertex, coords);
}

int
t8_default_scheme_line_c::t8_element_root_len (const t8_element_t * elem)
{
  return T8_DLINE_ROOT_LEN;
}

u_int64_t
  t8_default_scheme_line_c::t8_element_get_linear_id (const t8_element_t *
                                                      elem, int level)
{
  T8_ASSERT (0 <= level && level <= T8_DLINE_MAXLEVEL);

  return t8_dline_linear_id ((const t8_dline_t *) elem, level);
}

int
t8_default_scheme_line_c::t8_element_num_children (const t8_element_t * elem)
{
  return T8_DLINE_CHILDREN;
}

int
t8_default_scheme_line_c::t8_element_child_id (const t8_element_t * elem)
{
  return t8_dline_child_id ((const t8_dline_t *) elem);
}

void
t8_default_scheme_line_c::t8_element_children (const t8_element_t * elem,
                                               int length, t8_element_t * c[])
{
  T8_ASSERT (length == T8_DLINE_CHILDREN);

  t8_dline_childrenpv ((const t8_dline_t *) elem, (t8_dline_t **) c);
}

int
t8_default_scheme_line_c::t8_element_is_family (t8_element_t ** fam)
{
  return t8_dline_is_familypv ((const t8_dline_t **) fam);
}

/* Constructor */
t8_default_scheme_line_c::t8_default_scheme_line_c (void)
{
  eclass = T8_ECLASS_LINE;
  element_size = sizeof (t8_default_line_t);
  ts_context = sc_mempool_new (sizeof (t8_default_line_t));
}

t8_default_scheme_line_c::~t8_default_scheme_line_c ()
{
  /* This destructor is empty since the destructor of the
   * default_common scheme is called automatically and it
   * suffices to destroy the quad_scheme.
   * However we need to provide an implementation of the destructor
   * and hence this empty function. */
}
