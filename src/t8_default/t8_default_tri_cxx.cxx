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
#include "t8_default_tri_cxx.hxx"
#include "t8_dtri_bits.h"

typedef t8_dtri_t   t8_default_tri_t;

int
t8_default_scheme_tri_c::t8_element_maxlevel (void)
{
  return T8_DTRI_MAXLEVEL;
}

int
t8_default_scheme_tri_c::t8_element_level (const t8_element_t * elem)
{
  return t8_dtri_get_level ((t8_dtri_t *) elem);
}

void
t8_default_scheme_tri_c::t8_element_copy (const t8_element_t * source,
                                          t8_element_t * dest)
{
  t8_dtri_copy ((const t8_dtri_t *) source, (t8_dtri_t *) dest);
}

int
t8_default_scheme_tri_c::t8_element_compare (const t8_element_t * elem1,
                                             const t8_element_t * elem2)
{
  int                 maxlvl;
  u_int64_t           id1, id2;

  /* Compute the bigger level of the two */
  maxlvl = SC_MAX (t8_element_level (elem1), t8_element_level (elem2));
  /* Compute the linear ids of the elements */
  id1 = t8_default_scheme_tri_c::t8_element_get_linear_id (elem1, maxlvl);
  id2 = t8_default_scheme_tri_c::t8_element_get_linear_id (elem2, maxlvl);
  /* return negativ if id1 < id2, zero if id1 = id2, positive if id1 > id2 */
  return id1 < id2 ? -1 : id1 != id2;
}

void
t8_default_scheme_tri_c::t8_element_parent (const t8_element_t * elem,
                                            t8_element_t * parent)
{
  const t8_default_tri_t *t = (const t8_default_tri_t *) elem;
  t8_default_tri_t   *p = (t8_default_tri_t *) parent;

  t8_dtri_parent (t, p);
}

void
t8_default_scheme_tri_c::t8_element_sibling (const t8_element_t * elem,
                                             int sibid,
                                             t8_element_t * sibling)
{
  const t8_default_tri_t *t = (const t8_default_tri_t *) elem;
  t8_default_tri_t   *s = (t8_default_tri_t *) sibling;

  t8_dtri_sibling (t, sibid, s);
}

void
t8_default_scheme_tri_c::t8_element_child (const t8_element_t * elem,
                                           int childid, t8_element_t * child)
{
  const t8_default_tri_t *t = (const t8_default_tri_t *) elem;
  t8_default_tri_t   *c = (t8_default_tri_t *) child;

  t8_dtri_child (t, childid, c);
}

void
t8_default_scheme_tri_c::t8_element_children (const t8_element_t * elem,
                                              int length, t8_element_t * c[])
{
  T8_ASSERT (length == T8_DTRI_CHILDREN);

  t8_dtri_childrenpv ((const t8_dtri_t *) elem, (t8_dtri_t **) c);
}

int
t8_default_scheme_tri_c::t8_element_child_id (const t8_element_t * elem)
{
  return t8_dtri_child_id ((t8_dtri_t *) elem);
}

int
t8_default_scheme_tri_c::t8_element_is_family (t8_element_t ** fam)
{
  return t8_dtri_is_familypv ((const t8_dtri_t **) fam);
}

void
t8_default_scheme_tri_c::t8_element_nca (const t8_element_t * elem1,
                                         const t8_element_t * elem2,
                                         t8_element_t * nca)
{
  const t8_default_tri_t *t1 = (const t8_default_tri_t *) elem1;
  const t8_default_tri_t *t2 = (const t8_default_tri_t *) elem1;
  t8_default_tri_t   *c = (t8_default_tri_t *) nca;

  t8_dtri_nearest_common_ancestor (t1, t2, c);
}

void
t8_default_scheme_tri_c::t8_element_set_linear_id (t8_element_t * elem,
                                                   int level, uint64_t id)
{
  T8_ASSERT (0 <= level && level <= T8_DTRI_MAXLEVEL);
  T8_ASSERT (0 <= id && id < ((uint64_t) 1) << (2 * level));

  t8_dtri_init_linear_id ((t8_default_tri_t *) elem, id, level);
}

uint64_t
  t8_default_scheme_tri_c::t8_element_get_linear_id (const t8_element_t *
                                                     elem, int level)
{
  T8_ASSERT (0 <= level && level <= T8_DTRI_MAXLEVEL);

  return t8_dtri_linear_id ((t8_default_tri_t *) elem, level);
}

void
t8_default_scheme_tri_c::t8_element_first_descendant (const t8_element_t *
                                                      elem,
                                                      t8_element_t * desc)
{
  t8_dtri_first_descendant ((t8_dtri_t *) elem, (t8_dtri_t *) desc);
}

void
t8_default_scheme_tri_c::t8_element_last_descendant (const t8_element_t *
                                                     elem,
                                                     t8_element_t * desc)
{
  t8_dtri_last_descendant ((t8_dtri_t *) elem, (t8_dtri_t *) desc);
}

void
t8_default_scheme_tri_c::t8_element_successor (const t8_element_t *
                                               elem1,
                                               t8_element_t * elem2,
                                               int level)
{
  T8_ASSERT (0 <= level && level <= T8_DTRI_MAXLEVEL);

  t8_dtri_successor ((const t8_default_tri_t *) elem1,
                     (t8_default_tri_t *) elem2, level);
}

void
t8_default_scheme_tri_c::t8_element_anchor (const t8_element_t * elem,
                                            int anchor[3])
{
  t8_dtri_t          *tri = (t8_dtri_t *) elem;

  anchor[0] = tri->x;
  anchor[1] = tri->y;
  anchor[2] = 0;
}

int
t8_default_scheme_tri_c::t8_element_root_len (const t8_element_t * elem)
{
  return T8_DTRI_ROOT_LEN;
}

/* Constructor */
t8_default_scheme_tri_c::t8_default_scheme_tri_c (void)
{
  eclass = T8_ECLASS_TRIANGLE;
  element_size = sizeof (t8_dtri_t);
  ts_context = sc_mempool_new (sizeof (element_size));
}

/* Destructor */
t8_default_scheme_tri_c::~t8_default_scheme_tri_c ()
{
  /* This destructor is empty since the destructor of the
   * default_common scheme is called automatically and it
   * suffices to destroy the quad_scheme.
   * However we need to provide an implementation of the destructor
   * and hence this empty function. */
}
