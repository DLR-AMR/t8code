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
#include "t8_default_tet_cxx.hxx"
#include "t8_dtet_bits.h"

/* We want to export the whole implementation to be callable from "C" */
T8_EXTERN_C_BEGIN ();

typedef t8_dtet_t   t8_default_tet_t;

int
t8_default_scheme_tet_c::t8_element_maxlevel (void)
{
  return T8_DTET_MAXLEVEL;
}

int
t8_default_scheme_tet_c::t8_element_level (const t8_element_t * elem)
{
  return t8_dtet_get_level ((t8_dtet_t *) elem);
}

void
t8_default_scheme_tet_c::t8_element_copy (const t8_element_t * source,
                                          t8_element_t * dest)
{
  t8_dtet_copy ((const t8_dtet_t *) source, (t8_dtet_t *) dest);
}

int
t8_default_scheme_tet_c::t8_element_compare (const t8_element_t * elem1,
                                             const t8_element_t * elem2)
{
  int                 maxlvl;
  u_int64_t           id1, id2;

  /* Compute the bigger level of the two */
  maxlvl = SC_MAX (t8_default_scheme_tet_c::t8_element_level (elem1),
                   t8_default_scheme_tet_c::t8_element_level (elem2));
  /* Compute the linear ids of the elements */
  id1 = t8_default_scheme_tet_c::t8_element_get_linear_id (elem1, maxlvl);
  id2 = t8_default_scheme_tet_c::t8_element_get_linear_id (elem2, maxlvl);
  /* return negativ if id1 < id2, zero if id1 = id2, positive if id1 > id2 */
  return id1 < id2 ? -1 : id1 != id2;
}

void
t8_default_scheme_tet_c::t8_element_parent (const t8_element_t * elem,
                                            t8_element_t * parent)
{
  const t8_default_tet_t *t = (const t8_default_tet_t *) elem;
  t8_default_tet_t   *p = (t8_default_tet_t *) parent;

  t8_dtet_parent (t, p);
}

void
t8_default_scheme_tet_c::t8_element_sibling (const t8_element_t * elem,
                                             int sibid,
                                             t8_element_t * sibling)
{
  const t8_default_tet_t *t = (const t8_default_tet_t *) elem;
  t8_default_tet_t   *s = (t8_default_tet_t *) sibling;

  t8_dtet_sibling (t, sibid, s);
}

void
t8_default_scheme_tet_c::t8_element_child (const t8_element_t * elem,
                                           int childid, t8_element_t * child)
{
  const t8_default_tet_t *t = (const t8_default_tet_t *) elem;
  t8_default_tet_t   *c = (t8_default_tet_t *) child;

  t8_dtet_child (t, childid, c);
}

void
t8_default_scheme_tet_c::t8_element_children (const t8_element_t * elem,
                                              int length, t8_element_t * c[])
{
  T8_ASSERT (length == T8_DTET_CHILDREN);

  t8_dtet_childrenpv ((const t8_dtet_t *) elem, (t8_dtet_t **) c);
}

int
t8_default_scheme_tet_c::t8_element_child_id (const t8_element_t * elem)
{
  return t8_dtet_child_id ((t8_dtet_t *) elem);
}

int
t8_default_scheme_tet_c::t8_element_is_family (t8_element_t ** fam)
{
  return t8_dtet_is_familypv ((const t8_dtet_t **) fam);
}

void
t8_default_scheme_tet_c::t8_element_nca (const t8_element_t * elem1,
                                         const t8_element_t * elem2,
                                         t8_element_t * nca)
{
  const t8_default_tet_t *t1 = (const t8_default_tet_t *) elem1;
  const t8_default_tet_t *t2 = (const t8_default_tet_t *) elem2;
  t8_default_tet_t   *c = (t8_default_tet_t *) nca;

  t8_dtet_nearest_common_ancestor (t1, t2, c);
}

void
t8_default_scheme_tet_c::t8_element_set_linear_id (t8_element_t * elem,
                                                   int level, uint64_t id)
{
  T8_ASSERT (0 <= level && level <= T8_DTET_MAXLEVEL);
  T8_ASSERT (0 <= id && id < ((uint64_t) 1) << 3 * level);

  t8_dtet_init_linear_id ((t8_default_tet_t *) elem, id, level);
}

uint64_t
  t8_default_scheme_tet_c::t8_element_get_linear_id (const t8_element_t *
                                                     elem, int level)
{
  T8_ASSERT (0 <= level && level <= T8_DTET_MAXLEVEL);

  return t8_dtet_linear_id ((t8_default_tet_t *) elem, level);
}

void
t8_default_scheme_tet_c::t8_element_successor (const t8_element_t * elem1,
                                               t8_element_t * elem2,
                                               int level)
{
  T8_ASSERT (0 <= level && level <= T8_DTET_MAXLEVEL);

  t8_dtet_successor ((const t8_default_tet_t *) elem1,
                     (t8_default_tet_t *) elem2, level);
}

void
t8_default_scheme_tet_c::t8_element_first_descendant (const t8_element_t *
                                                      elem,
                                                      t8_element_t * desc)
{
  t8_dtet_first_descendant ((t8_dtet_t *) elem, (t8_dtet_t *) desc);
}

void
t8_default_scheme_tet_c::t8_element_last_descendant (const t8_element_t *
                                                     elem,
                                                     t8_element_t * desc)
{
  t8_dtet_last_descendant ((t8_dtet_t *) elem, (t8_dtet_t *) desc);
}

void
t8_default_scheme_tet_c::t8_element_anchor (const t8_element_t * elem,
                                            int anchor[3])
{
  t8_dtet_t          *tet = (t8_dtet_t *) elem;

  anchor[0] = tet->x;
  anchor[1] = tet->y;
  anchor[2] = tet->z;
}

int
t8_default_scheme_tet_c::t8_element_root_len (const t8_element_t * elem)
{
  return T8_DTET_ROOT_LEN;
}

void
t8_default_scheme_tet_c::t8_element_vertex_coords (const t8_element_t * t,
                                                   int vertex, int coords[])
{
  t8_dtet_compute_coords ((const t8_default_tet_t *) t, vertex, coords);
}

 /* Constructor */
t8_default_scheme_tet_c::t8_default_scheme_tet_c (void)
{
  eclass = T8_ECLASS_TET;
  element_size = sizeof (t8_dtet_t);
  ts_context = sc_mempool_new (sizeof (element_size));
}

 /* Destructor */
t8_default_scheme_tet_c::~t8_default_scheme_tet_c ()
{
  /* This destructor is empty since the destructor of the
   * default_common scheme is called automatically and it
   * suffices to destroy the quad_scheme.
   * However we need to provide an implementation of the destructor
   * and hence this empty function. */
}

T8_EXTERN_C_END ();
