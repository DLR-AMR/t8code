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
#include "t8_default_pyramid_cxx.hxx"
#include "t8_dpyramid_bits.h"
#include "t8_dpyramid.h"

typedef t8_dpyramid_t t8_default_pyramid_t;

T8_EXTERN_C_BEGIN ();

int
t8_default_scheme_pyramid_c::t8_element_maxlevel(void)
{
    return T8_DPYRAMID_MAXLEVEL;
}

void
t8_default_scheme_pyramid_c::t8_element_init(int length, t8_element_t *elem,
                                             int called_new)
{
#ifdef T8_ENABLE_DEBUG
  if (!called_new) {
    int                 i;
    t8_dpyramid_t        *pyramid = (t8_dpyramid_t *) elem;
    /* Set all values to 0 */
    for (i = 0; i < length; i++) {
      t8_dpyramid_init_linear_id (pyramid + i, 0, 0);
    }
  }
#endif
}

int
t8_default_scheme_pyramid_c::t8_element_num_vertices (const t8_element_t * elem)
{
    return t8_dpyramid_num_vertices((const t8_dpyramid_t *) elem);
}

int
t8_default_scheme_pyramid_c::t8_element_compare (const t8_element_t * elem1,
                                                 const t8_element_t * elem2)
{
    t8_dpyramid_compare((const t8_dpyramid_t *) elem1, (const t8_dpyramid_t *) elem2 );
}

void
t8_default_scheme_pyramid_c::t8_element_child (const t8_element_t * elem,
                                        int childid, t8_element_t * child)
{
    t8_dpyramid_child((t8_dpyramid_t *)elem, childid, (t8_dpyramid_t *) child);
}

void
t8_default_scheme_pyramid_c::t8_element_first_descendant(const t8_element_t *elem,
                                                         t8_element_t *desc, int level)
{
    t8_dpyramid_first_descendant((const t8_dpyramid_t *) elem, (t8_dpyramid_t *) desc, level );
}

int
t8_default_scheme_pyramid_c::t8_element_level (const t8_element_t * elem)
{
    t8_dpyramid_get_level((const t8_dpyramid_t *) elem);
}

int
t8_default_scheme_pyramid_c::t8_element_root_len (const t8_element_t * elem)
{
    return T8_DPYRAMID_ROOT_LEN;
}

void
t8_default_scheme_pyramid_c::t8_element_set_linear_id(t8_element_t *elem, int level, uint64_t id)
{
    t8_dpyramid_init_linear_id((t8_dpyramid_t *) elem, level, id);
}

t8_eclass_t
t8_default_scheme_pyramid_c::t8_element_shape(const t8_element_t *elem)
{
    return t8_dpyramid_shape((const t8_dpyramid_t *) elem);
}

u_int64_t
t8_default_scheme_pyramid_c::t8_element_get_linear_id (const t8_element_t * elem, int level)
{
    t8_dpyramid_linear_id((const t8_dpyramid_t *) elem, level);
}

void
t8_default_scheme_pyramid_c::t8_element_last_descendant(const t8_element_t *elem,
                                                        t8_element_t *desc, int level)
{
    t8_dpyramid_last_descendant((const t8_dpyramid_t *) elem, (t8_dpyramid_t *) desc, level);
}

void
t8_default_scheme_pyramid_c::t8_element_successor (const t8_element_t * elem,
                                                   t8_element_t * s, int level)
{
    t8_dpyramid_succesor((const t8_dpyramid_t *) elem, (t8_dpyramid_t *) s, level);
}

void
t8_default_scheme_pyramid_c::t8_element_vertex_coords (const t8_element_t * t,
                                                int vertex, int coords[])
{
    t8_dpyramid_compute_coords((const t8_dpyramid_t *) t, vertex, coords);
}

/* Constructor */
t8_default_scheme_pyramid_c::t8_default_scheme_pyramid_c (void)
{
  eclass = T8_ECLASS_PYRAMID;
  element_size = sizeof (t8_default_pyramid_t);
  ts_context = sc_mempool_new (sizeof (element_size));
}

t8_default_scheme_pyramid_c::~t8_default_scheme_pyramid_c ()
{
  /* This destructor is empty since the destructor of the
   * default_common scheme is called automatically and it
   * suffices to destroy the quad_scheme.
   * However we need to provide an implementation of the destructor
   * and hence this empty function. */
}

T8_EXTERN_C_END ();
