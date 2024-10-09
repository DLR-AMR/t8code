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

#include <t8_schemes/t8_default/t8_default_vertex/t8_dvertex_bits.h>

int
t8_dvertex_get_level (const t8_dvertex_t *v)
{
  return v->level;
}

void
t8_dvertex_copy (const t8_dvertex_t *v, t8_dvertex_t *dest)
{
  memcpy (dest, v, sizeof (t8_dvertex_t));
}

int
t8_dvertex_equal (const t8_dvertex_t *elem1, const t8_dvertex_t *elem2)
{
  return elem1->level == elem2->level;
}

int
t8_dvertex_compare (const t8_dvertex_t *v1, const t8_dvertex_t *v2)
{
  /* The vertex with the smaller level
   * is considered smaller */
  return v1->level - v2->level;
}

void
t8_dvertex_parent (const t8_dvertex_t *v, t8_dvertex_t *parent)
{
  T8_ASSERT (v->level > 0);

  /* Set the parent's level */
  parent->level = v->level - 1;
}

void
t8_dvertex_child (const t8_dvertex_t *v, t8_dvertex_t *child)
{
  T8_ASSERT (v->level < T8_DVERTEX_MAXLEVEL);

  /* The children level */
  child->level = v->level + 1;
}

void
t8_dvertex_nearest_common_ancestor (const t8_dvertex_t *v1, const t8_dvertex_t *v2, t8_dvertex_t *r)
{
  /* The nca is the one of the two vertices with smaller level */
  r->level = SC_MIN (v1->level, v2->level);
}

int
t8_dvertex_ancestor_id (const t8_dvertex_t *v, int level)
{
  /* There is only one possible child id, 0 */
  return 0;
}

int
t8_dvertex_child_id (const t8_dvertex_t *v)
{
  /* There is only one possible child id, 0 */
  return 0;
}

void
t8_dvertex_sibling (const t8_dvertex_t *v, int sibid, t8_dvertex_t *s)
{
  T8_ASSERT (sibid == 0);

  t8_dvertex_copy (v, s);
}

void
t8_dvertex_childrenpv (const t8_dvertex_t *v, t8_dvertex_t *c[T8_DVERTEX_CHILDREN])
{
  T8_ASSERT (v->level < T8_DVERTEX_MAXLEVEL);

  /* Set the Level, Level increases */
  c[0]->level = v->level + 1;
}

int
t8_dvertex_is_familypv (const t8_dvertex_t *f[])
{
  /* A vertex with level greater 0 is always a family */
  return f[0]->level > 0;
}

int
t8_dvertex_is_root_boundary (const t8_dvertex_t *v, int face)
{
  /* A vertex is always at the root boundary */
  return 1;
}

int
t8_dvertex_is_inside_root (const t8_dvertex_t *v)
{
  /* A vertex is always inside root */
  return 1;
}

void
t8_dvertex_init_linear_id (t8_dvertex_t *v, int level, t8_linearidx_t id)
{
  T8_ASSERT (0 <= level && level <= T8_DVERTEX_MAXLEVEL);
  T8_ASSERT (0 == id);

  /* Set the level */
  v->level = level;
}

void
t8_dvertex_transform_face (const t8_dvertex_t *vertex1, t8_dvertex_t *vertex2)
{
  /* The transformed vertex is the same vertex */
  vertex2->level = vertex1->level;
}

void
t8_dvertex_first_descendant (const t8_dvertex_t *v, t8_dvertex_t *s, int level)
{
  T8_ASSERT (level >= v->level && level <= T8_DVERTEX_MAXLEVEL);

  s->level = level;
}

void
t8_dvertex_last_descendant (const t8_dvertex_t *v, t8_dvertex_t *s, int level)
{
  T8_ASSERT (level >= v->level && level <= T8_DVERTEX_MAXLEVEL);

  s->level = level;
}

void
t8_dvertex_vertex_integer_coords (const t8_dvertex_t *elem, const int vertex, int coords[])
{
  T8_ASSERT (vertex == 0);

  coords[0] = 0;
}

void
t8_dvertex_vertex_ref_coords (const t8_dvertex_t *elem, const int vertex, double coords[])
{
  T8_ASSERT (vertex == 0);

  coords[0] = 0;
}

void
t8_dvertex_compute_reference_coords (const t8_dvertex_t *elem, const double *ref_coords, const size_t num_coords,
                                     const size_t padding, double *out_coords)
{
  T8_ASSERT (fabs (ref_coords[0]) <= T8_PRECISION_EPS);
  T8_ASSERT (t8_dvertex_is_valid (elem));
  for (size_t coord = 0; coord < num_coords; ++coord) {
    const size_t offset_out = coord * (1 + padding);
    out_coords[offset_out] = 0;
  }
}

t8_linearidx_t
t8_dvertex_linear_id (const t8_dvertex_t *elem, int level)
{
  T8_ASSERT (level <= T8_DVERTEX_MAXLEVEL && level >= 0);

  return 0;
}

int
t8_dvertex_is_valid (const t8_dvertex_t *v)
{
  /* A vertex is valid if its level is in the valid range */
  return 0 <= v->level && v->level <= T8_DVERTEX_MAXLEVEL;
}

void
t8_dvertex_debug_print (const t8_dvertex_t *v)
{
  t8_debugf ("level: %i\n", v->level);
}

void
t8_dvertex_init (t8_dvertex_t *v)
{
  v->level = 0;
}
