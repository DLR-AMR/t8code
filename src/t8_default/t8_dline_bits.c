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

#include "t8_dline_bits.h"

int
t8_dline_get_level (const t8_dline_t * l)
{
  return l->level;
}

void
t8_dline_copy (const t8_dline_t * l, t8_dline_t * dest)
{
  memcpy (dest, l, sizeof (t8_dline_t));
}

int
t8_dline_compare (const t8_dline_t * l1, const t8_dline_t * l2)
{
    int                 maxlvl;
    u_int64_t           id1, id2;

    maxlvl = SC_MAX (l1->level, l2->level);
    /* Compute the linear ids of the elements */
    id1 = t8_dline_linear_id (l1, maxlvl);
    id2 = t8_dline_linear_id (l2, maxlvl);
    if (id1 == id2) {
    /* The linear ids are the same, the triangle with the smaller level
    * is considered smaller */
    return l1->level - l2->level;
    }
    /* return negativ if id1 < id2, zero if id1 = id2, positive if id1 >
    id2 */
    return id1 < id2 ? -1 : id1 != id2;
}

void
t8_dline_parent (const t8_dline_t * l, t8_dline_t * parent)
{
  t8_dline_coord_t    h;

  T8_ASSERT (l->level > 0);

  /* Get the length of l */
  h = T8_DLINE_LEN (l->level);

  /* Set coordinates of parent */
  parent->x = l->x & ~h;
  /* Set the parent's level */
  parent->level = l->level - 1;
}

void
t8_dline_child (const t8_dline_t * l, int childid, t8_dline_t * child)
{
  t8_dline_coord_t    h;

  T8_ASSERT (l->level < T8_DLINE_MAXLEVEL);
  T8_ASSERT (childid == 0 || childid == 1);

  /* Compute the length of the child */
  h = T8_DLINE_LEN (l->level + 1);
  /* If childid = 0 then the childs x coord is the same as l's,
   * if childid = 1 then it is x + h.
   */
  child->x = l->x + (childid == 0 ? 0 : h);
  /* The childs level */
  child->level = l->level + 1;
}

int
t8_dline_child_id (const t8_dline_t * elem)
{
  T8_ASSERT (elem->level < T8_DLINE_MAXLEVEL);
  /* bitshifting the Levelbit to first position & check if it is 1 or 0 */
  return ((elem->x >> (T8_DLINE_MAXLEVEL - elem->level)) & 1);
}

void
t8_dline_childrenpv (const t8_dline_t * elem,
                     t8_dline_t * c[T8_DLINE_CHILDREN])
{
  const int8_t        level = elem->level;

  T8_ASSERT (elem->level < T8_DLINE_MAXLEVEL);

  /* Set the Level, Level increases */
  c[0]->level = level + 1;
  c[1]->level = level + 1;
  /* Set the coordinates of the children */
  c[0]->x = elem->x;
  c[1]->x = elem->x + T8_DLINE_LEN (c[1]->level);
}

int
t8_dline_is_familypv (const t8_dline_t * f[])
{
  const int8_t        level = f[0]->level;
  t8_dline_coord_t    len = T8_DLINE_LEN (level);

  /*Check the level */
  if (level == 0 || level != f[1]->level) {
    return 0;
  }                             /* Check the parent */
  else if ((f[0]->x >> (T8_DLINE_MAXLEVEL - level + 1)) !=
           (f[1]->x >> (T8_DLINE_MAXLEVEL - level + 1))) {
    return 0;
  }

  /*Check the coordinate */
  return (f[0]->x + len == f[1]->x);
}

void
t8_dline_init_linear_id (t8_dline_t * l, int level, uint64_t id)
{
  T8_ASSERT (0 <= level && level <= T8_DLINE_MAXLEVEL);
  T8_ASSERT (0 <= id && id <= ((uint64_t) 1) << level);

  /* Set the level */
  l->level = level;
  /* Set the new x coordinate */
  l->x = id << (T8_DLINE_MAXLEVEL - level);
}

void
t8_dline_successor (const t8_dline_t * l, t8_dline_t * succ, int level)
{
  t8_dline_coord_t    h = 0;
  int                 i;

  T8_ASSERT (1 <= level && level <= l->level);

  /* To compute the successor we zero out all bits in places bigger
   * than level and then we add the length of a line of level. */
  for (i = level + 1; i <= T8_DLINE_MAXLEVEL; i++) {
    h |= T8_DLINE_LEN (i);
  }
  succ->x = l->x & ~h;
  succ->x += T8_DLINE_LEN (level);
  succ->level = level;
}

void
t8_dline_first_descendant (const t8_dline_t * l, t8_dline_t * s, int level)
{
  T8_ASSERT (level >= l->level && level <= T8_DLINE_MAXLEVEL);

  s->level = level;
  s->x = l->x;
}

void
t8_dline_last_descendant (const t8_dline_t * l, t8_dline_t * s, int level)
{
  T8_ASSERT (level >= l->level && level <= T8_DLINE_MAXLEVEL);

  s->level = level;
  s->x = l->x + T8_DLINE_LEN (l->level) - T8_DLINE_LEN (level);
}

void
t8_dline_vertex_coords (const t8_dline_t * elem, int vertex, int coords[])
{
  T8_ASSERT (vertex == 0 || vertex == 1);
  if (vertex == 0) {
    coords[0] = elem->x;
  }
  if (vertex == 1) {
    coords[0] = elem->x + T8_DLINE_LEN (elem->level);
  }
}

uint64_t
t8_dline_linear_id (const t8_dline_t * elem, int level)
{
  uint64_t            id;

  T8_ASSERT (level <= T8_DLINE_MAXLEVEL && level >= 0);

  /* this preserves the high bits from negative numbers */
  id = elem->x >> (T8_DLINE_MAXLEVEL - level);

  return id;
}
