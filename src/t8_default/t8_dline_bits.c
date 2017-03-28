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
  t8_dline_coord_t    h;
  int                 i;

  T8_ASSERT (1 <= level && level <= l->level);

  /* To compute the successor we zero out all bits in places bigger
   * than level and then we add the length of a line of level. */
  for (i = level + 1; i <= l->level; i++) {
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
