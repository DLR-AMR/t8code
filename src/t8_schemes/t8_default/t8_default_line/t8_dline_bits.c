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

#include <t8_schemes/t8_default/t8_default_line/t8_dline_bits.h>

int
t8_dline_get_level (const t8_dline_t *l)
{
  return l->level;
}

void
t8_dline_copy (const t8_dline_t *l, t8_dline_t *dest)
{
  memcpy (dest, l, sizeof (t8_dline_t));
}

int
t8_dline_compare (const t8_dline_t *l1, const t8_dline_t *l2)
{
  t8_linearidx_t id1, id2;

  /* Compute the linear ids of the elements */
  id1 = l1->x;
  id2 = l2->x;
  if (id1 == id2) {
    /* The linear ids are the same, the line with the smaller level
     * is considered smaller */
    return l1->level - l2->level;
  }
  /* return negative if id1 < id2, zero if id1 = id2, positive if id1 > id2 */
  return id1 < id2 ? -1 : id1 != id2;
}

int
t8_dline_equal (const t8_dline_t *l1, const t8_dline_t *l2)
{
  return (l1->level == l2->level && l1->x == l2->x);
}

void
t8_dline_parent (const t8_dline_t *l, t8_dline_t *parent)
{
  t8_dline_coord_t h;

  T8_ASSERT (l->level > 0);

  /* Get the length of l */
  h = T8_DLINE_LEN (l->level);

  /* Set coordinates of parent */
  parent->x = l->x & ~h;
  /* Set the parent's level */
  parent->level = l->level - 1;
}

void
t8_dline_ancestor (const t8_dline_t *l, int level, t8_dline_t *ancestor)
{
  /* Compute the new x-coordinate by setting all bits in positions
   * greater than level to zero. */
  ancestor->x = l->x & ~(T8_DLINE_LEN (level) - 1);
  /* set the level */
  ancestor->level = level;
}

void
t8_dline_child (const t8_dline_t *l, int childid, t8_dline_t *child)
{
  t8_dline_coord_t h;

  T8_ASSERT (l->level < T8_DLINE_MAXLEVEL);
  T8_ASSERT (childid == 0 || childid == 1);

  /* Compute the length of the child */
  h = T8_DLINE_LEN (l->level + 1);
  /* If childid = 0 then the children x coord is the same as l's,
   * if childid = 1 then it is x + h.
   */
  child->x = l->x + (childid == 0 ? 0 : h);
  /* The children level */
  child->level = l->level + 1;
}

void
t8_dline_face_neighbour (const t8_dline_t *l, t8_dline_t *neigh, int face, int *dual_face)
{
  T8_ASSERT (0 <= face && face < T8_DLINE_FACES);

  neigh->level = l->level;
  switch (face) {
  case 0:
    neigh->x = l->x - T8_DLINE_LEN (l->level);
    break;
  case 1:
    neigh->x = l->x + T8_DLINE_LEN (l->level);
    break;
  }
  if (dual_face != NULL) {
    /* The dual face is 1 if face=0 and 0 if face=1 */
    *dual_face = 1 - face;
  }
}

void
t8_dline_nearest_common_ancestor (const t8_dline_t *t1, const t8_dline_t *t2, t8_dline_t *r)
{
  int level;
  t8_dline_coord_t exclusive_or;

  /* At first we compute the first level at which the two x-coordinates differ */
  exclusive_or = t1->x ^ t2->x;
  level = SC_LOG2_32 (exclusive_or) + 1;

  T8_ASSERT (level <= T8_DLINE_MAXLEVEL);

  /* Compute the new x-coordinate by zeroing out all higher levels */
  r->x = t1->x & ~((1 << level) - 1);
  r->level = SC_MIN (T8_DLINE_MAXLEVEL - level, SC_MIN (t1->level, t2->level));
}

int
t8_dline_ancestor_id (const t8_dline_t *l, int level)
{
  t8_dline_coord_t h;

  T8_ASSERT (0 <= level && level <= T8_DLINE_MAXLEVEL);
  h = T8_DLINE_LEN (level);

  if (level == 0) {
    /* The root line as id 0 */
    return 0;
  }

  /* If the h bit is set in l's x coordinate, then acestor id is 1
   * else 0. */
  return (l->x & h) != 0;
}

int
t8_dline_face_parent_face (const t8_dline_t *l, int face)
{
  T8_ASSERT (0 <= face && face < T8_DLINE_FACES);

  if (l->level == 0) {
    return face;
  }
  /* If the child id is 0 and face is 0, then the parent's face is 0.
   * If the child id is 1 and face is 1, then the parent's face is 1.
   * In the other cases, this is an inner face. */
  return t8_dline_child_id (l) == face ? face : -1;
}

int
t8_dline_child_id (const t8_dline_t *elem)
{
  T8_ASSERT (elem->level < T8_DLINE_MAXLEVEL);
  /* bitshifting the Levelbit to first position & check if it is 1 or 0 */
  return ((elem->x >> (T8_DLINE_MAXLEVEL - elem->level)) & 1);
}

void
t8_dline_childrenpv (const t8_dline_t *elem, t8_dline_t *c[T8_DLINE_CHILDREN])
{
  const int8_t level = elem->level;

  T8_ASSERT (elem->level < T8_DLINE_MAXLEVEL);

  /* Set the Level, Level increases */
  c[0]->level = level + 1;
  c[1]->level = level + 1;
  /* Set the coordinates of the children */
  c[0]->x = elem->x;
  c[1]->x = elem->x + T8_DLINE_LEN (c[1]->level);
}

int
t8_dline_extrude_face (const t8_dvertex_t *face, int root_face, t8_dline_t *line)
{
  T8_ASSERT (root_face == 0 || root_face == 1);

  /* The level of the line is the same as the level of the boundary vertex */
  line->level = face->level;
  /* The x-coord of the line is either 0 (face = 0) or the length of the root
   * tree minus the length of the line */
  line->x = root_face == 0 ? 0 : T8_DLINE_ROOT_LEN - T8_DLINE_LEN (line->level);

  return root_face;
}

int
t8_dline_is_familypv (const t8_dline_t *f[])
{
  const int8_t level = f[0]->level;
  /* Check the level */
  if (level == 0 || level != f[1]->level) {
    return 0;
  } /* Check the parent */
  else if ((f[0]->x >> (T8_DLINE_MAXLEVEL - level + 1)) != (f[1]->x >> (T8_DLINE_MAXLEVEL - level + 1))) {
    return 0;
  }

  const t8_dline_coord_t len = T8_DLINE_LEN (level);
  /* Check the coordinate */
  return (f[0]->x + len == f[1]->x);
}

int
t8_dline_is_root_boundary (const t8_dline_t *p, int face)
{
  /* A line is at the boundary if and only if
   * face = 0 and x = 0
   * or
   * face = 1 and x = Maximum x - length of line
   */
  if (face == 0) {
    return p->x == 0;
  }
  else {
    return p->x == T8_DLINE_ROOT_LEN - T8_DLINE_LEN (p->level);
  }
}

int
t8_dline_is_inside_root (const t8_dline_t *l)
{
  return (l->x >= 0 && l->x < T8_DLINE_ROOT_LEN);
}

void
t8_dline_init_linear_id (t8_dline_t *l, int level, t8_linearidx_t id)
{
  T8_ASSERT (0 <= level && level <= T8_DLINE_MAXLEVEL);
  T8_ASSERT (0 <= id && id < ((t8_linearidx_t) 1) << level);

  /* Set the level */
  l->level = level;
  /* Set the new x coordinate */
  l->x = id << (T8_DLINE_MAXLEVEL - level);
}

void
t8_dline_successor (const t8_dline_t *l, t8_dline_t *succ, int level)
{
  t8_dline_coord_t h = 0;

  T8_ASSERT (1 <= level && level <= l->level);
  T8_ASSERT (l->x < T8_DLINE_ROOT_LEN - T8_DLINE_LEN (level));

  /* To compute the successor we zero out all bits in places bigger
   * than level and then we add the length of a line of level. */
  h = T8_DLINE_LEN (level) - 1;
  succ->x = l->x & ~h;
  succ->x += T8_DLINE_LEN (level);
  succ->level = level;
}

void
t8_dline_transform_face (const t8_dline_t *line1, t8_dline_t *line2, int orientation)
{
  T8_ASSERT (orientation == 0 || orientation == 1);

  /* The transformed line has the same level */
  line2->level = line1->level;
  if (orientation == 0) {
    /* If orientation is zero then the transformed line is the same
     * as the original line. */
    line2->x = line1->x;
  }
  else {
    t8_dline_coord_t h;
    T8_ASSERT (orientation == 1);
    /* Otherwise the lines are placed like this;
     *
     * 0 ---|_|------- N
     * N ---|_|------- 0
     *
     * With N = 2^maxlvl the root length, |_| marks the line elements within their trees.
     * Thus, the x-coordinate of line2 is given as N-line1.x - h.
     * Where h is the length of the line element.
     */
    h = T8_DLINE_LEN (line1->level);
    line2->x = T8_DLINE_ROOT_LEN - line1->x - h;
  }
}

void
t8_dline_first_descendant (const t8_dline_t *l, t8_dline_t *s, int level)
{
  T8_ASSERT (level >= l->level && level <= T8_DLINE_MAXLEVEL);

  s->level = level;
  s->x = l->x;
}

void
t8_dline_last_descendant (const t8_dline_t *l, t8_dline_t *s, int level)
{
  T8_ASSERT (level >= l->level && level <= T8_DLINE_MAXLEVEL);

  s->level = level;
  s->x = l->x + T8_DLINE_LEN (l->level) - T8_DLINE_LEN (level);
}

void
t8_dline_vertex_integer_coords (const t8_dline_t *elem, const int vertex, int coords[])
{
  T8_ASSERT (vertex == 0 || vertex == 1);
  if (vertex == 0) {
    coords[0] = elem->x;
  }
  else if (vertex == 1) {
    coords[0] = elem->x + T8_DLINE_LEN (elem->level);
  }
}

void
t8_dline_vertex_ref_coords (const t8_dline_t *elem, const int vertex, double coordinates[1])
{
  /* we need to set and initial value to prevent compiler warning. */
  int coords_int = -1;
  T8_ASSERT (vertex == 0 || vertex == 1);

  /* Compute integer coordinates and divide by root length. */
  t8_dline_vertex_integer_coords (elem, vertex, &coords_int);
  coordinates[0] = coords_int / (double) T8_DLINE_ROOT_LEN;
}

void
t8_dline_compute_reference_coords (const t8_dline_t *elem, const double *ref_coords, const size_t num_coords,
                                   const size_t skip_coords, double *out_coords)
{
  T8_ASSERT (t8_dline_is_valid (elem));
  for (size_t coord = 0; coord < num_coords; ++coord) {
    const size_t offset = coord * (1 + skip_coords);
    const size_t offset_3d = coord * 3;
    out_coords[offset] = elem->x;
    out_coords[offset] += T8_DLINE_LEN (elem->level) * ref_coords[offset_3d];
    out_coords[offset] /= (double) T8_DLINE_ROOT_LEN;
  }
}

t8_linearidx_t
t8_dline_linear_id (const t8_dline_t *elem, int level)
{
  t8_linearidx_t id;

  T8_ASSERT (level <= T8_DLINE_MAXLEVEL && level >= 0);

  /* this preserves the high bits from negative numbers */
  id = elem->x >> (T8_DLINE_MAXLEVEL - level);

  return id;
}

int
t8_dline_is_valid (const t8_dline_t *l)
{
  t8_dline_coord_t max_coord;

  max_coord = ((int64_t) 2 * T8_DLINE_ROOT_LEN) - 1;
  /* A line is valid if its level and its x coordinates are in the
   * correct bounds of the root three and its left and right neighbor */
  return 0 <= l->level && l->level <= T8_DLINE_MAXLEVEL && -T8_DLINE_ROOT_LEN <= l->x && l->x <= max_coord;
}

void
t8_dline_init (t8_dline_t *l)
{
  l->level = l->x = 0;
}
