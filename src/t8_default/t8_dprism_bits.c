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

#include "t8_dprism_bits.h"
#include "t8_dline_bits.h"
#include "t8_dtri_bits.h"

int
t8_dprism_get_level (const t8_dprism_t * p)
{
  T8_ASSERT (p->line.level == p->tri.level);
  return p->line.level;
}

void
t8_dprism_copy (const t8_dprism_t * p, t8_dprism_t * dest)
{
  T8_ASSERT (p->line.level == p->tri.level);
  memcpy (dest, p, sizeof (t8_dprism_t));
}

void
t8_dprism_init_linear_id (t8_dprism_t * p, int level, uint64_t id)
{
  uint64_t            tri_id = 0;
  uint64_t            line_id = 0;
  int                 i;
  int                 triangles_of_size_i = 1;

  T8_ASSERT (0 <= level && level <= T8_DPRISM_MAXLEVEL);
  T8_ASSERT (id < 1 << 3 * level);

  for (i = 0; i <= level; i++) {
    /*Get the number of the i-th prism and get the related triangle number
     * then multiplicate it by the number of triangles of level size.*/
    tri_id += ((id % 8) % 4) * triangles_of_size_i;

    /*If id % 8 is larger than 3, the prism is in the upper part of the
     * parent prism. => line_id + 2^i*/
    line_id += (id % 8 > 3) ? 1 << i : 0;

    /*Each Prism divides into 8 children */
    id /= 8;
    /*Each triangle divides into 4 children */
    triangles_of_size_i *= 4;
  }
  t8_dtri_init_linear_id (&p->tri, tri_id, level);
  t8_dline_init_linear_id (&p->line, level, line_id);
  T8_ASSERT (p->line.level == p->tri.level);
}

void
t8_dprism_parent (const t8_dprism_t * p, t8_dprism_t * parent)
{
  T8_ASSERT (p->line.level > 0);
  T8_ASSERT (p->line.level == p->tri.level);
  t8_dtri_parent (&p->tri, &parent->tri);
  t8_dline_parent (&p->line, &parent->line);
}

int
t8_dprism_child_id (const t8_dprism_t * p)
{
  int                 tri_child_id = t8_dtri_child_id (&p->tri);
  int                 line_child_id = t8_dline_child_id (&p->line);
  T8_ASSERT (p->line.level == p->tri.level);
  /*Prism in lower plane has the same id as the triangle, in the upper plane
   * it's a shift by 4*/
  return tri_child_id + 4 * line_child_id;
}

int
t8_dprism_is_familypv (t8_dprism_t ** fam)
{
  int                 i;
  t8_dtri_t         **tri_fam = malloc (4 * sizeof (t8_dtri_t *));
  t8_dline_t        **line_fam = malloc (2 * sizeof (t8_dline_t *));
  int                 is_family = 1;
  for (i = 0; i < 2; i++) {
    tri_fam[0] = &fam[0 + i * 4]->tri;
    tri_fam[1] = &fam[1 + i * 4]->tri;
    tri_fam[2] = &fam[2 + i * 4]->tri;
    tri_fam[3] = &fam[3 + i * 4]->tri;
    is_family = is_family && t8_dtri_is_familypv ((const t8_dtri_t **) tri_fam);
  }
  for (i = 0; i < 4; i++) {
    line_fam[0] = &fam[i]->line;
    line_fam[1] = &fam[i + 4]->line;
    is_family =  is_family && t8_dline_is_familypv ((const t8_dline_t **) line_fam);
  }
  free (tri_fam);
  free (line_fam);
  return is_family;
}

void
t8_dprism_child (const t8_dprism_t * p, int childid, t8_dprism_t * child)
{
  T8_ASSERT (0 <= childid && childid < 8);
  T8_ASSERT (p->line.level == p->tri.level);
  t8_dtri_child (&p->tri, childid % 4, &child->tri);
  t8_dline_child (&p->line, childid / 4, &child->line);
}

void
t8_dprism_childrenpv (const t8_dprism_t * p, int length, t8_dprism_t * c[])
{
  int                 i;
  T8_ASSERT (length == T8_DPRISM_CHILDREN);
  T8_ASSERT (p->line.level < T8_DPRISM_MAXLEVEL &&
             p->tri.level == p->line.level);
  for (i = 7; i >= 0; i--) {
    t8_dprism_child (p, i, c[i]);
  }
}

void
t8_dprism_successor (const t8_dprism_t * p, t8_dprism_t * succ, int level)
{
  int                 prism_child_id;
  t8_dprism_copy (p, succ);
  /*update the level */
  succ->line.level = level;
  succ->tri.level = level;
  prism_child_id = t8_dprism_child_id(succ);

  T8_ASSERT (1 <= level && level <= T8_DPRISM_MAXLEVEL);

  /*The next prism is the child with local ID 0 of the next parent prism */
  if (prism_child_id == 7) {
    t8_dprism_successor (p, succ, level - 1);
    /*Zero out the bits of higher level, caused by recursion */
#if 1
    succ->tri.x =
      (succ->tri.x >> (T8_DTRI_MAXLEVEL - level + 1)) << (T8_DTRI_MAXLEVEL -
                                                          level + 1);
    succ->tri.y =
      (succ->tri.y >> (T8_DTRI_MAXLEVEL - level + 1)) << (T8_DTRI_MAXLEVEL -
                                                          level + 1);
#endif
    succ->line.x =
      (succ->line.x >> (T8_DTRI_MAXLEVEL - level + 1)) << (T8_DTRI_MAXLEVEL -
                                                           level + 1);
    /*Set the level to the actual level */
    succ->line.level = level;
    succ->tri.level = level;
  }
  /*The next prism is one plane up, local_tri_id = 0 */
  else if (prism_child_id == 3) {
    /*parent is computed with succ, cause there are the updated datas */
    t8_dprism_parent(succ, succ);
    t8_dprism_child(succ, 4, succ);
  }
  /*The next Prism is in the same plane, but has the next base-triangle */
  else {
    t8_dtri_successor (&p->tri, &succ->tri, level);
  }

}

void
t8_dprism_first_descendant (const t8_dprism_t * p, t8_dprism_t * s, int level)
{
  T8_ASSERT (level >= p->line.level && level <= T8_DPRISM_MAXLEVEL);
  T8_ASSERT (p->line.level == p->tri.level);
  /*First prism descendant = first triangle desc x first line desc */

  t8_dtri_first_descendant (&l->tri, &s->tri, level);
  t8_dline_first_descendant (&l->line, &s->line, level);
}

void
t8_dprism_last_descendant (const t8_dprism_t * p, t8_dprism_t * s, int level)
{
  T8_ASSERT (level >= p->line.level && level <= T8_DPRISM_MAXLEVEL);
  /*Last prism descendant = last triangle desc x last line desc */
  T8_ASSERT (level == T8_DTRI_MAXLEVEL);
<<<<<<< HEAD
  T8_ASSERT (p->line.level == p->tri.level);
  /*TODO: MERGE mit johannes und dann level zufügen*/
  t8_dtri_last_descendant (&p->tri, &s->tri);
  t8_dline_last_descendant (&p->line, &s->line, level);
=======
  t8_dtri_last_descendant (&l->tri, &s->tri, level);
  t8_dline_last_descendant (&l->line, &s->line, level);
>>>>>>> johannes/feature-prism_elements_cxx
}

void
t8_dprism_vertex_coords (const t8_dprism_t * p, int vertex, int coords[3])
{
  T8_ASSERT (vertex >= 0 && vertex < 6);
  T8_ASSERT (p->line.level == p->tri.level);
  /*Compute x and y coordinate */
  t8_dtri_compute_coords (&p->tri, vertex % 3, coords);
  /*Compute z coordinate */
  t8_dline_vertex_coords (&p->line, vertex / 3, &coords[2]);
}

uint64_t
t8_dprism_linear_id (const t8_dprism_t * p, int level)
{
  uint64_t            id = 0;
  uint64_t            tri_id;
  uint64_t            line_id;
  int                 i;
  int                 prisms_of_size_i = 1;
  /*line_level = 2 ^ (level - 1) */
  int                 line_level = 1 << (level - 1);
  /*prism_shift = 2 * 8 ^ (level - 1) */
  int                 prism_shift = 4 * 1 << (3 * (level - 1));

  T8_ASSERT (0 <= level && level <= T8_DPRISM_MAXLEVEL);
  T8_ASSERT (p->line.level == p->tri.level);

  tri_id = t8_dtri_linear_id (&p->tri, level);
  line_id = t8_dline_linear_id (&p->line, level);
  for (i = 0; i < level; i++) {
    /*Compute via getting the local id of each ancestor triangle in which
     *elem->tri lies, the prism id, that elem would have, if it lies on the
     * lowest plane of the prism of level 0*/
    id += (tri_id % 4) * prisms_of_size_i;

    tri_id /= 4;
    prisms_of_size_i *= 8;
  }
  /*Now add the actual plane in which the prism is, which is computed via
   * line_id*/
  for (i = level - 1; i >= 0; i--) {
    /*The number to add to the id computed via the tri_id is 4*8^(level-i)
     *for each upper half in a prism of size i*/
    id += (line_id / line_level > 0) ? prism_shift : 0;
    line_id = (uint64_t) (line_id % line_level);
    prism_shift /= 8;
    line_level /= 2;
  }
  return id;
}
