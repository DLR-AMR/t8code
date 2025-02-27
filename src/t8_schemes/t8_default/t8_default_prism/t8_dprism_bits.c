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

#include <p4est_bits.h>
#include <sc_functions.h>
#include <t8_schemes/t8_default/t8_default_line/t8_dline_bits.h>
#include <t8_schemes/t8_default/t8_default_prism/t8_dprism_bits.h>
#include <t8_schemes/t8_default/t8_default_tri/t8_dtri_bits.h>

int t8_dprism_face_corners[5][4] = { { 1, 2, 4, 5 }, { 0, 2, 3, 5 }, { 0, 1, 3, 4 }, { 0, 1, 2, -1 }, { 3, 4, 5, -1 } };

int
t8_dprism_get_level (const t8_dprism_t *p)
{
  T8_ASSERT (p->line.level == p->tri.level);
  return p->line.level;
}

void
t8_dprism_copy (const t8_dprism_t *p, t8_dprism_t *dest)
{
  T8_ASSERT (p->line.level == p->tri.level);
  memcpy (dest, p, sizeof (t8_dprism_t));
  T8_ASSERT (dest->line.level == dest->tri.level);
}

int
t8_dprism_equal (const t8_dprism_t *elem1, const t8_dprism_t *elem2)
{
  return t8_dline_equal (&elem1->line, &elem2->line) && t8_dtri_equal (&elem1->tri, &elem2->tri);
}

int
t8_dprism_compare (const t8_dprism_t *p1, const t8_dprism_t *p2)
{
  int maxlvl;
  t8_linearidx_t id1, id2;
  T8_ASSERT (p1->line.level == p1->tri.level);
  T8_ASSERT (p2->line.level == p2->tri.level);

  maxlvl = SC_MAX (p1->line.level, p2->line.level);
  /* Compute the linear ids of the elements */
  id1 = t8_dprism_linear_id (p1, maxlvl);
  id2 = t8_dprism_linear_id (p2, maxlvl);
  if (id1 == id2) {
    /* The linear ids are the same, the prism with the smaller level
     * is considered smaller */
    return p1->line.level - p2->line.level;
  }
  /* return negative if id1 < id2, zero if id1 = id2, positive if id1 >
     id2 */
  return id1 < id2 ? -1 : id1 != id2;
}

void
t8_dprism_init_linear_id (t8_dprism_t *p, int level, t8_linearidx_t id)
{
  t8_linearidx_t tri_id = 0;
  t8_linearidx_t line_id = 0;
  int i;
  int triangles_of_size_i = 1;

  T8_ASSERT (0 <= level && level <= T8_DPRISM_MAXLEVEL);
  T8_ASSERT (id < sc_intpow64u (T8_DPRISM_CHILDREN, level));
  for (i = 0; i <= level; i++) {
    /*Get the number of the i-th prism and get the related triangle number
     * then multiplicate it by the number of triangles of level size.*/
    tri_id += (id % T8_DTRI_CHILDREN) * triangles_of_size_i;

    /*If id % 8 is larger than 3, the prism is in the upper part of the
     * parent prism. => line_id + 2^i*/
    line_id += (id % T8_DPRISM_CHILDREN) / T8_DTRI_CHILDREN * sc_intpow64u (T8_DLINE_CHILDREN, i);

    /*Each Prism divides into 8 children */
    id /= T8_DPRISM_CHILDREN;
    /*Each triangle divides into 4 children */
    triangles_of_size_i *= T8_DTRI_CHILDREN;
  }
  t8_dtri_init_linear_id (&p->tri, tri_id, level);
  t8_dline_init_linear_id (&p->line, level, line_id);

  T8_ASSERT (p->line.level == p->tri.level);
}

void
t8_dprism_parent (const t8_dprism_t *p, t8_dprism_t *parent)
{
  T8_ASSERT (p->line.level > 0);
  T8_ASSERT (p->line.level == p->tri.level);

  t8_dtri_parent (&p->tri, &parent->tri);
  t8_dline_parent (&p->line, &parent->line);

  T8_ASSERT (parent->line.level == parent->tri.level);
}

int
t8_dprism_child_id (const t8_dprism_t *p)
{
  int tri_child_id = t8_dtri_child_id (&p->tri);
  int line_child_id = t8_dline_child_id (&p->line);
  T8_ASSERT (p->line.level == p->tri.level);
  /*Prism in lower plane has the same id as the triangle, in the upper plane
   * it's a shift by the number of children a triangle has*/
  return tri_child_id + T8_DTRI_CHILDREN * line_child_id;
}

int
t8_dprism_is_familypv (t8_dprism_t **fam)
{
  t8_dtri_t *tri_fam[T8_DPRISM_CHILDREN];
  t8_dline_t *line_fam[T8_DLINE_CHILDREN];

  for (int i = 0; i < T8_DLINE_CHILDREN; i++) {
    for (int j = 0; j < T8_DTRI_CHILDREN; j++) {
      tri_fam[j] = &fam[j + i * T8_DTRI_CHILDREN]->tri;
    }
    if (!t8_dtri_is_familypv ((const t8_dtri_t **) tri_fam)) {
      return 0;
    }
  }

  for (int i = 0; i < T8_DTRI_CHILDREN; i++) {
    for (int j = 0; j < T8_DLINE_CHILDREN; j++) {
      line_fam[j] = &fam[j * T8_DTRI_CHILDREN + i]->line;
    }
    /* Proof for line_family and equality of triangles in both planes */
    if (!(t8_dline_is_familypv ((const t8_dline_t **) line_fam)
          && (fam[i]->tri.level == fam[i + T8_DTRI_CHILDREN]->tri.level)
          && (fam[i]->tri.type == fam[i + T8_DTRI_CHILDREN]->tri.type)
          && (fam[i]->tri.x == fam[i + T8_DTRI_CHILDREN]->tri.x)
          && (fam[i]->tri.y == fam[i + T8_DTRI_CHILDREN]->tri.y))) {
      return 0;
    }
  }

  for (int i = 0; i < T8_DPRISM_CHILDREN; i++) {
    if (fam[i]->line.level != fam[i]->tri.level) {
      return 0;
    }
  }
  return 1;
}

void
t8_dprism_nearest_common_ancestor (const t8_dprism_t *p1, const t8_dprism_t *p2, t8_dprism_t *r)
{
  int level;
  t8_dtri_nearest_common_ancestor (&p1->tri, &p2->tri, &r->tri);
  t8_dline_nearest_common_ancestor (&p1->line, &p2->line, &r->line);

  /* if the line and the triangle don't have the same level,
   * we compute the ancestor of the one with larger level at the
   * minimum level of the two. */
  if (r->tri.level != r->line.level) {
    level = SC_MIN (r->tri.level, r->line.level);
    if (r->tri.level > r->line.level) {
      /* triangle has larger level, compute its ancestor */
      t8_dtri_ancestor (&r->tri, level, &r->tri);
    }
    else {
      /* line has larger level, compute its ancestor */
      t8_dline_ancestor (&r->line, level, &r->line);
    }
  }
  T8_ASSERT (r->tri.level == r->line.level);
}

void
t8_dprism_boundary_face (const t8_dprism_t *p, int face, t8_element_t *boundary)
{
  T8_ASSERT (0 <= face && face < T8_DPRISM_FACES);
  p4est_quadrant_t *q = (p4est_quadrant_t *) boundary;
  if (face >= 3) {
    t8_dtri_t *t = (t8_dtri_t *) boundary;
    t8_dtri_copy (&p->tri, t);
    return;
  }
  switch (face) {
  case 0:
    q->x = ((int64_t) p->tri.y * P4EST_ROOT_LEN) / T8_DTRI_ROOT_LEN;
    q->y = ((int64_t) p->line.x * P4EST_ROOT_LEN) / T8_DLINE_ROOT_LEN;
    q->level = p->tri.level;
    break;
  case 1:
    q->x = ((int64_t) p->tri.x * P4EST_ROOT_LEN) / T8_DTRI_ROOT_LEN;
    q->y = ((int64_t) p->line.x * P4EST_ROOT_LEN) / T8_DLINE_ROOT_LEN;
    q->level = p->tri.level;
    break;
  case 2:
    q->x = ((int64_t) p->tri.x * P4EST_ROOT_LEN) / T8_DTRI_ROOT_LEN;
    q->y = ((int64_t) p->line.x * P4EST_ROOT_LEN) / T8_DLINE_ROOT_LEN;
    q->level = p->tri.level;
    break;
  default:
    SC_ABORT_NOT_REACHED ();
  }
}

int
t8_dprism_is_root_boundary (const t8_dprism_t *p, int face)
{
  T8_ASSERT (0 <= face && face < T8_DPRISM_FACES);
  /*face is not the bottom or top face of a prism */
  if (face < 3) {
    return t8_dtri_is_root_boundary (&p->tri, face);
  }
  else {
    /* Check whether the line is at the boundary, either at
     * its 0 face (face = 3) or
     * its 1 face (face = 4) */
    return t8_dline_is_root_boundary (&p->line, face != 3);
  }
}

int
t8_dprism_is_inside_root (t8_dprism_t *p)
{
  return t8_dtri_is_inside_root (&p->tri) && t8_dline_is_inside_root (&p->line);
}

void
t8_dprism_child (const t8_dprism_t *p, int childid, t8_dprism_t *child)
{
  T8_ASSERT (0 <= childid && childid < T8_DPRISM_CHILDREN);
  T8_ASSERT (p->line.level == p->tri.level);
  t8_dtri_child (&p->tri, childid % T8_DTRI_CHILDREN, &child->tri);
  t8_dline_child (&p->line, childid / T8_DTRI_CHILDREN, &child->line);
  T8_ASSERT (child->line.level == child->tri.level);
}

t8_element_shape_t
t8_dprism_face_shape (int face)
{
  T8_ASSERT (0 <= face && face < T8_DPRISM_FACES);
  if (face < 3) {
    return T8_ECLASS_QUAD;
  }
  return T8_ECLASS_TRIANGLE;
}

int
t8_dprism_num_face_children (int face)
{
  T8_ASSERT (0 <= face && face < T8_DPRISM_FACES);
  /*Bottom and top have T8_DTRI_CHILDREN, the other three faces depend on
     the children the triangle face has */
  return (face >= 3 ? T8_DTRI_CHILDREN : T8_DTRI_FACE_CHILDREN * T8_DLINE_CHILDREN);
}

int
t8_dprism_face_neighbour (const t8_dprism_t *p, int face, t8_dprism_t *neigh)
{
  T8_ASSERT (0 <= face && face < T8_DPRISM_FACES);
  if (face < 3) {
    t8_dline_copy (&p->line, &neigh->line);
    t8_dtri_face_neighbour (&p->tri, face, &neigh->tri);
    /* face neighbors face number:
     *  0 -> 2
     *  1 -> 1
     *  2 -> 0 */
    return 2 - face;
  }
  else if (face == 3) {
    t8_dtri_copy (&p->tri, &neigh->tri);
    t8_dline_face_neighbour (&p->line, &neigh->line, 0, NULL);
    return 4;
  }
  else {
    t8_dtri_copy (&p->tri, &neigh->tri);
    t8_dline_face_neighbour (&p->line, &neigh->line, 1, NULL);
    return 3;
  }
}

int
t8_dprism_get_face_corner (int face, int corner)
{
  T8_ASSERT (0 <= face && face < T8_DPRISM_FACES);
  T8_ASSERT (0 <= corner);
  T8_ASSERT ((corner <= 3 && face <= 2) || (corner <= 2));
  T8_ASSERT (t8_dprism_face_corners[face][corner] >= 0);

  return t8_dprism_face_corners[face][corner];
}

void
t8_dprism_childrenpv (const t8_dprism_t *p, int length, t8_dprism_t *c[])
{
  int i;
  T8_ASSERT (length == T8_DPRISM_CHILDREN);
  T8_ASSERT (p->line.level < T8_DPRISM_MAXLEVEL && p->tri.level == p->line.level);
  for (i = 7; i >= 0; i--) {
    t8_dprism_child (p, i, c[i]);
  }
}

int
t8_dprism_ancestor_id (t8_dprism_t *p, int level)
{
  return 4 * t8_dline_ancestor_id (&p->line, level) + t8_dtri_ancestor_id (&p->tri, level);
}

const int children_at_face[2][12] = { { 1, 3, 5, 7, 0, 3, 4, 7, 0, 1, 4, 5 }, { 2, 3, 6, 7, 0, 3, 4, 7, 0, 2, 4, 6 } };

void
t8_dprism_children_at_face (const t8_dprism_t *p, int face, t8_dprism_t **children, int num_children,
                            int *child_indices)
{
  T8_ASSERT (num_children == t8_dprism_num_face_children (face));
  T8_ASSERT (0 <= face && face < T8_DPRISM_FACES);
  if (face < 3) {
    for (int ichild = 3; ichild >= 0; ichild--) {
      t8_dprism_child (p, children_at_face[p->tri.type][face * 4 + ichild], children[ichild]);
    }
  }
  else {
    for (int ichild = 3; ichild >= 0; ichild--) {
      t8_dprism_child (p, (face % 3) * 4 + ichild, children[ichild]);
    }
  }
  /* Fill child-indices array */
  if (child_indices != NULL) {
    if (face < 3) {
      for (int ichild = 3; ichild >= 0; ichild--) {
        child_indices[ichild] = children_at_face[p->tri.type][face * 4 + ichild];
      }
    }
    else {
      for (int ichild = 3; ichild >= 0; ichild--) {
        child_indices[ichild] = (face % 3) * 4 + ichild;
      }
    }
  }
}

int
t8_dprism_face_child_face (int face)
{
  T8_ASSERT (0 <= face && face < T8_DPRISM_FACES);
  /* For prisms the face number of the children is the same as the one
   * of the parent */
  return face;
}

int
t8_dprism_face_parent_face (const t8_dprism_t *prism, int face)
{
  T8_ASSERT (0 <= face && face < T8_DPRISM_FACES);
  if (face < 3 && t8_dtri_face_parent_face (&prism->tri, face) != -1) {
    return face;
  }
  /*prism_face 3 = line_face 0, prism_face 4 = line_face 1 */
  if (face >= 3 && t8_dline_face_parent_face (&prism->line, face - 3) != -1) {
    return face;
  }
  else {
    return -1;
  }
}

int
t8_dprism_tree_face (int face)
{
  T8_ASSERT (0 <= face && face < T8_DPRISM_FACES);
  /*For prisms, the face number coincides with the number of the root
     tree face */
  return face;
}

void
t8_dprism_extrude_face (const t8_element_t *face, t8_element_t *elem, const int root_face)
{
  t8_dprism_t *p = (t8_dprism_t *) elem;
  const t8_dtri_t *t = (const t8_dtri_t *) face;
  const p4est_quadrant_t *q = (const p4est_quadrant_t *) face;
  T8_ASSERT (0 <= root_face && root_face < T8_DPRISM_FACES);
  switch (root_face) {
  case 0:
    p->tri.type = 0;
    p->line.level = q->level;
    p->tri.level = q->level;
    p->tri.x = T8_DTRI_ROOT_LEN - T8_DTRI_LEN (p->tri.level);
    p->tri.y = ((int64_t) q->x * T8_DTRI_ROOT_LEN) / P4EST_ROOT_LEN;
    p->line.x = ((int64_t) q->y * T8_DLINE_ROOT_LEN) / P4EST_ROOT_LEN;
    break;
  case 1:
    p->tri.type = 0;
    p->line.level = q->level;
    p->tri.level = q->level;
    p->tri.x = ((int64_t) q->x * T8_DTRI_ROOT_LEN) / P4EST_ROOT_LEN;
    p->tri.y = ((int64_t) q->x * T8_DTRI_ROOT_LEN) / P4EST_ROOT_LEN;
    p->line.x = ((int64_t) q->y * T8_DLINE_ROOT_LEN) / P4EST_ROOT_LEN;
    break;
  case 2:
    p->tri.type = 0;
    p->line.level = q->level;
    p->tri.level = q->level;
    p->tri.x = ((int64_t) q->x * T8_DTRI_ROOT_LEN) / P4EST_ROOT_LEN;
    p->tri.y = 0;
    p->line.x = ((int64_t) q->y * T8_DLINE_ROOT_LEN) / P4EST_ROOT_LEN;
    break;
  case 3:
    p->tri.type = t->type;
    p->line.level = t->level;
    p->tri.level = t->level;
    p->tri.x = t->x;
    p->tri.y = t->y;
    p->line.x = 0;
    break;
  case 4:
    p->tri.type = t->type;
    p->line.level = t->level;
    p->tri.level = t->level;
    p->tri.x = t->x;
    p->tri.y = t->y;
    p->line.x = T8_DLINE_ROOT_LEN - T8_DLINE_LEN (p->line.level);
    break;
  default:
    SC_ABORT_NOT_REACHED ();
  }
}

void
t8_dprism_successor (const t8_dprism_t *p, t8_dprism_t *succ, int level)
{
  int prism_child_id;
  t8_dprism_copy (p, succ);
  /*update the level */
  succ->line.level = level;
  succ->tri.level = level;
  prism_child_id = t8_dprism_child_id (succ);
  T8_ASSERT (1 <= level && level <= T8_DPRISM_MAXLEVEL);
  T8_ASSERT (p->line.level == p->tri.level);
  /*The next prism is the child with local ID 0 of the next parent prism */
  if (prism_child_id == T8_DPRISM_CHILDREN - 1) {
    t8_dprism_successor (p, succ, level - 1);
    /*Zero out the bits of higher level, caused by recursion */
    succ->tri.x = (succ->tri.x >> (T8_DTRI_MAXLEVEL - level + 1)) << (T8_DTRI_MAXLEVEL - level + 1);
    succ->tri.y = (succ->tri.y >> (T8_DTRI_MAXLEVEL - level + 1)) << (T8_DTRI_MAXLEVEL - level + 1);
    succ->line.x = (succ->line.x >> (T8_DLINE_MAXLEVEL - level + 1)) << (T8_DLINE_MAXLEVEL - level + 1);
    /*Set the level to the actual level */
    succ->line.level = level;
    succ->tri.level = level;
  }
  /*The next prism is one plane up, local_tri_id = 0 */
  else if ((prism_child_id + 1) % T8_DTRI_CHILDREN == 0) {
    /*parent is computed with succ, cause there are the updated data */
    t8_dprism_parent (succ, succ);
    t8_dprism_child (succ, prism_child_id + 1, succ);
  }
  /*The next Prism is in the same plane, but has the next base-triangle */
  else {
    t8_dtri_successor (&p->tri, &succ->tri, level);
  }
  T8_ASSERT (succ->line.level == succ->tri.level);
}

void
t8_dprism_first_descendant (const t8_dprism_t *p, t8_dprism_t *s, int level)
{
  T8_ASSERT (level >= p->line.level && level <= T8_DPRISM_MAXLEVEL);
  T8_ASSERT (p->line.level == p->tri.level);
  /*First prism descendant = first triangle desc x first line desc */
  t8_dtri_first_descendant (&p->tri, &s->tri, level);
  t8_dline_first_descendant (&p->line, &s->line, level);
#ifdef T8_ENABLE_DEBUG
  {
    t8_linearidx_t id;
    id = t8_dprism_linear_id (p, level);
    T8_ASSERT (t8_dprism_linear_id (s, level) == id);
  }
#endif

  T8_ASSERT (s->line.level == s->tri.level);
  T8_ASSERT (s->line.level == level);
}

void
t8_dprism_last_descendant (const t8_dprism_t *p, t8_dprism_t *s, int level)
{
  T8_ASSERT (level >= p->line.level && level <= T8_DPRISM_MAXLEVEL);
  T8_ASSERT (p->line.level == p->tri.level);
  /*Last prism descendant = last triangle desc x last line desc */
  T8_ASSERT (level <= T8_DTRI_MAXLEVEL);
  t8_dtri_last_descendant (&p->tri, &s->tri, level);
  t8_dline_last_descendant (&p->line, &s->line, level);
  T8_ASSERT (s->line.level == s->tri.level);
}

void
t8_dprism_corner_descendant (const t8_dprism_t *p, t8_dprism_t *s, int corner, int level)
{
  T8_ASSERT (p->tri.level <= level && level <= T8_DPRISM_MAXLEVEL);
  T8_ASSERT (0 <= corner && corner < T8_DPRISM_CORNERS);
  t8_dtri_corner_descendant (&p->tri, &s->tri, corner % 3, level);
  if (corner < 3) {
    t8_dline_first_descendant (&p->line, &s->line, level);
  }
  else {
    t8_dline_last_descendant (&p->line, &s->line, level);
  }
}

void
t8_dprism_vertex_integer_coords (const t8_dprism_t *elem, const int vertex, int coords[3])
{
  T8_ASSERT (vertex >= 0 && vertex < 6);
  T8_ASSERT (elem->line.level == elem->tri.level);
  /*Compute x and y coordinate */
  t8_dtri_compute_integer_coords (&elem->tri, vertex % 3, coords);
  /* Compute z coordinate coords[0] *= T8_DPRISM_ROOT_BY_DTRI_ROOT; */
  t8_dline_vertex_integer_coords (&elem->line, vertex / 3, &coords[2]);
  coords[0] /= T8_DPRISM_ROOT_BY_DTRI_ROOT;
  coords[1] /= T8_DPRISM_ROOT_BY_DTRI_ROOT;
  coords[2] /= T8_DPRISM_ROOT_BY_DLINE_ROOT;
}

void
t8_dprism_vertex_ref_coords (const t8_dprism_t *elem, const int vertex, double coords[3])
{
  int coords_int[3];
  T8_ASSERT (t8_dprism_is_valid (elem));
  T8_ASSERT (vertex >= 0 && vertex < 6);

  /* Compute the integer coordinates in [0, root_len]^3 */
  t8_dprism_vertex_integer_coords (elem, vertex, coords_int);

  /* Divide by the root length. */
  coords[0] = coords_int[0] / (double) T8_DPRISM_ROOT_LEN;
  coords[1] = coords_int[1] / (double) T8_DPRISM_ROOT_LEN;
  coords[2] = coords_int[2] / (double) T8_DPRISM_ROOT_LEN;
}

void
t8_dprism_compute_reference_coords (const t8_dprism_t *elem, const double *ref_coords, const size_t num_coords,
                                    double *out_coords)
{
  T8_ASSERT (t8_dprism_is_valid (elem));
  T8_ASSERT (elem->line.level == elem->tri.level);
  /*Compute x and y coordinate */
  t8_dtri_compute_reference_coords (&elem->tri, ref_coords, num_coords, 1, out_coords);
  /*Compute z coordinate */
  t8_dline_compute_reference_coords (&elem->line, ref_coords + 2, num_coords, 2, out_coords + 2);
}

t8_linearidx_t
t8_dprism_linear_id (const t8_dprism_t *p, int level)
{
  t8_linearidx_t id = 0;
  t8_linearidx_t tri_id;
  t8_linearidx_t line_id;
  int i;
  t8_linearidx_t prisms_of_size_i = 1;
  /*lines_at_level = Num_of_Line_children ^ (level - 1) */
  t8_linearidx_t lines_at_level;
  /*prism_shift = Num_of_Prism_children / Num_of_line-Children * 8 ^ (level - 1) */
  t8_linearidx_t prism_shift;
  T8_ASSERT (0 <= level && level <= T8_DPRISM_MAXLEVEL);
  T8_ASSERT (p->line.level == p->tri.level);
  /*id = 0 for root element */
  if (level == 0) {
    return 0;
  }

  lines_at_level = sc_intpow64u (T8_DLINE_CHILDREN, level - 1);
  prism_shift = (T8_DPRISM_CHILDREN / T8_DLINE_CHILDREN) * sc_intpow64u (T8_DPRISM_CHILDREN, level - 1);
  tri_id = t8_dtri_linear_id (&p->tri, level);
  line_id = t8_dline_linear_id (&p->line, level);
  for (i = 0; i < level; i++) {
    /*Compute via getting the local id of each ancestor triangle in which
     *elem->tri lies, the prism id, that elem would have, if it lies on the
     * lowest plane of the prism of level 0 */
    id += (tri_id % T8_DTRI_CHILDREN) * prisms_of_size_i;
    tri_id /= T8_DTRI_CHILDREN;
    prisms_of_size_i *= T8_DPRISM_CHILDREN;
  }
  /*Now add the actual plane in which the prism is, which is computed via
   * line_id*/
  for (i = level - 1; i >= 0; i--) {
    /*The number to add to the id computed via the tri_id is 4*8^(level-i)
     *for each plane in a prism of size i*/
    id += line_id / lines_at_level * prism_shift;
    line_id = (line_id % lines_at_level);
    prism_shift /= T8_DPRISM_CHILDREN;
    lines_at_level /= T8_DLINE_CHILDREN;
  }

  return id;
}

/* Returns true if and only if p is a valid prism,
 * that is its triangle and line part are valid, and they have the same
 * refinement level. */
int
t8_dprism_is_valid (const t8_dprism_t *p)
{
  const t8_dline_t *line = &p->line;
  const t8_dtri_t *tri = &p->tri;
  const int same_level = line->level == tri->level;
  return t8_dtri_is_valid (tri) && t8_dline_is_valid (line) && same_level;
}
