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
#include "t8_dline_bits.h"
#include "t8_default_common_cxx.hxx"
#include "t8_default_quad_cxx.hxx"

/* We want to export the whole implementation to be callable from "C" */
T8_EXTERN_C_BEGIN ();

/* This function is used by other element functions and we thus need to
 * declare it up here */
uint64_t            t8_element_get_linear_id (const t8_element_t * elem,
                                              int level);

#ifdef T8_ENABLE_DEBUG

int
t8_element_surround_matches (const p4est_quadrant_t * q,
                             const p4est_quadrant_t * r)
{
  return T8_QUAD_GET_TDIM (q) == T8_QUAD_GET_TDIM (r) &&
    (T8_QUAD_GET_TDIM (q) == -1 ||
     (T8_QUAD_GET_TNORMAL (q) == T8_QUAD_GET_TNORMAL (r) &&
      T8_QUAD_GET_TCOORD (q) == T8_QUAD_GET_TCOORD (r)));
}

#endif /* T8_ENABLE_DEBUG */

int
t8_default_scheme_quad_c::t8_element_maxlevel (void)
{
  return P4EST_QMAXLEVEL;
}

/* *INDENT-OFF* */
t8_eclass_t
t8_default_scheme_quad_c::t8_element_child_eclass (int childid)
/* *INDENT-ON* */

{
  T8_ASSERT (0 <= childid && childid < P4EST_CHILDREN);

  return T8_ECLASS_QUAD;
}

int
t8_default_scheme_quad_c::t8_element_level (const t8_element_t * elem)
{
  return (int) ((const p4est_quadrant_t *) elem)->level;
}

static void
t8_element_copy_surround (const p4est_quadrant_t * q, p4est_quadrant_t * r)
{
  T8_QUAD_SET_TDIM (r, T8_QUAD_GET_TDIM (q));
  if (T8_QUAD_GET_TDIM (q) == 3) {
    T8_QUAD_SET_TNORMAL (r, T8_QUAD_GET_TNORMAL (q));
    T8_QUAD_SET_TCOORD (r, T8_QUAD_GET_TCOORD (q));
  }
}

void
t8_default_scheme_quad_c::t8_element_copy (const t8_element_t * source,
                                           t8_element_t * dest)
{
  const p4est_quadrant_t *q = (const p4est_quadrant_t *) source;
  p4est_quadrant_t   *r = (p4est_quadrant_t *) dest;

  if (r == q) {
    /* Do nothing if they are already the same quadrant. */
    return;
  }
  *r = *q;
  t8_element_copy_surround (q, r);
}

int
t8_default_scheme_quad_c::t8_element_compare (const t8_element_t * elem1,
                                              const t8_element_t * elem2)
{
  int                 maxlvl;
  u_int64_t           id1, id2;

  /* Compute the bigger level of the two */
  maxlvl = SC_MAX (t8_element_level (elem1), t8_element_level (elem2));
  /* Compute the linear ids of the elements */
  id1 = t8_element_get_linear_id (elem1, maxlvl);
  id2 = t8_element_get_linear_id (elem2, maxlvl);
  /* return negativ if id1 < id2, zero if id1 = id2, positive if id1 > id2 */
  return id1 < id2 ? -1 : id1 != id2;
}

void
t8_default_scheme_quad_c::t8_element_parent (const t8_element_t * elem,
                                             t8_element_t * parent)
{
  const p4est_quadrant_t *q = (const p4est_quadrant_t *) elem;
  p4est_quadrant_t   *r = (p4est_quadrant_t *) parent;

  p4est_quadrant_parent (q, r);
  t8_element_copy_surround (q, r);
}

void
t8_default_scheme_quad_c::t8_element_sibling (const t8_element_t * elem,
                                              int sibid,
                                              t8_element_t * sibling)
{
  const p4est_quadrant_t *q = (const p4est_quadrant_t *) elem;
  p4est_quadrant_t   *r = (p4est_quadrant_t *) sibling;

  p4est_quadrant_sibling (q, r, sibid);
  t8_element_copy_surround (q, r);
}

int
t8_default_scheme_quad_c::t8_element_num_faces (const t8_element_t * elem)
{
  return P4EST_FACES;
}

int
t8_default_scheme_quad_c::t8_element_num_children (const t8_element_t * elem)
{
  return P4EST_CHILDREN;
}

int
t8_default_scheme_quad_c::t8_element_num_face_children (const t8_element_t *
                                                        elem, int face)
{
  return 2;
}

void
t8_default_scheme_quad_c::t8_element_child (const t8_element_t * elem,
                                            int childid, t8_element_t * child)
{
  const p4est_quadrant_t *q = (const p4est_quadrant_t *) elem;
  const p4est_qcoord_t shift = P4EST_QUADRANT_LEN (q->level + 1);
  p4est_quadrant_t   *r = (p4est_quadrant_t *) child;

  T8_ASSERT (p4est_quadrant_is_extended (q));
  T8_ASSERT (q->level < P4EST_QMAXLEVEL);
  T8_ASSERT (childid >= 0 && childid < P4EST_CHILDREN);

  r->x = childid & 0x01 ? (q->x | shift) : q->x;
  r->y = childid & 0x02 ? (q->y | shift) : q->y;
  r->level = q->level + 1;
  T8_ASSERT (p4est_quadrant_is_parent (q, r));

  t8_element_copy_surround (q, r);
}

void
t8_default_scheme_quad_c::t8_element_children (const t8_element_t * elem,
                                               int length, t8_element_t * c[])
{
  const p4est_quadrant_t *q = (const p4est_quadrant_t *) elem;
  int                 i;

  T8_ASSERT (length == P4EST_CHILDREN);

  p4est_quadrant_childrenpv (q, (p4est_quadrant_t **) c);
  for (i = 0; i < P4EST_CHILDREN; ++i) {
    t8_element_copy_surround (q, (p4est_quadrant_t *) c[i]);
  }
}

int
t8_default_scheme_quad_c::t8_element_child_id (const t8_element_t * elem)
{
  return p4est_quadrant_child_id ((p4est_quadrant_t *) elem);
}

int
t8_default_scheme_quad_c::t8_element_is_family (t8_element_t ** fam)
{
  return p4est_quadrant_is_familypv ((p4est_quadrant_t **) fam);
}

void
t8_default_scheme_quad_c::t8_element_set_linear_id (t8_element_t * elem,
                                                    int level, uint64_t id)
{
  T8_ASSERT (0 <= level && level <= P4EST_QMAXLEVEL);
  T8_ASSERT (0 <= id && id < ((uint64_t) 1) << P4EST_DIM * level);

  p4est_quadrant_set_morton ((p4est_quadrant_t *) elem, level, id);
  T8_QUAD_SET_TDIM ((p4est_quadrant_t *) elem, 2);
}

uint64_t
  t8_default_scheme_quad_c::t8_element_get_linear_id (const t8_element_t *
                                                      elem, int level)
{
  T8_ASSERT (0 <= level && level <= P4EST_QMAXLEVEL);

  return p4est_quadrant_linear_id ((p4est_quadrant_t *) elem, level);
}

void
t8_default_scheme_quad_c::t8_element_first_descendant (const t8_element_t *
                                                       elem,
                                                       t8_element_t * desc)
{
  p4est_quadrant_first_descendant ((p4est_quadrant_t *) elem,
                                   (p4est_quadrant_t *) desc,
                                   P4EST_QMAXLEVEL);
}

void
t8_default_scheme_quad_c::t8_element_last_descendant (const t8_element_t *
                                                      elem,
                                                      t8_element_t * desc)
{
  p4est_quadrant_last_descendant ((p4est_quadrant_t *) elem,
                                  (p4est_quadrant_t *) desc, P4EST_QMAXLEVEL);
}

void
t8_default_scheme_quad_c::t8_element_successor (const t8_element_t * elem1,
                                                t8_element_t * elem2,
                                                int level)
{
  uint64_t            id;
  T8_ASSERT (0 <= level && level <= P4EST_QMAXLEVEL);

  id = p4est_quadrant_linear_id ((const p4est_quadrant_t *) elem1, level);
  T8_ASSERT (id + 1 < ((uint64_t) 1) << P4EST_DIM * level);
  p4est_quadrant_set_morton ((p4est_quadrant_t *) elem2, level, id + 1);
  t8_element_copy_surround ((const p4est_quadrant_t *) elem1,
                            (p4est_quadrant_t *) elem2);
}

void
t8_default_scheme_quad_c::t8_element_nca (const t8_element_t * elem1,
                                          const t8_element_t * elem2,
                                          t8_element_t * nca)
{
  const p4est_quadrant_t *q1 = (const p4est_quadrant_t *) elem1;
  const p4est_quadrant_t *q2 = (const p4est_quadrant_t *) elem2;
  p4est_quadrant_t   *r = (p4est_quadrant_t *) nca;

  T8_ASSERT (t8_element_surround_matches (q1, q2));

  p4est_nearest_common_ancestor (q1, q2, r);
  t8_element_copy_surround (q1, r);
}

t8_eclass_t
  t8_default_scheme_quad_c::t8_element_face_class (const t8_element_t * elem,
                                                   int face)
{
  return T8_ECLASS_LINE;
}

void
t8_default_scheme_quad_c::t8_element_children_at_face (const t8_element_t *
                                                       elem, int face,
                                                       t8_element_t *
                                                       children[],
                                                       int num_children)
{
  int                 first_child, second_child;

  T8_ASSERT (0 <= face && face < P4EST_FACES);
  T8_ASSERT (num_children == t8_element_num_face_children (elem, face));

  /*
   * Compute the child id of the first and second child at the face.
   *
   *            3
   *
   *      x - - x - - x           This picture shows a refined quadrant
   *      |     |     |           with child_ids and the label for the faces.
   *      | 2   | 3   |           For examle for face 2 (bottom face) we see
   * 0    x - - x - - x   1       first_child = 0 and second_child = 1.
   *      |     |     |
   *      | 0   | 1   |
   *      x - - x - - x
   *
   *            2
   */
  /* TODO: Think about a short and easy bitwise formula. */
  switch (face) {
  case 0:
    first_child = 0;
    second_child = 2;
    break;
  case 1:
    first_child = 1;
    second_child = 3;
    break;
  case 2:
    first_child = 0;
    second_child = 1;
    break;
  case 3:
    first_child = 2;
    second_child = 3;
    break;
  default:
    SC_ABORT_NOT_REACHED ();
  }

  /* From the child ids we now construct the children at the faces. */
  /* We have to revert the order and compute second child first, since
   * the usage allows for elem == children[0].
   */
  t8_element_child (elem, second_child, children[1]);
  t8_element_child (elem, first_child, children[0]);
}

int
t8_default_scheme_quad_c::t8_element_face_child_face (const t8_element_t *
                                                      elem, int face,
                                                      int face_child)
{
  /* For quadrants the face enumeration of children is the same as for the parent. */
  return face;
}

void
t8_default_scheme_quad_c::t8_element_transform_face (const t8_element_t *
                                                     elem1,
                                                     t8_element_t * elem2,
                                                     int orientation,
                                                     int is_smaller_face)
{
  const p4est_quadrant_t *q = (const p4est_quadrant_t *) elem1;
  p4est_quadrant_t   *p = (p4est_quadrant_t *) elem2;
  p4est_qcoord_t      h = P4EST_QUADRANT_LEN (q->level);
  T8_ASSERT (0 <= orientation && orientation < P4EST_FACES);

  p->level = q->level;
  /*
   * The faces of the root quadrant are enumerated like this:
   *
   *   v_2      v_3
   *     x -->-- x
   *     |       |
   *     ^       ^
   *     |       |
   *     x -->-- x
   *   v_0      v_1
   *
   * Orientation is the corner number of the bigger face that coincides
   * with the corner v_0 of the smaller face.
   */
  switch (orientation) {
  case 0:                      /* Nothing to do */
    break;
  case 1:
    p->x = P4EST_ROOT_LEN - q->x - h;
    /* p->y remains q->y */
    break;
  case 2:
    /* p->x remains q->x */
    p->y = P4EST_ROOT_LEN - q->y - h;
    break;
  case 3:
    p->x = P4EST_ROOT_LEN - q->y - h;
    p->y = P4EST_ROOT_LEN - q->x - h;
    break;
  default:
    SC_ABORT_NOT_REACHED ();
  }
}

void
t8_default_scheme_quad_c::t8_element_extrude_face (const t8_element_t * face,
                                                   const t8_eclass_scheme_c *
                                                   face_scheme,
                                                   t8_element_t * elem,
                                                   int root_face)
{
  const t8_dline_t   *l = (const t8_dline_t *) face;
  p4est_quadrant_t   *q = (p4est_quadrant_t *) elem;

  T8_ASSERT (face_scheme->eclass == T8_ECLASS_LINE);
  T8_ASSERT (face_scheme->t8_element_is_valid (elem));
  T8_ASSERT (0 <= root_face && root_face < P4EST_FACES);
  /*
   * The faces of the root quadrant are enumerated like this:
   *
   *        f_2
   *     x -->-- x
   *     |       |
   *     ^       ^
   * f_0 |       | f_1
   *     x -->-- x
   *        f_3
   *
   * The arrows >,^ denote the orientation of the faces.
   * We need to scale the coordinates since a root line may have a different
   * length than a root quad.
   */
  q->level = l->level;
  switch (root_face) {
  case 0:
    q->x = 0;
    q->y = ((int64_t) l->x * P4EST_ROOT_LEN) / T8_DLINE_ROOT_LEN;
    break;
  case 1:
    q->x = P4EST_LAST_OFFSET (q->level);
    q->y = ((int64_t) l->x * P4EST_ROOT_LEN) / T8_DLINE_ROOT_LEN;
    break;
  case 2:
    q->x = ((int64_t) l->x * P4EST_ROOT_LEN) / T8_DLINE_ROOT_LEN;
    q->y = 0;
    break;
  case 3:
    q->x = ((int64_t) l->x * P4EST_ROOT_LEN) / T8_DLINE_ROOT_LEN;
    q->y = P4EST_LAST_OFFSET (q->level);
    break;
  default:
    SC_ABORT_NOT_REACHED ();
  }
}

int
t8_default_scheme_quad_c::t8_element_tree_face (const t8_element_t * elem,
                                                int face)
{
  T8_ASSERT (0 <= face && face < P4EST_FACES);
  /* For quadrants the face and the tree face number are the same. */
  return face;
}

void
t8_default_scheme_quad_c::t8_element_boundary_face (const t8_element_t * elem,
                                                    int face,
                                                    t8_element_t * boundary,
                                                    const t8_eclass_scheme_c *
                                                    boundary_scheme)
{
  const p4est_quadrant_t *q = (const p4est_quadrant_t *) elem;
  t8_dline_t         *l = (t8_dline_t *) boundary;

  T8_ASSERT (boundary_scheme->eclass == T8_ECLASS_LINE);
  T8_ASSERT (boundary_scheme->t8_element_is_valid (boundary));
  T8_ASSERT (0 <= face && face < P4EST_FACES);
  /* The level of the boundary element is the same as the quadrant's level */
  l->level = q->level;
  /*
   * The faces of the quadrant are enumerated like this:
   *        f_2
   *     x ---- x
   *     |      |
   * f_0 |      | f_1
   *     x ---- x
   *        f_3
   *
   * If face = 0 or face = 1 then l->x = q->y
   * if face = 2 or face = 3 then l->x = q->x
   */
  l->x = ((face >> 1 ? q->x : q->y) *
          ((int64_t) T8_DLINE_ROOT_LEN) / P4EST_ROOT_LEN);
}

void
t8_default_scheme_quad_c::t8_element_boundary (const t8_element_t * elem,
                                               int min_dim, int length,
                                               t8_element_t ** boundary)
{
  SC_ABORT ("Not implemented\n");
#if 0
#ifdef T8_ENABLE_DEBUG
  int                 per_eclass[T8_ECLASS_COUNT];
#endif
  int                 iface;

  T8_ASSERT (length ==
             t8_eclass_count_boundary (T8_ECLASS_QUAD, min_dim, per_eclass));

  /* TODO: write this function */

  T8_ASSERT (length == P4EST_FACES);
  for (iface = 0; iface < P4EST_FACES; iface++) {
    t8_element_boundary_face (elem, iface, boundary[iface]);
  }
#endif
}

int
t8_default_scheme_quad_c::t8_element_is_root_boundary (const t8_element_t *
                                                       elem, int face)
{
  const p4est_quadrant_t *q = (const p4est_quadrant_t *) elem;
  p4est_qcoord_t      coord;

  T8_ASSERT (0 <= face && face < P4EST_FACES);

  /* if face is 0 or 1 q->x
   *            2 or 3 q->y
   */
  coord = face >> 1 ? q->y : q->x;
  /* If face is 0 or 2 check against 0.
   * If face is 1 or 3  check against LAST_OFFSET */
  return coord == (face & 1 ? P4EST_LAST_OFFSET (q->level) : 0);
}

int
t8_default_scheme_quad_c::t8_element_face_neighbor_inside (const t8_element_t
                                                           * elem,
                                                           t8_element_t *
                                                           neigh, int face)
{
  const p4est_quadrant_t *q = (const p4est_quadrant_t *) elem;
  p4est_quadrant_t   *n = (p4est_quadrant_t *) neigh;

  T8_ASSERT (0 <= face && face < P4EST_FACES);
  p4est_quadrant_face_neighbor (q, face, n);
  /* return true if neigh is inside the root */
  return p4est_quadrant_is_inside_root (n);
}

void
t8_default_scheme_quad_c::t8_element_anchor (const t8_element_t * elem,
                                             int coord[3])
{
  p4est_quadrant_t   *q;

  q = (p4est_quadrant_t *) elem;
  coord[0] = q->x;
  coord[1] = q->y;
  coord[2] = 0;
}

int
t8_default_scheme_quad_c::t8_element_root_len (const t8_element_t * elem)
{
  return P4EST_ROOT_LEN;
}

void
t8_default_scheme_quad_c::t8_element_vertex_coords (const t8_element_t * t,
                                                    int vertex, int coords[])
{
  const p4est_quadrant_t *q1 = (const p4est_quadrant_t *) t;
  int                 len;

  T8_ASSERT (0 <= vertex && vertex < 4);
  /* Get the length of the quadrant */
  len = P4EST_QUADRANT_LEN (q1->level);
  /* Compute the x and y coordinates of the vertex depending on the
   * vertex number */
  coords[0] = q1->x + (vertex & 1 ? 1 : 0) * len;
  coords[1] = q1->y + (vertex & 2 ? 1 : 0) * len;
}

void
t8_default_scheme_quad_c::t8_element_new (int length, t8_element_t ** elem)
{
  /* allocate memory for a tet */
  t8_default_scheme_common_c::t8_element_new (length, elem);

  /* in debug mode, set sensible default values. */
#ifdef T8_ENABLE_DEBUG
  {
    int                 i;
    for (i = 0; i < length; i++) {
      t8_element_init (1, elem[i], 0);
    }
  }
#endif
}

void
t8_default_scheme_quad_c::t8_element_init (int length, t8_element_t * elem,
                                           int new_called)
{
#ifdef T8_ENABLE_DEBUG
  if (!new_called) {
    int                 i;
    p4est_quadrant_t   *quads = (p4est_quadrant_t *) elem;
    for (i = 0; i < length; i++) {
      P4EST_QUADRANT_INIT (quads + i);
    }
  }
#endif
}

#ifdef T8_ENABLE_DEBUG
/* *INDENT-OFF* */
/* indent bug, indent adds a second "const" modifier */
int
t8_default_scheme_quad_c::t8_element_is_valid (const t8_element_t * elem) const
/* *INDENT-ON* */
{
  /* TODO: additional checks? do we set pad8 or similar?
   */
  return p4est_quadrant_is_valid ((const p4est_quadrant_t *) elem);
}
#endif

/* Constructor */
t8_default_scheme_quad_c::t8_default_scheme_quad_c (void)
{
  eclass = T8_ECLASS_QUAD;
  element_size = sizeof (t8_pquad_t);
  ts_context = sc_mempool_new (sizeof (element_size));
}

t8_default_scheme_quad_c::~t8_default_scheme_quad_c ()
{
  /* This destructor is empty since the destructor of the
   * default_common scheme is called automatically and it
   * suffices to destroy the quad_scheme.
   * However we need to provide an implementation of the destructor
   * and hence this empty function. */
}

T8_EXTERN_C_END ();
