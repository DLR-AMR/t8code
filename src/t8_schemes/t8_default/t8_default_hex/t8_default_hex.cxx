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

#include <p8est_bits.h>
#include <p4est_bits.h>
#include <t8_schemes/t8_default/t8_default_common/t8_default_common.hxx>
#include <t8_schemes/t8_default/t8_default_hex/t8_default_hex.hxx>

#define HEX_LINEAR_MAXLEVEL P8EST_OLD_QMAXLEVEL
#define HEX_REFINE_MAXLEVEL P8EST_OLD_QMAXLEVEL
/* This is the ideal but currently not exposed in the overall interface. */
/* #define HEX_REFINE_MAXLEVEL P8EST_QMAXLEVEL */

/* We want to export the whole implementation to be callable from "C" */
T8_EXTERN_C_BEGIN ();

int
t8_default_scheme_hex_c::get_maxlevel (void) const
{
  return HEX_REFINE_MAXLEVEL;
}

int
t8_default_scheme_hex_c::element_get_level (const t8_element_t *elem) const
{
  T8_ASSERT (element_is_valid (elem));
  return (int) ((const p8est_quadrant_t *) elem)->level;
}

void
t8_default_scheme_hex_c::element_copy (const t8_element_t *source, t8_element_t *dest) const
{
  T8_ASSERT (element_is_valid (source));
  T8_ASSERT (element_is_valid (dest));
  *(p8est_quadrant_t *) dest = *(const p8est_quadrant_t *) source;
}

int
t8_default_scheme_hex_c::element_compare (const t8_element_t *elem1, const t8_element_t *elem2) const
{
  T8_ASSERT (element_is_valid (elem1));
  T8_ASSERT (element_is_valid (elem2));

  return p8est_quadrant_compare ((const p8est_quadrant_t *) elem1, (const p8est_quadrant_t *) elem2);
}

int
t8_default_scheme_hex_c::element_is_equal (const t8_element_t *elem1, const t8_element_t *elem2) const
{
  return p8est_quadrant_is_equal ((const p8est_quadrant_t *) elem1, (const p8est_quadrant_t *) elem2);
}

void
t8_default_scheme_hex_c::element_get_parent (const t8_element_t *elem, t8_element_t *parent) const
{
  T8_ASSERT (element_is_valid (elem));
  T8_ASSERT (element_is_valid (parent));
  p8est_quadrant_parent ((const p8est_quadrant_t *) elem, (p8est_quadrant_t *) parent);
}

void
t8_default_scheme_hex_c::element_get_sibling (const t8_element_t *elem, int sibid, t8_element_t *sibling) const
{
  T8_ASSERT (element_is_valid (elem));
  T8_ASSERT (element_is_valid (sibling));
  p8est_quadrant_sibling ((const p8est_quadrant_t *) elem, (p8est_quadrant_t *) sibling, sibid);
}

int
t8_default_scheme_hex_c::element_get_num_faces (const t8_element_t *elem) const
{
  T8_ASSERT (element_is_valid (elem));
  return P8EST_FACES;
}

int
t8_default_scheme_hex_c::element_get_max_num_faces (const t8_element_t *elem) const
{
  return P8EST_FACES;
}

int
t8_default_scheme_hex_c::element_get_num_children (const t8_element_t *elem) const
{
  T8_ASSERT (element_is_valid (elem));
  return P8EST_CHILDREN;
}

int
t8_default_scheme_hex_c::element_get_num_face_children (const t8_element_t *elem, int face) const
{
  T8_ASSERT (element_is_valid (elem));
  return 4;
}

int
t8_default_scheme_hex_c::element_get_face_corner (const t8_element_t *element, int face, int corner) const
{
  T8_ASSERT (element_is_valid (element));
  T8_ASSERT (0 <= face && face < P8EST_FACES);
  T8_ASSERT (0 <= corner && corner < 4);

  return p8est_face_corners[face][corner];
}

void
t8_default_scheme_hex_c::element_get_child (const t8_element_t *elem, int childid, t8_element_t *child) const
{
  const p8est_quadrant_t *q = (const p8est_quadrant_t *) elem;
  const p4est_qcoord_t shift = P8EST_QUADRANT_LEN (q->level + 1);
  p8est_quadrant_t *r = (p8est_quadrant_t *) child;

  T8_ASSERT (element_is_valid (elem));
  T8_ASSERT (element_is_valid (child));
  T8_ASSERT (p8est_quadrant_is_extended (q));
  T8_ASSERT (q->level < HEX_REFINE_MAXLEVEL);
  T8_ASSERT (0 <= childid && childid < P8EST_CHILDREN);

  r->x = childid & 0x01 ? (q->x | shift) : q->x;
  r->y = childid & 0x02 ? (q->y | shift) : q->y;
  r->z = childid & 0x04 ? (q->z | shift) : q->z;
  r->level = q->level + 1;
  if (q != r) {
    T8_ASSERT (p8est_quadrant_is_parent (q, r));
  }
}

void
t8_default_scheme_hex_c::element_get_children (const t8_element_t *elem, int length, t8_element_t *c[]) const
{
  T8_ASSERT (element_is_valid (elem));
#ifdef T8_ENABLE_DEBUG
  {
    int i;
    for (i = 0; i < P8EST_CHILDREN; i++) {
      T8_ASSERT (element_is_valid (c[i]));
    }
  }
#endif
  T8_ASSERT (length == P8EST_CHILDREN);

  p8est_quadrant_childrenpv ((const p8est_quadrant_t *) elem, (p8est_quadrant_t **) c);
}

int
t8_default_scheme_hex_c::element_get_child_id (const t8_element_t *elem) const
{
  T8_ASSERT (element_is_valid (elem));
  return p8est_quadrant_child_id ((const p8est_quadrant_t *) elem);
}

int
t8_default_scheme_hex_c::element_get_ancestor_id (const t8_element_t *elem, int level) const
{
  return p8est_quadrant_ancestor_id ((p8est_quadrant_t *) elem, level);
}

int
t8_default_scheme_hex_c::elements_are_family (t8_element_t *const *fam) const
{
#ifdef T8_ENABLE_DEBUG
  {
    int i;
    for (i = 0; i < P8EST_CHILDREN; i++) {
      T8_ASSERT (element_is_valid (fam[i]));
    }
  }
#endif
  return p8est_quadrant_is_familypv ((p8est_quadrant_t **) fam);
}

void
t8_default_scheme_hex_c::element_get_nca (const t8_element_t *elem1, const t8_element_t *elem2, t8_element_t *nca) const
{
  T8_ASSERT (element_is_valid (elem1));
  T8_ASSERT (element_is_valid (elem2));
  p8est_nearest_common_ancestor ((const p8est_quadrant_t *) elem1, (const p8est_quadrant_t *) elem2,
                                 (p8est_quadrant_t *) nca);
}

t8_element_shape_t
t8_default_scheme_hex_c::element_get_face_shape (const t8_element_t *elem, int face) const
{
  T8_ASSERT (element_is_valid (elem));
  return T8_ECLASS_QUAD;
}

void
t8_default_scheme_hex_c::element_get_children_at_face (const t8_element_t *elem, int face, t8_element_t *children[],
                                                       int num_children, int *child_indices) const
{
  int child_ids_local[4], i, *child_ids;

  T8_ASSERT (element_is_valid (elem));
#ifdef T8_ENABLE_DEBUG
  {
    int j;
    for (j = 0; j < P4EST_CHILDREN; j++) {
      T8_ASSERT (element_is_valid (children[j]));
    }
  }
#endif
  T8_ASSERT (0 <= face && face < P8EST_FACES);
  T8_ASSERT (num_children == element_get_num_face_children (elem, face));

  if (child_indices != NULL) {
    child_ids = child_indices;
  }
  else {
    child_ids = child_ids_local;
  }
  /*
   * Compute the child id of the first and second child at the face.
   *
   * The faces of the quadrant are enumerated like this:
   *
   *          f_3
   *       x ---- x
   *      /  f_5 /|          z y
   *     x ---- x |          |/
   * f_0 |      | x f_1       -- x
   *     |  f_2 |/
   *     x ---- x
   *        f_4
   */
  for (i = 0; i < P8EST_HALF; ++i) {
    child_ids[i] = p8est_face_corners[face][i];
  }

  /* Create the four face children */
  /* We have to revert the order and compute the zeroth child last, since
   * the usage allows for elem == children[0].
   */
  for (i = 3; i >= 0; i--) {
    element_get_child (elem, child_ids[i], children[i]);
  }
}

int
t8_default_scheme_hex_c::element_face_get_child_face (const t8_element_t *elem, int face, int face_child) const
{
  T8_ASSERT (element_is_valid (elem));
  /* For octants the face enumeration of children is the same as for the parent. */
  return face;
}

int
t8_default_scheme_hex_c::element_face_get_parent_face (const t8_element_t *elem, int face) const
{
  int child_id;
  const p8est_quadrant_t *q = (const p8est_quadrant_t *) elem;

  T8_ASSERT (element_is_valid (elem));
  if (q->level == 0) {
    return face;
  }
  /* Determine whether face is a subface of the parent.
   * This is the case if the child_id matches one of the faces corners */
  child_id = p8est_quadrant_child_id (q);
  if (child_id == p8est_face_corners[face][0] || child_id == p8est_face_corners[face][1]
      || child_id == p8est_face_corners[face][2] || child_id == p8est_face_corners[face][3]) {
    return face;
  }
  return -1;
}

int
t8_default_scheme_hex_c::element_get_tree_face (const t8_element_t *elem, int face) const
{
  T8_ASSERT (element_is_valid (elem));
  T8_ASSERT (0 <= face && face < P8EST_FACES);
  /* For hexahedra the face and the tree face number are the same. */
  return face;
}

int
t8_default_scheme_hex_c::element_extrude_face (const t8_element_t *face, const t8_scheme *face_scheme,
                                               t8_element_t *elem, int root_face) const
{
  const p4est_quadrant_t *b = (const p4est_quadrant_t *) face;
  p8est_quadrant_t *q = (p8est_quadrant_t *) elem;

  T8_ASSERT (T8_COMMON_IS_TYPE (face_scheme, const t8_default_scheme_quad_c *));
  T8_ASSERT (element_is_valid (elem));
  T8_ASSERT (face_scheme->eclass == T8_ECLASS_QUAD);
  T8_ASSERT (face_scheme->element_is_valid (face));
  T8_ASSERT (0 <= root_face && root_face < P8EST_FACES);
  q->level = b->level;
  /*
   * The faces of the root quadrant are enumerated like this:
   *
   *       x ---- x
   *      /  f_5 /|
   *     x ---- x |
   * f_0 |      | x f_1
   *     |  f_2 |/
   *     x ---- x
   *        f_4
   *
   * We need to rescale the coordinates since a quadrant may have a different
   * root length than an octant.
   */
  switch (root_face) {
  case 0:
    q->x = 0;
    q->y = ((int64_t) b->x * P8EST_ROOT_LEN) / P4EST_ROOT_LEN;
    q->z = ((int64_t) b->y * P8EST_ROOT_LEN) / P4EST_ROOT_LEN;
    break;
  case 1:
    q->x = P8EST_LAST_OFFSET (q->level);
    q->y = ((int64_t) b->x * P8EST_ROOT_LEN) / P4EST_ROOT_LEN;
    q->z = ((int64_t) b->y * P8EST_ROOT_LEN) / P4EST_ROOT_LEN;
    break;
  case 2:
    q->x = ((int64_t) b->x * P8EST_ROOT_LEN) / P4EST_ROOT_LEN;
    q->y = 0;
    q->z = ((int64_t) b->y * P8EST_ROOT_LEN) / P4EST_ROOT_LEN;
    break;
  case 3:
    q->x = ((int64_t) b->x * P8EST_ROOT_LEN) / P4EST_ROOT_LEN;
    q->y = P8EST_LAST_OFFSET (q->level);
    q->z = ((int64_t) b->y * P8EST_ROOT_LEN) / P4EST_ROOT_LEN;
    break;
  case 4:
    q->x = ((int64_t) b->x * P8EST_ROOT_LEN) / P4EST_ROOT_LEN;
    q->y = ((int64_t) b->y * P8EST_ROOT_LEN) / P4EST_ROOT_LEN;
    q->z = 0;
    break;
  case 5:
    q->x = ((int64_t) b->x * P8EST_ROOT_LEN) / P4EST_ROOT_LEN;
    q->y = ((int64_t) b->y * P8EST_ROOT_LEN) / P4EST_ROOT_LEN;
    q->z = P8EST_LAST_OFFSET (q->level);
    break;
  }
  /* We return the face of q at which we extruded. This is the same number
   * as root_face. */
  return root_face;
}

/** Construct the first descendant of an element that touches a given face. */
void
t8_default_scheme_hex_c::element_construct_first_descendant_face (const t8_element_t *elem, int face,
                                                                  t8_element_t *first_desc, int level) const
{
  const p8est_quadrant_t *q = (const p8est_quadrant_t *) elem;
  p8est_quadrant_t *desc = (p8est_quadrant_t *) first_desc;
  int first_face_corner;

  T8_ASSERT (0 <= face && face < P8EST_FACES);
  T8_ASSERT (0 <= level && level <= HEX_REFINE_MAXLEVEL);

  /* Get the first corner of q that belongs to face */
  first_face_corner = p8est_face_corners[face][0];
  /* Construct the descendant of q in this corner */
  p8est_quadrant_corner_descendant (q, desc, first_face_corner, level);
}

/** Construct the last descendant of an element that touches a given face. */
void
t8_default_scheme_hex_c::element_construct_last_descendant_face (const t8_element_t *elem, int face,
                                                                 t8_element_t *last_desc, int level) const
{
  const p8est_quadrant_t *q = (const p8est_quadrant_t *) elem;
  p8est_quadrant_t *desc = (p8est_quadrant_t *) last_desc;
  int last_face_corner;

  T8_ASSERT (0 <= face && face < P8EST_FACES);
  T8_ASSERT (0 <= level && level <= HEX_REFINE_MAXLEVEL);

  /* Get the last corner of q that belongs to face */
  last_face_corner = p8est_face_corners[face][3];
  /* Construct the descendant of q in this corner */
  p8est_quadrant_corner_descendant (q, desc, last_face_corner, level);
}

void
t8_default_scheme_hex_c::element_construct_boundary_face (const t8_element_t *elem, int face, t8_element_t *boundary,
                                                          const t8_scheme *boundary_scheme) const
{
  const p8est_quadrant_t *q = (const p8est_quadrant_t *) elem;
  p4est_quadrant_t *b = (p4est_quadrant_t *) boundary;

  T8_ASSERT (element_is_valid (elem));
  T8_ASSERT (T8_COMMON_IS_TYPE (boundary_scheme, const t8_default_scheme_quad_c *));
  T8_ASSERT (boundary_scheme->eclass == T8_ECLASS_QUAD);
  T8_ASSERT (boundary_scheme->element_is_valid (boundary));
  T8_ASSERT (0 <= face && face < P8EST_FACES);

  /* The level of the boundary element is the same as the quadrant's level */
  b->level = q->level;
  /*
   * The faces of the quadrant are enumerated like this:
   *
   *       x ---- x
   *      /  f_5 /|
   *     x ---- x |
   * f_0 |      | x f_1
   *     |  f_2 |/
   *     x ---- x
   *        f_4
   *
   * If face = 0 or face = 1 then b->x = q->y, b->y = q->z
   * if face = 2 or face = 3 then b->x = q->x, b->y = q->z
   * if face = 4 or face = 5 then b->x = q->x, b->y = q->y
   *
   * We have to scale the coordinates since a root quadrant may have
   * different length than a root hex.
   */
  b->x = (face >> 1 ? q->x : q->y) * ((t8_linearidx_t) P4EST_ROOT_LEN / P8EST_ROOT_LEN); /* true if face >= 2 */
  b->y = (face >> 2 ? q->y : q->z) * ((t8_linearidx_t) P4EST_ROOT_LEN / P8EST_ROOT_LEN); /* true if face >= 4 */
  T8_ASSERT (!p8est_quadrant_is_extended (q) || p4est_quadrant_is_extended (b));
}

int
t8_default_scheme_hex_c::element_is_root_boundary (const t8_element_t *elem, int face) const
{
  const p8est_quadrant_t *q = (const p8est_quadrant_t *) elem;
  p4est_qcoord_t coord;

  T8_ASSERT (element_is_valid (elem));
  T8_ASSERT (0 <= face && face < P8EST_FACES);

  /* if face is 0 or 1 q->x
   *            2 or 3 q->y
   *            4 or 5 q->z
   */
  coord = face >> 2 ? q->z : face >> 1 ? q->y : q->x;
  /* If face is 0,2 or 4 check against 0.
   * If face is 1,3 or 5 check against LAST_OFFSET */
  return coord == (face & 1 ? P8EST_LAST_OFFSET (q->level) : 0);
}

int
t8_default_scheme_hex_c::element_construct_face_neighbor_inside (const t8_element_t *elem, t8_element_t *neigh,
                                                                 int face, int *neigh_face) const
{
  const p8est_quadrant_t *q = (const p8est_quadrant_t *) elem;
  p8est_quadrant_t *n = (p8est_quadrant_t *) neigh;

  T8_ASSERT (element_is_valid (elem));
  T8_ASSERT (element_is_valid (neigh));
  T8_ASSERT (0 <= face && face < P8EST_FACES);
  /* Compute the face neighbor */
  p8est_quadrant_face_neighbor (q, face, n);

  /* Compute the face of q that coincides with face.
   * face   neigh_face    face      neigh_face
   *   0        1           4           5
   *   1        0           5           4
   *   2        3
   *   3        2
   */
  T8_ASSERT (neigh_face != NULL);
  *neigh_face = p8est_face_dual[face];
  /* return true if neigh is inside the root */
  return p8est_quadrant_is_inside_root (n);
}

void
t8_default_scheme_hex_c::element_set_linear_id (t8_element_t *elem, int level, t8_linearidx_t id) const
{
  T8_ASSERT (element_is_valid (elem));
  T8_ASSERT (0 <= level && level <= HEX_LINEAR_MAXLEVEL);
  T8_ASSERT (0 <= id && id < ((t8_linearidx_t) 1) << P8EST_DIM * level);

  p8est_quadrant_set_morton ((p8est_quadrant_t *) elem, level, id);
}

t8_linearidx_t
t8_default_scheme_hex_c::element_get_linear_id (const t8_element_t *elem, int level) const
{
  T8_ASSERT (element_is_valid (elem));
  T8_ASSERT (0 <= level && level <= HEX_LINEAR_MAXLEVEL);

  return p8est_quadrant_linear_id ((p8est_quadrant_t *) elem, level);
}

void
t8_default_scheme_hex_c::element_construct_first_descendant (const t8_element_t *elem, t8_element_t *desc,
                                                             int level) const
{
  T8_ASSERT (element_is_valid (elem));
  T8_ASSERT (element_is_valid (desc));
  T8_ASSERT (0 <= level && level <= HEX_REFINE_MAXLEVEL);
  p8est_quadrant_first_descendant ((p8est_quadrant_t *) elem, (p8est_quadrant_t *) desc, level);
}

void
t8_default_scheme_hex_c::element_construct_last_descendant (const t8_element_t *elem, t8_element_t *desc,
                                                            int level) const
{
  T8_ASSERT (element_is_valid (elem));
  T8_ASSERT (element_is_valid (desc));
  T8_ASSERT (0 <= level && level <= HEX_REFINE_MAXLEVEL);
  p8est_quadrant_last_descendant ((p8est_quadrant_t *) elem, (p8est_quadrant_t *) desc, level);
}

void
t8_default_scheme_hex_c::element_construct_successor (const t8_element_t *elem1, t8_element_t *elem2) const
{
  T8_ASSERT (element_is_valid (elem1));
  T8_ASSERT (element_is_valid (elem2));
  T8_ASSERT (0 <= element_get_level (elem1) && element_get_level (elem1) <= HEX_REFINE_MAXLEVEL);
  p8est_quadrant_successor ((p8est_quadrant_t *) elem1, (p8est_quadrant_t *) elem2);
}

void
t8_default_scheme_hex_c::element_get_anchor (const t8_element_t *elem, int coord[3]) const
{
  p8est_quadrant_t *q;

  T8_ASSERT (element_is_valid (elem));
  q = (p8est_quadrant_t *) elem;
  coord[0] = q->x;
  coord[1] = q->y;
  coord[2] = q->z;
}

void
t8_default_scheme_hex_c::element_get_vertex_integer_coords (const t8_element_t *elem, int vertex, int coords[]) const
{
  const p8est_quadrant_t *q1 = (const p8est_quadrant_t *) elem;

  T8_ASSERT (0 <= vertex && vertex < 8);
  /* Get the length of the quadrant */
  const int len = P8EST_QUADRANT_LEN (q1->level);
  /* Compute the x, y and z coordinates of the vertex depending on the
   * vertex number */
  coords[0] = q1->x + (vertex & 1 ? 1 : 0) * len;
  coords[1] = q1->y + (vertex & 2 ? 1 : 0) * len;
  coords[2] = q1->z + (vertex & 4 ? 1 : 0) * len;
}

void
t8_default_scheme_hex_c::element_get_vertex_reference_coords (const t8_element_t *elem, const int vertex,
                                                              double coords[]) const
{
  T8_ASSERT (element_is_valid (elem));
  T8_ASSERT (0 <= vertex && vertex < 8);

  int coords_int[3];
  element_get_vertex_integer_coords (elem, vertex, coords_int);

  /* We divide the integer coordinates by the root length of the hex
   * to obtain the reference coordinates. */
  coords[0] = coords_int[0] / (double) P8EST_ROOT_LEN;
  coords[1] = coords_int[1] / (double) P8EST_ROOT_LEN;
  coords[2] = coords_int[2] / (double) P8EST_ROOT_LEN;
}

void
t8_default_scheme_hex_c::element_get_reference_coords (const t8_element_t *elem, const double *ref_coords,
                                                       const size_t num_coords, double *out_coords) const
{
  T8_ASSERT (element_is_valid (elem));
  t8_dhex_compute_reference_coords ((const t8_dhex_t *) elem, ref_coords, num_coords, out_coords);
}

int
t8_default_scheme_hex_c::refines_irregular () const
{
  /* Hex always refine regularly */
  return 0;
}

void
t8_default_scheme_hex_c::element_new (int length, t8_element_t **elem) const
{
  /* allocate memory for a hex */
  t8_default_scheme_common_c::element_new (length, elem);

  /* in debug mode, set sensible default values. */
#ifdef T8_ENABLE_DEBUG
  {
    int i;
    for (i = 0; i < length; i++) {
      get_root (elem[i]);
      T8_QUAD_SET_TDIM ((p8est_quadrant_t *) elem[i], 3);
    }
  }
#endif
}

void
t8_default_scheme_hex_c::element_init (int length, t8_element_t *elem) const
{
#ifdef T8_ENABLE_DEBUG
  p8est_quadrant_t *quads = (p8est_quadrant_t *) elem;
  for (int i = 0; i < length; i++) {
    p8est_quadrant_set_morton (quads + i, 0, 0);
    T8_QUAD_SET_TDIM (quads + i, 3);
    T8_ASSERT (p8est_quadrant_is_extended (quads + i));
  }
#endif
}

#ifdef T8_ENABLE_DEBUG
int
t8_default_scheme_hex_c::element_is_valid (const t8_element_t *elem) const
{
  /* TODO: additional checks? do we set pad8 or similar?
   */
  return p8est_quadrant_is_extended ((const p8est_quadrant_t *) elem)
         && T8_QUAD_GET_TDIM ((const p8est_quadrant_t *) elem) == 3;
}

void
t8_default_scheme_hex_c::element_to_string (const t8_element_t *elem, char *debug_string, const int string_size) const
{
  T8_ASSERT (element_is_valid (elem));
  T8_ASSERT (debug_string != NULL);
  p8est_quadrant_t *hex = (p8est_quadrant_t *) elem;
  snprintf (debug_string, string_size, "x: %i, y: %i, z: %i, level: %i", hex->x, hex->y, hex->z, hex->level);
}
#endif

/* Constructor */
t8_default_scheme_hex_c::t8_default_scheme_hex_c (void)
{
  eclass = T8_ECLASS_HEX;
  element_size = sizeof (t8_phex_t);
  ts_context = sc_mempool_new (element_size);
}

t8_default_scheme_hex_c::~t8_default_scheme_hex_c ()
{
  /* This destructor is empty since the destructor of the
   * default_common scheme is called automatically and it
   * suffices to destroy the quad_scheme.
   * However we need to provide an implementation of the destructor
   * and hence this empty function. */
}
void
t8_default_scheme_hex_c::get_root (t8_element_t *elem) const
{
  p8est_quadrant_t *hex = (p8est_quadrant_t *) elem;
  p8est_quadrant_set_morton (hex, 0, 0);
  T8_ASSERT (p8est_quadrant_is_extended (hex));
}

/* each hex is packed as x,y,z coordinates and the level */
void
t8_default_scheme_hex_c::element_MPI_Pack (t8_element_t **const elements, const unsigned int count, void *send_buffer,
                                           const int buffer_size, int *position, sc_MPI_Comm comm) const
{
  int mpiret;
  p8est_quadrant_t **quads = (p8est_quadrant_t **) elements;
  for (unsigned int ielem = 0; ielem < count; ielem++) {
    mpiret = sc_MPI_Pack (&(quads[ielem]->x), 1, sc_MPI_INT, send_buffer, buffer_size, position, comm);
    SC_CHECK_MPI (mpiret);
    mpiret = sc_MPI_Pack (&quads[ielem]->y, 1, sc_MPI_INT, send_buffer, buffer_size, position, comm);
    SC_CHECK_MPI (mpiret);
    mpiret = sc_MPI_Pack (&quads[ielem]->z, 1, sc_MPI_INT, send_buffer, buffer_size, position, comm);
    SC_CHECK_MPI (mpiret);
    mpiret = sc_MPI_Pack (&quads[ielem]->level, 1, sc_MPI_INT8_T, send_buffer, buffer_size, position, comm);
    SC_CHECK_MPI (mpiret);
  }
}

/* each hex is packed as x,y,z coordinates and the level */
void
t8_default_scheme_hex_c::element_MPI_Pack_size (const unsigned int count, sc_MPI_Comm comm, int *pack_size) const
{
  int singlesize = 0;
  int datasize = 0;
  int mpiret;

  /* x,y,z */
  mpiret = sc_MPI_Pack_size (1, sc_MPI_INT, comm, &datasize);
  SC_CHECK_MPI (mpiret);
  singlesize += 3 * datasize;

  /* level */
  mpiret = sc_MPI_Pack_size (1, sc_MPI_INT8_T, comm, &datasize);
  SC_CHECK_MPI (mpiret);
  singlesize += datasize;

  *pack_size = count * singlesize;
}

/* each hex is packed as x,y,z coordinates and the level */
void
t8_default_scheme_hex_c::element_MPI_Unpack (void *recvbuf, const int buffer_size, int *position,
                                             t8_element_t **elements, const unsigned int count, sc_MPI_Comm comm) const
{
  int mpiret;
  p8est_quadrant_t **quads = (p8est_quadrant_t **) elements;
  for (unsigned int ielem = 0; ielem < count; ielem++) {
    mpiret = sc_MPI_Unpack (recvbuf, buffer_size, position, &(quads[ielem]->x), 1, sc_MPI_INT, comm);
    SC_CHECK_MPI (mpiret);
    mpiret = sc_MPI_Unpack (recvbuf, buffer_size, position, &(quads[ielem]->y), 1, sc_MPI_INT, comm);
    SC_CHECK_MPI (mpiret);
    mpiret = sc_MPI_Unpack (recvbuf, buffer_size, position, &(quads[ielem]->z), 1, sc_MPI_INT, comm);
    SC_CHECK_MPI (mpiret);
    mpiret = sc_MPI_Unpack (recvbuf, buffer_size, position, &(quads[ielem]->level), 1, sc_MPI_INT8_T, comm);
    SC_CHECK_MPI (mpiret);
  }
}

T8_EXTERN_C_END ();
