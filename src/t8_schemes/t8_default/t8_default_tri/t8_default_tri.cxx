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

#include <t8_schemes/t8_default/t8_default_common/t8_default_common.hxx>
#include <t8_schemes/t8_default/t8_default_tri/t8_default_tri.hxx>
#include <t8_schemes/t8_default/t8_default_tri/t8_dtri_bits.h>
#include <t8_schemes/t8_default/t8_default_line/t8_dline_bits.h>
#include <t8_schemes/t8_default/t8_default_tet/t8_dtet.h>
#include <t8_schemes/t8_default/t8_default_tri/t8_dtri_connectivity.h>
#include <t8_schemes/t8_scheme.hxx>

/* We want to export the whole implementation to be callable from "C" */
T8_EXTERN_C_BEGIN ();

size_t
t8_default_scheme_tri::get_element_size (void) const
{
  return sizeof (t8_dtri_t);
}

int
t8_default_scheme_tri::get_maxlevel (void) const
{
  return T8_DTRI_MAXLEVEL;
}

int
t8_default_scheme_tri::element_get_level (const t8_element_t *elem) const
{
  T8_ASSERT (element_is_valid (elem));
  return t8_dtri_get_level ((t8_dtri_t *) elem);
}

void
t8_default_scheme_tri::element_copy (const t8_element_t *source, t8_element_t *dest) const
{
  T8_ASSERT (element_is_valid (source));
  T8_ASSERT (element_is_valid (dest));
  t8_dtri_copy ((const t8_dtri_t *) source, (t8_dtri_t *) dest);
}

int
t8_default_scheme_tri::element_compare (const t8_element_t *elem1, const t8_element_t *elem2) const
{
  T8_ASSERT (element_is_valid (elem1));
  T8_ASSERT (element_is_valid (elem2));

  return t8_dtri_compare ((const t8_dtri_t *) elem1, (const t8_dtri_t *) elem2);
}

int
t8_default_scheme_tri::element_is_equal (const t8_element_t *elem1, const t8_element_t *elem2) const
{
  return t8_dtri_equal ((const t8_dtri_t *) elem1, (const t8_dtri_t *) elem2);
}

void
t8_default_scheme_tri::element_get_parent (const t8_element_t *elem, t8_element_t *parent) const
{
  const t8_dtri_t *t = (const t8_dtri_t *) elem;
  t8_dtri_t *p = (t8_dtri_t *) parent;

  T8_ASSERT (element_is_valid (elem));
  T8_ASSERT (element_is_valid (parent));
  t8_dtri_parent (t, p);
}

void
t8_default_scheme_tri::element_get_sibling (const t8_element_t *elem, int sibid, t8_element_t *sibling) const
{
  const t8_dtri_t *t = (const t8_dtri_t *) elem;
  t8_dtri_t *s = (t8_dtri_t *) sibling;

  T8_ASSERT (element_is_valid (elem));
  t8_dtri_sibling (t, sibid, s);
}

int
t8_default_scheme_tri::element_get_num_faces ([[maybe_unused]] const t8_element_t *elem) const
{
  T8_ASSERT (element_is_valid (elem));
  return T8_DTRI_FACES;
}

int
t8_default_scheme_tri::element_get_max_num_faces ([[maybe_unused]] const t8_element_t *elem) const
{
  return T8_DTRI_FACES;
}

int
t8_default_scheme_tri::element_get_num_children ([[maybe_unused]] const t8_element_t *elem) const
{
  T8_ASSERT (element_is_valid (elem));
  return T8_DTRI_CHILDREN;
}

int
t8_default_scheme_tri::element_get_num_face_children ([[maybe_unused]] const t8_element_t *elem, [[maybe_unused]] int face) const
{
  T8_ASSERT (element_is_valid (elem));
  return T8_DTRI_FACE_CHILDREN;
}

int
t8_default_scheme_tri::element_get_face_corner ([[maybe_unused]] const t8_element_t *element, int face, int corner) const
{
  T8_ASSERT (element_is_valid (element));
  T8_ASSERT (0 <= face && face < T8_DTRI_FACES);
  T8_ASSERT (0 <= corner && corner < 2);
  return t8_dtri_face_corner[face][corner];
}

int
t8_default_scheme_tri::element_get_corner_face ([[maybe_unused]] const t8_element_t *element, int corner, int face) const
{
  T8_ASSERT (element_is_valid (element));
  T8_ASSERT (0 <= corner && corner < T8_DTRI_CORNERS);
  T8_ASSERT (0 <= face && face < 2);
  return t8_dtri_corner_face[corner][face];
}

void
t8_default_scheme_tri::element_get_child (const t8_element_t *elem, int childid, t8_element_t *child) const
{
  const t8_dtri_t *t = (const t8_dtri_t *) elem;
  t8_dtri_t *c = (t8_dtri_t *) child;

  T8_ASSERT (element_is_valid (elem));
  T8_ASSERT (element_is_valid (child));
  t8_dtri_child (t, childid, c);
}

void
t8_default_scheme_tri::element_get_children (const t8_element_t *elem, [[maybe_unused]] int length, t8_element_t *c[]) const
{
  T8_ASSERT (element_is_valid (elem));
#ifdef T8_ENABLE_DEBUG
  {
    for (int j = 0; j < length; j++) {
      T8_ASSERT (element_is_valid (c[j]));
    }
  }
#endif
  T8_ASSERT (length == T8_DTRI_CHILDREN);

  t8_dtri_childrenpv ((const t8_dtri_t *) elem, (t8_dtri_t **) c);
}

int
t8_default_scheme_tri::element_get_child_id (const t8_element_t *elem) const
{
  T8_ASSERT (element_is_valid (elem));
  return t8_dtri_child_id ((t8_dtri_t *) elem);
}

int
t8_default_scheme_tri::element_get_ancestor_id (const t8_element_t *elem, int level) const
{
  return t8_dtri_ancestor_id ((t8_dtri_t *) elem, level);
}

int
t8_default_scheme_tri::elements_are_family (t8_element_t *const *fam) const
{
#ifdef T8_ENABLE_DEBUG
  {
    for (int j = 0; j < T8_DTRI_CHILDREN; j++) {
      T8_ASSERT (element_is_valid (fam[j]));
    }
  }
#endif
  return t8_dtri_is_familypv ((const t8_dtri_t **) fam);
}

void
t8_default_scheme_tri::element_get_nca (const t8_element_t *elem1, const t8_element_t *elem2, t8_element_t *nca) const
{
  const t8_dtri_t *t1 = (const t8_dtri_t *) elem1;
  const t8_dtri_t *t2 = (const t8_dtri_t *) elem2;
  t8_dtri_t *c = (t8_dtri_t *) nca;

  T8_ASSERT (element_is_valid (elem1));
  T8_ASSERT (element_is_valid (elem2));
  t8_dtri_nearest_common_ancestor (t1, t2, c);
}

t8_element_shape_t
t8_default_scheme_tri::element_get_face_shape ([[maybe_unused]] const t8_element_t *elem, [[maybe_unused]] int face) const
{
  T8_ASSERT (element_is_valid (elem));
  return T8_ECLASS_LINE;
}

void
t8_default_scheme_tri::element_get_children_at_face (const t8_element_t *elem, int face, t8_element_t *children[],
                                                     int num_children, int *child_indices) const
{
  const t8_dtri_t *t = (const t8_dtri_t *) elem;
  t8_dtri_t **c = (t8_dtri_t **) children;

  T8_ASSERT (element_is_valid (elem));
#ifdef T8_ENABLE_DEBUG
  {
    for (int j = 0; j < num_children; j++) {
      T8_ASSERT (element_is_valid (children[j]));
    }
  }
#endif
  T8_ASSERT (0 <= face && face < T8_DTRI_FACES);
  T8_ASSERT (num_children == T8_DTRI_FACE_CHILDREN);

  t8_dtri_children_at_face (t, face, c, num_children, child_indices);
}

int
t8_default_scheme_tri::element_face_get_child_face (const t8_element_t *elem, int face, int face_child) const
{
  T8_ASSERT (element_is_valid (elem));
  T8_ASSERT (0 <= face && face < T8_DTRI_FACES);
  T8_ASSERT (0 <= face_child && face_child < T8_DTRI_FACE_CHILDREN);
  return t8_dtri_face_child_face ((const t8_dtri_t *) elem, face, face_child);
}

int
t8_default_scheme_tri::element_face_get_parent_face (const t8_element_t *elem, int face) const
{
  T8_ASSERT (0 <= face && face < T8_DTRI_FACES);

  return t8_dtri_face_parent_face ((const t8_dtri_t *) elem, face);
}

int
t8_default_scheme_tri::element_get_tree_face (const t8_element_t *elem, int face) const
{
  T8_ASSERT (element_is_valid (elem));
  T8_ASSERT (0 <= face && face < T8_DTRI_FACES);
  return t8_dtri_tree_face ((t8_dtri_t *) elem, face);
}

void
t8_default_scheme_tri::element_transform_face (const t8_element_t *elem1, t8_element_t *elem2, int orientation,
                                               int sign, int is_smaller_face) const
{
  T8_ASSERT (element_is_valid (elem1));
  T8_ASSERT (element_is_valid (elem2));
  t8_dtri_transform_face ((const t8_dtri_t *) elem1, (t8_dtri_t *) elem2, orientation, sign, is_smaller_face);
}

/* Construct the inner element from a boundary element. */
/* This function is defined here instead of in t8_dri_bits.c since
 * the compile logic does not allow for t8_dtri_t and t8_dtet_t to exist
 * both in t8_dtri_bits.c. This would be needed by an implementation, at least
 * for tets. */
int
t8_default_scheme_tri::element_extrude_face (const t8_element_t *face, t8_element_t *elem, int root_face,
    [[maybe_unused]] const t8_scheme *scheme) const
{
  const t8_dline_t *l = (const t8_dline_t *) face;
  t8_dtri_t *t = (t8_dtri_t *) elem;

  T8_ASSERT (element_is_valid (elem));
  T8_ASSERT (scheme->element_is_valid (T8_ECLASS_LINE, face));
  T8_ASSERT (0 <= root_face && root_face < T8_DTRI_FACES);
  /*
   * The faces of the root triangle are enumerated like this
   *
   *
   *         x
   *    f_1 /|
   *       / | f_0
   *      x--x
   *      f_2
   *
   * Boundary triangles are always of type 0.
   * We have to scale the coordinates since the root triangle may
   * have a different scale than the root line.
   */
  t->level = l->level;
  t->type = 0;
  switch (root_face) {
  case 0:
    t->x = T8_DTRI_ROOT_LEN - T8_DTRI_LEN (t->level);
    t->y = ((int64_t) l->x * T8_DTRI_ROOT_LEN) / T8_DLINE_ROOT_LEN;
    break;
  case 1:
    t->x = ((int64_t) l->x * T8_DTRI_ROOT_LEN) / T8_DLINE_ROOT_LEN;
    t->y = ((int64_t) l->x * T8_DTRI_ROOT_LEN) / T8_DLINE_ROOT_LEN;
    break;
  case 2:
    t->x = ((int64_t) l->x * T8_DTRI_ROOT_LEN) / T8_DLINE_ROOT_LEN;
    t->y = 0;
    break;
  default:
    SC_ABORT_NOT_REACHED ();
  }
  /* We return the face number of t of the face at which we extruded.
   * Since in all cases t has type 0, the face number coincides with
   * the root face number. */
  return t8_dtri_root_face_to_face (t, root_face);
}

void
t8_default_scheme_tri::element_get_first_descendant_face (const t8_element_t *elem, int face, t8_element_t *first_desc,
                                                          int level) const
{
  int corner;
  T8_ASSERT (0 <= face && face < T8_DTRI_FACES);
  T8_ASSERT (0 <= level && level <= T8_DTRI_MAXLEVEL);

  /* Compute the first corner of this face */
  corner = t8_dtri_face_corner[face][0];
  /* Compute the descendant in this corner */
  t8_dtri_corner_descendant ((const t8_dtri_t *) elem, (t8_dtri_t *) first_desc, corner, level);
}

void
t8_default_scheme_tri::element_get_last_descendant_face (const t8_element_t *elem, int face, t8_element_t *last_desc,
                                                         int level) const
{
  T8_ASSERT (0 <= face && face < T8_DTRI_FACES);
  T8_ASSERT (0 <= level && level <= T8_DTRI_MAXLEVEL);
  /* Compute the last corner of this face */
  const int corner = t8_dtri_face_corner[face][1];
  /* Compute the descendant in this corner */
  t8_dtri_corner_descendant ((const t8_dtri_t *) elem, (t8_dtri_t *) last_desc, corner, level);
}

/* Construct the boundary element at a specific face. */
void
t8_default_scheme_tri::element_get_boundary_face (const t8_element_t *elem, int face, t8_element_t *boundary,
    [[maybe_unused]] const t8_scheme *scheme) const
{
  const t8_dtri_t *t = (const t8_dtri_t *) elem;
  t8_dline_t *l = (t8_dline_t *) boundary;

  T8_ASSERT (element_is_valid (elem));
  T8_ASSERT (scheme->element_is_valid (T8_ECLASS_LINE, boundary));
  T8_ASSERT (0 <= face && face < T8_DTRI_FACES);
  /* The level of the boundary element is the same as the quadrant's level */
  l->level = t->level;
  /*
   * The faces of the triangle are enumerated like this
   *        type 0                    type 1
   *
   *                                    f_0
   *         x                         x--x
   *    f_1 /|                         | /
   *       / | f_0                 f_2 |/ f_1
   *      x--x                         x
   *      f_2
   *
   *  If face = 0, type = 0 or face = 2, type = 1 then l->x = t->y
   *  If face = 1, face = 2, type = 0 or
   *     face = 0, face = 1, type = 1             then l->x = t->x
   */
  if ((face == 0 && t->type == 0) || (face == 2 && t->type == 1)) {
    l->x = t->y * T8_DLINE_ROOT_BY_DTRI_ROOT;
  }
  else {
    l->x = t->x * T8_DLINE_ROOT_BY_DTRI_ROOT;
  }
  /* TODO: Take the level into account! */
}

int
t8_default_scheme_tri::element_is_root_boundary (const t8_element_t *elem, int face) const
{
  const t8_dtri_t *t = (const t8_dtri_t *) elem;
  T8_ASSERT (element_is_valid (elem));

  return t8_dtri_is_root_boundary (t, face);
}

int
t8_default_scheme_tri::element_get_face_neighbor_inside (const t8_element_t *elem, t8_element_t *neigh, int face,
                                                         int *neigh_face) const
{
  const t8_dtri_t *t = (const t8_dtri_t *) elem;
  t8_dtri_t *n = (t8_dtri_t *) neigh;

  T8_ASSERT (element_is_valid (elem));
  T8_ASSERT (element_is_valid (neigh));
  T8_ASSERT (0 <= face && face < T8_DTRI_FACES);
  T8_ASSERT (neigh_face != NULL);

  *neigh_face = t8_dtri_face_neighbour (t, face, n);
  /* return true if neigh is inside the root */
  return t8_dtri_is_inside_root (n);
}

void
t8_default_scheme_tri::element_set_linear_id (t8_element_t *elem, int level, t8_linearidx_t id) const
{
  T8_ASSERT (element_is_valid (elem));
  T8_ASSERT (0 <= level && level <= T8_DTRI_MAXLEVEL);
  T8_ASSERT (0 <= id && id < ((t8_linearidx_t) 1) << (2 * level));

  t8_dtri_init_linear_id ((t8_dtri_t *) elem, id, level);
}

t8_linearidx_t
t8_default_scheme_tri::element_get_linear_id (const t8_element_t *elem, int level) const
{
  T8_ASSERT (element_is_valid (elem));
  T8_ASSERT (0 <= level && level <= T8_DTRI_MAXLEVEL);

  return t8_dtri_linear_id ((t8_dtri_t *) elem, level);
}

void
t8_default_scheme_tri::element_get_first_descendant (const t8_element_t *elem, t8_element_t *desc, int level) const
{
  T8_ASSERT (element_is_valid (elem));
  T8_ASSERT (element_is_valid (desc));
  T8_ASSERT (0 <= level && level <= T8_DTRI_MAXLEVEL);
  t8_dtri_first_descendant ((t8_dtri_t *) elem, (t8_dtri_t *) desc, level);
}

void
t8_default_scheme_tri::element_get_last_descendant (const t8_element_t *elem, t8_element_t *desc, int level) const
{
  T8_ASSERT (element_is_valid (elem));
  T8_ASSERT (element_is_valid (desc));
  T8_ASSERT (0 <= level && level <= T8_DTRI_MAXLEVEL);
  t8_dtri_last_descendant ((t8_dtri_t *) elem, (t8_dtri_t *) desc, level);
}

void
t8_default_scheme_tri::element_construct_successor (const t8_element_t *elem1, t8_element_t *elem2) const
{
  T8_ASSERT (element_is_valid (elem1));
  T8_ASSERT (element_is_valid (elem2));
  T8_ASSERT (0 <= element_get_level (elem1) && element_get_level (elem1) <= T8_DTRI_MAXLEVEL);

  t8_dtri_successor ((const t8_dtri_t *) elem1, (t8_dtri_t *) elem2, element_get_level (elem1));
}

void
t8_default_scheme_tri::element_get_anchor (const t8_element_t *elem, int anchor[3]) const
{
  t8_dtri_t *tri = (t8_dtri_t *) elem;

  T8_ASSERT (element_is_valid (elem));
  anchor[0] = tri->x;
  anchor[1] = tri->y;
  anchor[2] = 0;
}

void
t8_default_scheme_tri::element_get_vertex_integer_coords (const t8_element_t *elem, int vertex, int coords[]) const
{
  T8_ASSERT (element_is_valid (elem));
  t8_dtri_compute_integer_coords ((const t8_dtri_t *) elem, vertex, coords);
}

void
t8_default_scheme_tri::element_get_vertex_reference_coords (const t8_element_t *elem, const int vertex,
                                                            double coords[]) const
{
  T8_ASSERT (element_is_valid (elem));
  t8_dtri_compute_vertex_ref_coords ((const t8_dtri_t *) elem, vertex, coords);
}

void
t8_default_scheme_tri::element_get_reference_coords (const t8_element_t *elem, const double *ref_coords,
                                                     const size_t num_coords, double *out_coords) const
{
  T8_ASSERT (element_is_valid (elem));
  t8_dtri_compute_reference_coords ((const t8_dtri_t *) elem, ref_coords, num_coords, 0, out_coords);
}

int
t8_default_scheme_tri::refines_irregular () const
{
  /*tris refine regularly */
  return 0;
}

#ifdef T8_ENABLE_DEBUG
int
t8_default_scheme_tri::element_is_valid (const t8_element_t *t) const
{
  return t8_dtri_is_valid ((const t8_dtri_t *) t);
}

void
t8_default_scheme_tri::element_to_string (const t8_element_t *elem, char *debug_string, const int string_size) const
{
  T8_ASSERT (element_is_valid (elem));
  T8_ASSERT (debug_string != NULL);
  t8_dtri_t *tri = (t8_dtri_t *) elem;
  snprintf (debug_string, string_size, "x: %i, y: %i, type: %i, level: %i", tri->x, tri->y, tri->type, tri->level);
}
#endif

void
t8_default_scheme_tri::element_new (int length, t8_element_t **elem) const
{
  /* allocate memory for a tet */
  t8_default_scheme_common::element_new (length, elem);

  /* in debug mode, set sensible default values. */
#ifdef T8_ENABLE_DEBUG
  {
    int i;
    for (i = 0; i < length; i++) {
      get_root (elem[i]);
    }
  }
#endif
}

void
t8_default_scheme_tri::element_init ([[maybe_unused]] int length, [[maybe_unused]] t8_element_t *elem) const
{
#ifdef T8_ENABLE_DEBUG
  t8_dtri_t *tris = (t8_dtri_t *) elem;
  for (int i = 0; i < length; i++) {
    t8_dtri_init (tris + i);
  }
#endif
}

void
t8_default_scheme_tri::get_root (t8_element_t *elem) const
{
  t8_dtri_t *tri = (t8_dtri_t *) elem;
  tri->level = 0;
  tri->x = 0;
  tri->y = 0;
  tri->type = 0;
}
/* use macro tri functionality */
void
t8_default_scheme_tri::element_MPI_Pack (t8_element_t **const elements, const unsigned int count, void *send_buffer,
                                         const int buffer_size, int *position, sc_MPI_Comm comm) const
{
  t8_dtri_element_pack ((t8_dtri_t **) elements, count, send_buffer, buffer_size, position, comm);
}

/* use macro tri functionality */
void
t8_default_scheme_tri::element_MPI_Pack_size (const unsigned int count, sc_MPI_Comm comm, int *pack_size) const
{
  t8_dtri_element_pack_size (count, comm, pack_size);
}

/* use macro tri functionality */
void
t8_default_scheme_tri::element_MPI_Unpack (void *recvbuf, const int buffer_size, int *position, t8_element_t **elements,
                                           const unsigned int count, sc_MPI_Comm comm) const
{
  t8_dtri_element_unpack (recvbuf, buffer_size, position, (t8_dtri_t **) elements, count, comm);
}

T8_EXTERN_C_END ();
