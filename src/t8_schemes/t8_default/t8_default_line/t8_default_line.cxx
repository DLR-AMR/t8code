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
#include <t8_schemes/t8_default/t8_default_vertex/t8_default_vertex.hxx>
#include <t8_schemes/t8_default/t8_default_line/t8_dline_bits.h>
#include <t8_schemes/t8_default/t8_default_line/t8_dline.h>
#include <t8_schemes/t8_scheme.hxx>

typedef t8_dline_t t8_default_line_t;

T8_EXTERN_C_BEGIN ();

size_t
t8_default_scheme_line::get_element_size (void) const
{
  return sizeof (t8_dline_t);
}

int
t8_default_scheme_line::get_maxlevel (void) const
{
  return T8_DLINE_MAXLEVEL;
}

int
t8_default_scheme_line::element_get_level (const t8_element_t *elem) const
{
  T8_ASSERT (element_is_valid (elem));
  return t8_dline_get_level ((const t8_dline_t *) elem);
}

void
t8_default_scheme_line::element_copy (const t8_element_t *source, t8_element_t *dest) const
{
  T8_ASSERT (element_is_valid (source));
  T8_ASSERT (element_is_valid (dest));
  t8_dline_copy ((const t8_dline_t *) source, (t8_dline_t *) dest);
}

int
t8_default_scheme_line::element_compare (const t8_element_t *elem1, const t8_element_t *elem2) const
{
  T8_ASSERT (element_is_valid (elem1));
  T8_ASSERT (element_is_valid (elem2));
  return t8_dline_compare ((const t8_dline_t *) elem1, (const t8_dline_t *) elem2);
}

int
t8_default_scheme_line::element_is_equal (const t8_element_t *elem1, const t8_element_t *elem2) const
{
  return t8_dline_equal ((const t8_dline_t *) elem1, (const t8_dline_t *) elem2);
}

void
t8_default_scheme_line::element_get_parent (const t8_element_t *elem, t8_element_t *parent) const
{
  const t8_default_line_t *l = (const t8_default_line_t *) elem;
  t8_default_line_t *p = (t8_default_line_t *) parent;

  T8_ASSERT (element_is_valid (elem));
  T8_ASSERT (element_is_valid (parent));
  t8_dline_parent (l, p);
}

void
t8_default_scheme_line::element_get_child (const t8_element_t *elem, int childid, t8_element_t *child) const
{
  const t8_default_line_t *l = (const t8_default_line_t *) elem;
  t8_default_line_t *c = (t8_default_line_t *) child;

  T8_ASSERT (element_is_valid (elem));
  T8_ASSERT (element_is_valid (child));
  t8_dline_child (l, childid, c);
}

void
t8_default_scheme_line::element_get_nca (const t8_element_t *elem1, const t8_element_t *elem2, t8_element_t *nca) const
{
  T8_ASSERT (element_is_valid (elem1));
  T8_ASSERT (element_is_valid (elem2));
  T8_ASSERT (element_is_valid (nca));
  t8_dline_nearest_common_ancestor ((const t8_dline_t *) elem1, (const t8_dline_t *) elem2, (t8_dline_t *) nca);
}

t8_element_shape_t
t8_default_scheme_line::element_get_face_shape ([[maybe_unused]] const t8_element_t *elem,
                                                [[maybe_unused]] int face) const
{
  T8_ASSERT (element_is_valid (elem));
  return T8_ECLASS_VERTEX;
}

void
t8_default_scheme_line::element_get_children_at_face (const t8_element_t *elem, int face, t8_element_t *children[],
                                                      [[maybe_unused]] int num_children, int *child_indices) const
{
  T8_ASSERT (element_is_valid (elem));
  T8_ASSERT (0 <= face && face < T8_DLINE_FACES);
  T8_ASSERT (num_children == 1);
  T8_ASSERT (element_is_valid (children[0]));

  /* We have exactly one child at a face and this is child 0 if face = 0
   * and child 1 if face = 1 */
  if (child_indices != NULL) {
    *child_indices = face;
  }
  t8_dline_child ((const t8_dline_t *) elem, face, (t8_dline_t *) children[0]);
}

int
t8_default_scheme_line::element_face_get_child_face ([[maybe_unused]] const t8_element_t *elem, int face,
                                                     [[maybe_unused]] int face_child) const
{
  T8_ASSERT (element_is_valid (elem));
  T8_ASSERT (0 <= face && face < T8_DLINE_FACES);
  T8_ASSERT (face_child == 0);

  /* The face id of the child is the same as the face */
  return face;
}

int
t8_default_scheme_line::element_face_get_parent_face (const t8_element_t *elem, int face) const
{
  /* The number of faces does not change from parent to child */
  T8_ASSERT (element_is_valid (elem));
  T8_ASSERT (0 <= face && face < T8_DLINE_FACES);
  return t8_dline_face_parent_face ((const t8_dline_t *) elem, face);
}

int
t8_default_scheme_line::element_get_tree_face ([[maybe_unused]] const t8_element_t *elem, int face) const
{
  /* The number of faces does not change from tree to element */
  T8_ASSERT (element_is_valid (elem));
  T8_ASSERT (0 <= face && face < T8_DLINE_FACES);
  return face;
}

void
t8_default_scheme_line::element_transform_face (const t8_element_t *elem1, t8_element_t *elem2, int orientation,
                                                [[maybe_unused]] int sign, [[maybe_unused]] int is_smaller_face) const
{
  T8_ASSERT (element_is_valid (elem1));
  T8_ASSERT (element_is_valid (elem2));
  T8_ASSERT (orientation == 0 || orientation == 1);

  /* We can ignore is_smaller_face, since for lines the orientation is independent
   * of the face. */
  t8_dline_transform_face ((const t8_dline_t *) elem1, (t8_dline_t *) elem2, orientation);
}

/** Given a boundary face inside a root tree's face construct
 *  the element inside the root tree that has the given face as a
 *  face. */
int
t8_default_scheme_line::element_extrude_face (const t8_element_t *face, t8_element_t *elem, int root_face,
                                              [[maybe_unused]] const t8_scheme *scheme) const
{
  T8_ASSERT (element_is_valid (elem));
  T8_ASSERT (scheme->element_is_valid (T8_ECLASS_VERTEX, face));

  return t8_dline_extrude_face ((const t8_dvertex_t *) face, root_face, (t8_dline_t *) elem);
}

/** Construct the boundary element at a specific face. */
void
t8_default_scheme_line::element_get_boundary_face (const t8_element_t *elem, [[maybe_unused]] int face,
                                                   t8_element_t *boundary,
                                                   [[maybe_unused]] const t8_scheme *scheme) const
{
  T8_ASSERT (element_is_valid (elem));
  T8_ASSERT (scheme->element_is_valid (T8_ECLASS_VERTEX, boundary));
  T8_ASSERT (0 <= face && face < T8_DLINE_FACES);

  /* Since each vertex is the same, we just construct a vertex of the same level
   * as elem. */
  t8_default_scheme_vertex::element_set_linear_id (boundary, element_get_level (elem), 0);
}

/** Construct the first descendant of an element that touches a given face.   */
void
t8_default_scheme_line::element_get_first_descendant_face (const t8_element_t *elem, int face, t8_element_t *first_desc,
                                                           int level) const
{
  T8_ASSERT (element_is_valid (elem));
  T8_ASSERT (element_is_valid (first_desc));
  T8_ASSERT (0 <= face && face < T8_DLINE_FACES);

  T8_ASSERT (0 <= level && level <= T8_DLINE_MAXLEVEL);
  /* The first descendant at the face is the first desc of elem if face = 0.
   * If face = 1 it is the last desc of elem. */
  if (face == 0) {
    t8_dline_first_descendant ((const t8_dline_t *) elem, (t8_dline_t *) first_desc, level);
  }
  else {
    T8_ASSERT (face == 1);
    t8_dline_last_descendant ((const t8_dline_t *) elem, (t8_dline_t *) first_desc, level);
  }
}

void
t8_default_scheme_line::element_get_last_descendant_face (const t8_element_t *elem, int face, t8_element_t *last_desc,
                                                          int level) const
{
  T8_ASSERT (element_is_valid (elem));
  T8_ASSERT (element_is_valid (last_desc));
  T8_ASSERT (0 <= face && face < T8_DLINE_FACES);
  T8_ASSERT (0 <= level && level <= T8_DLINE_MAXLEVEL);

  /* The last descendant is the same as the first descendant. */
  element_get_first_descendant_face (elem, face, last_desc, level);
}

int
t8_default_scheme_line::element_is_root_boundary (const t8_element_t *elem, int face) const
{
  T8_ASSERT (element_is_valid (elem));
  T8_ASSERT (0 <= face && face < T8_DLINE_FACES);

  return t8_dline_is_root_boundary ((const t8_dline_t *) elem, face);
}

int
t8_default_scheme_line::element_get_face_neighbor_inside (const t8_element_t *elem, t8_element_t *neigh, int face,
                                                          int *neigh_face) const
{
  T8_ASSERT (element_is_valid (elem));
  T8_ASSERT (element_is_valid (neigh));
  T8_ASSERT (0 <= face && face < T8_DLINE_FACES);

  t8_dline_face_neighbour ((const t8_dline_t *) elem, (t8_dline_t *) neigh, face, neigh_face);
  return t8_dline_is_inside_root ((t8_dline_t *) neigh);
}

void
t8_default_scheme_line::element_set_linear_id (t8_element_t *elem, int level, t8_linearidx_t id) const
{
  T8_ASSERT (element_is_valid (elem));
  T8_ASSERT (0 <= level && level <= T8_DLINE_MAXLEVEL);
  T8_ASSERT (id < ((t8_linearidx_t) 1) << level);

  t8_dline_init_linear_id ((t8_default_line_t *) elem, level, id);
}

void
t8_default_scheme_line::element_construct_successor (const t8_element_t *elem1, t8_element_t *elem2) const
{
  T8_ASSERT (element_is_valid (elem1));
  T8_ASSERT (element_is_valid (elem2));
  T8_ASSERT (1 <= element_get_level (elem1) && element_get_level (elem1) <= T8_DLINE_MAXLEVEL);

  t8_dline_successor ((const t8_default_line_t *) elem1, (t8_default_line_t *) elem2, element_get_level (elem1));
}

void
t8_default_scheme_line::element_get_first_descendant (const t8_element_t *elem, t8_element_t *desc, int level) const
{
  T8_ASSERT (element_is_valid (elem));
  T8_ASSERT (element_is_valid (desc));

  T8_ASSERT (0 <= level && level <= T8_DLINE_MAXLEVEL);
  t8_dline_first_descendant ((const t8_dline_t *) elem, (t8_dline_t *) desc, level);
}

void
t8_default_scheme_line::element_get_last_descendant (const t8_element_t *elem, t8_element_t *desc, int level) const
{
  T8_ASSERT (element_is_valid (elem));
  T8_ASSERT (element_is_valid (desc));
  T8_ASSERT (0 <= level && level <= T8_DLINE_MAXLEVEL);
  t8_dline_last_descendant ((const t8_dline_t *) elem, (t8_dline_t *) desc, level);
}

t8_linearidx_t
t8_default_scheme_line::element_get_linear_id (const t8_element_t *elem, int level) const
{
  T8_ASSERT (element_is_valid (elem));
  T8_ASSERT (0 <= level && level <= T8_DLINE_MAXLEVEL);

  return t8_dline_linear_id ((const t8_dline_t *) elem, level);
}

int
t8_default_scheme_line::element_get_num_faces ([[maybe_unused]] const t8_element_t *elem) const
{
  T8_ASSERT (element_is_valid (elem));
  return T8_DLINE_FACES;
}

int
t8_default_scheme_line::element_get_max_num_faces ([[maybe_unused]] const t8_element_t *elem) const
{
  T8_ASSERT (element_is_valid (elem));
  return T8_DLINE_FACES;
}

int
t8_default_scheme_line::element_get_num_children ([[maybe_unused]] const t8_element_t *elem) const
{
  T8_ASSERT (element_is_valid (elem));
  return T8_DLINE_CHILDREN;
}

int
t8_default_scheme_line::element_get_num_face_children ([[maybe_unused]] const t8_element_t *elem,
                                                       [[maybe_unused]] int face) const
{
  T8_ASSERT (element_is_valid (elem));
  T8_ASSERT (0 <= face && face < T8_DLINE_FACES);

  return T8_DLINE_FACE_CHILDREN;
}

int
t8_default_scheme_line::element_get_child_id (const t8_element_t *elem) const
{
  T8_ASSERT (element_is_valid (elem));
  return t8_dline_child_id ((const t8_dline_t *) elem);
}

void
t8_default_scheme_line::element_get_children (const t8_element_t *elem, [[maybe_unused]] int length,
                                              t8_element_t *c[]) const
{
  T8_ASSERT (element_is_valid (elem));
  T8_ASSERT (length == T8_DLINE_CHILDREN);

  t8_dline_childrenpv ((const t8_dline_t *) elem, (t8_dline_t **) c);
}

int
t8_default_scheme_line::element_get_ancestor_id (const t8_element_t *elem, int level) const
{
  T8_ASSERT (element_is_valid (elem));
  return t8_dline_ancestor_id ((const t8_dline_t *) elem, level);
}

int
t8_default_scheme_line::elements_are_family (t8_element_t *const *fam) const
{
#if T8_ENABLE_DEBUG
  int i;
  for (i = 0; i < T8_DLINE_CHILDREN; i++) {
    T8_ASSERT (element_is_valid (fam[i]));
  }
#endif
  return t8_dline_is_familypv ((const t8_dline_t **) fam);
}

int
t8_default_scheme_line::refines_irregular () const
{
  /*lines always refine regularly */
  return 0;
}

#if T8_ENABLE_DEBUG

void
t8_default_scheme_line::element_to_string (const t8_element_t *elem, char *debug_string, const int string_size) const
{
  T8_ASSERT (element_is_valid (elem));
  T8_ASSERT (debug_string != NULL);
  t8_dline_t *line = (t8_dline_t *) elem;
  snprintf (debug_string, string_size, "x: %i, level: %i", line->x, line->level);
}
#endif

void
t8_default_scheme_line::element_new (int length, t8_element_t **elem) const
{
  /* allocate memory for a line */
  t8_default_scheme_common::element_new (length, elem);

  /* in debug mode, set sensible default values. */
#if T8_ENABLE_DEBUG
  {
    for (int i = 0; i < length; i++) {
      set_to_root (elem[i]);
    }
  }
#endif
}

void
t8_default_scheme_line::element_init ([[maybe_unused]] int length, [[maybe_unused]] t8_element_t *elem) const
{
#if T8_ENABLE_DEBUG
  t8_dline_t *lines = (t8_dline_t *) elem;
  for (int i = 0; i < length; i++) {
    t8_dline_init (lines + i);
  }
#endif
}

/* each line is packed as an x coordinate and the level */
void
t8_default_scheme_line::element_MPI_Pack (t8_element_t **const elements, const unsigned int count, void *send_buffer,
                                          const int buffer_size, int *position, sc_MPI_Comm comm) const
{
  t8_default_line_t **lines = (t8_default_line_t **) elements;
  int mpiret;
  for (unsigned int ielem = 0; ielem < count; ielem++) {
    mpiret = sc_MPI_Pack (&(lines[ielem]->x), 1, sc_MPI_INT, send_buffer, buffer_size, position, comm);
    SC_CHECK_MPI (mpiret);
    mpiret = sc_MPI_Pack (&lines[ielem]->level, 1, sc_MPI_INT8_T, send_buffer, buffer_size, position, comm);
    SC_CHECK_MPI (mpiret);
  }
}

/* each line is packed as an x coordinate and the level */
void
t8_default_scheme_line::element_MPI_Pack_size (const unsigned int count, sc_MPI_Comm comm, int *pack_size) const
{
  int singlesize = 0;
  int datasize = 0;
  int mpiret;

  mpiret = sc_MPI_Pack_size (1, sc_MPI_INT, comm, &datasize);
  SC_CHECK_MPI (mpiret);
  singlesize += datasize;

  mpiret = sc_MPI_Pack_size (1, sc_MPI_INT8_T, comm, &datasize);
  SC_CHECK_MPI (mpiret);
  singlesize += datasize;

  *pack_size = count * singlesize;
}

/* each line is packed as an x coordinate and the level */
void
t8_default_scheme_line::element_MPI_Unpack (void *recvbuf, const int buffer_size, int *position,
                                            t8_element_t **elements, const unsigned int count, sc_MPI_Comm comm) const
{
  int mpiret;
  t8_default_line_t **lines = (t8_default_line_t **) elements;
  for (unsigned int ielem = 0; ielem < count; ielem++) {
    mpiret = sc_MPI_Unpack (recvbuf, buffer_size, position, &(lines[ielem]->x), 1, sc_MPI_INT, comm);
    SC_CHECK_MPI (mpiret);
    mpiret = sc_MPI_Unpack (recvbuf, buffer_size, position, &(lines[ielem]->level), 1, sc_MPI_INT8_T, comm);
    SC_CHECK_MPI (mpiret);
  }
}

void
t8_default_scheme_line::set_to_root (t8_element_t *elem) const
{
  t8_dline_t *line = (t8_dline_t *) elem;
  line->level = 0;
  line->x = 0;
}
T8_EXTERN_C_END ();
