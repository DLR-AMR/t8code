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

#include "t8_default_common_cxx.hxx"
#include "t8_default_vertex_cxx.hxx"
#include "t8_dvertex_bits.h"
#include "t8_dline_bits.h"
#include "t8_dline.h"

typedef t8_dline_t  t8_default_line_t;

T8_EXTERN_C_BEGIN ();

int
t8_default_scheme_line_c::t8_element_maxlevel (void)
{
  return T8_DLINE_MAXLEVEL;
}

int
t8_default_scheme_line_c::t8_element_level (const t8_element_t * elem)
{
  T8_ASSERT (t8_element_is_valid (elem));
  return t8_dline_get_level ((const t8_dline_t *) elem);
}

void
t8_default_scheme_line_c::t8_element_copy (const t8_element_t * source,
                                           t8_element_t * dest)
{
  T8_ASSERT (t8_element_is_valid (source));
  T8_ASSERT (t8_element_is_valid (dest));
  t8_dline_copy ((const t8_dline_t *) source, (t8_dline_t *) dest);
}

int
t8_default_scheme_line_c::t8_element_compare (const t8_element_t * elem1,
                                              const t8_element_t * elem2)
{
  T8_ASSERT (t8_element_is_valid (elem1));
  T8_ASSERT (t8_element_is_valid (elem2));
  return t8_dline_compare ((const t8_dline_t *) elem1,
                           (const t8_dline_t *) elem2);
}

void
t8_default_scheme_line_c::t8_element_parent (const t8_element_t * elem,
                                             t8_element_t * parent)
{
  const t8_default_line_t *l = (const t8_default_line_t *) elem;
  t8_default_line_t  *p = (t8_default_line_t *) parent;

  T8_ASSERT (t8_element_is_valid (elem));
  T8_ASSERT (t8_element_is_valid (parent));
  t8_dline_parent (l, p);
}

void
t8_default_scheme_line_c::t8_element_child (const t8_element_t * elem,
                                            int childid, t8_element_t * child)
{
  const t8_default_line_t *l = (const t8_default_line_t *) elem;
  t8_default_line_t  *c = (t8_default_line_t *) child;

  T8_ASSERT (t8_element_is_valid (elem));
  T8_ASSERT (t8_element_is_valid (child));
  t8_dline_child (l, childid, c);
}

void
t8_default_scheme_line_c::t8_element_nca (const t8_element_t * elem1,
                                          const t8_element_t * elem2,
                                          t8_element_t * nca)
{
  T8_ASSERT (t8_element_is_valid (elem1));
  T8_ASSERT (t8_element_is_valid (elem2));
  T8_ASSERT (t8_element_is_valid (nca));
  t8_dline_nearest_common_ancestor ((const t8_dline_t *) elem1,
                                    (const t8_dline_t *) elem2,
                                    (t8_dline_t *) nca);
}

t8_element_shape_t
  t8_default_scheme_line_c::t8_element_face_shape (const t8_element_t * elem,
                                                   int face)
{
  T8_ASSERT (t8_element_is_valid (elem));
  return T8_ECLASS_VERTEX;
}

void
t8_default_scheme_line_c::t8_element_children_at_face (const t8_element_t *
                                                       elem, int face,
                                                       t8_element_t *
                                                       children[],
                                                       int num_children,
                                                       int *child_indices)
{
  T8_ASSERT (t8_element_is_valid (elem));
  T8_ASSERT (0 <= face && face < T8_DLINE_FACES);
  T8_ASSERT (num_children == 1);
  T8_ASSERT (t8_element_is_valid (children[0]));

  /* We have exactly one child at a face and this is child 0 if face = 0
   * and child 1 if face = 1 */
  if (child_indices != NULL) {
    *child_indices = face;
  }
  t8_dline_child ((const t8_dline_t *) elem, face,
                  (t8_dline_t *) children[0]);
}

int
t8_default_scheme_line_c::t8_element_face_child_face (const t8_element_t *
                                                      elem, int face,
                                                      int face_child)
{
  T8_ASSERT (t8_element_is_valid (elem));
  T8_ASSERT (0 <= face && face < T8_DLINE_FACES);
  T8_ASSERT (face_child == 0);

  /* The face id of the child is the same as the face */
  return face;
}

int
t8_default_scheme_line_c::t8_element_face_parent_face (const t8_element_t *
                                                       elem, int face)
{
  /* The number of faces does not change from parent to child */
  T8_ASSERT (t8_element_is_valid (elem));
  T8_ASSERT (0 <= face && face < T8_DLINE_FACES);
  return t8_dline_face_parent_face ((const t8_dline_t *) elem, face);
}

int
t8_default_scheme_line_c::t8_element_tree_face (const t8_element_t * elem,
                                                int face)
{
  /* The number of faces does not change from tree to element */
  T8_ASSERT (t8_element_is_valid (elem));
  T8_ASSERT (0 <= face && face < T8_DLINE_FACES);
  return face;
}

void
t8_default_scheme_line_c::t8_element_transform_face (const t8_element_t *
                                                     elem1,
                                                     t8_element_t * elem2,
                                                     int orientation,
                                                     int sign,
                                                     int is_smaller_face)
{
  T8_ASSERT (t8_element_is_valid (elem1));
  T8_ASSERT (t8_element_is_valid (elem2));
  T8_ASSERT (orientation == 0 || orientation == 1);

  /* We can ignore is_smaller_face, since for lines the orientation is independent
   * of the face. */
  t8_dline_transform_face ((const t8_dline_t *) elem1, (t8_dline_t *) elem2,
                           orientation);
}

/** Given a boundary face inside a root tree's face construct
 *  the element inside the root tree that has the given face as a
 *  face. */
int
t8_default_scheme_line_c::t8_element_extrude_face (const t8_element_t * face,
                                                   const t8_eclass_scheme_c *
                                                   face_scheme,
                                                   t8_element_t * elem,
                                                   int root_face)
{
  T8_ASSERT (t8_element_is_valid (elem));
  T8_ASSERT (T8_COMMON_IS_TYPE
             (face_scheme, const t8_default_scheme_vertex_c *));
  T8_ASSERT (face_scheme->t8_element_is_valid (face));

  return t8_dline_extrude_face ((const t8_dvertex_t *) face, root_face,
                                (t8_dline_t *) elem);
}

/** Construct the boundary element at a specific face. */
void
t8_default_scheme_line_c::t8_element_boundary_face (const t8_element_t * elem,
                                                    int face,
                                                    t8_element_t * boundary,
                                                    const t8_eclass_scheme_c *
                                                    boundary_scheme)
{
  T8_ASSERT (t8_element_is_valid (elem));
  T8_ASSERT (T8_COMMON_IS_TYPE
             (boundary_scheme, const t8_default_scheme_vertex_c *));
  T8_ASSERT (boundary_scheme->eclass == T8_ECLASS_VERTEX);
  T8_ASSERT (boundary_scheme->t8_element_is_valid (boundary));
  T8_ASSERT (0 <= face && face < T8_DLINE_FACES);

  /* Since each vertex is the same, we just construc a vertex of the same level
   * as elem. */
  t8_dvertex_init_linear_id ((t8_dvertex_t *) boundary,
                             t8_element_level (elem), 0);
}

/** Construct the first descendant of an element that touches a given face.   */
void
t8_default_scheme_line_c::t8_element_first_descendant_face (const t8_element_t
                                                            * elem, int face,
                                                            t8_element_t *
                                                            first_desc,
                                                            int level)
{
  T8_ASSERT (t8_element_is_valid (elem));
  T8_ASSERT (t8_element_is_valid (first_desc));
  T8_ASSERT (0 <= face && face < T8_DLINE_FACES);

  T8_ASSERT (0 <= level && level <= T8_DLINE_MAXLEVEL);
  /* The first descandant at the face is the first desc of elem if face = 0.
   * If face = 1 it is the last desc of elem. */
  if (face == 0) {
    t8_dline_first_descendant ((const t8_dline_t *) elem,
                               (t8_dline_t *) first_desc, level);
  }
  else {
    T8_ASSERT (face == 1);
    t8_dline_last_descendant ((const t8_dline_t *) elem,
                              (t8_dline_t *) first_desc, level);

  }
}

void
t8_default_scheme_line_c::t8_element_last_descendant_face (const t8_element_t
                                                           * elem, int face,
                                                           t8_element_t *
                                                           last_desc,
                                                           int level)
{
  T8_ASSERT (t8_element_is_valid (elem));
  T8_ASSERT (t8_element_is_valid (last_desc));
  T8_ASSERT (0 <= face && face < T8_DLINE_FACES);
  T8_ASSERT (0 <= level && level <= T8_DLINE_MAXLEVEL);

  /* The last descendant is the same as the first descendant. */
  t8_element_first_descendant_face (elem, face, last_desc, level);
}

int
t8_default_scheme_line_c::t8_element_is_root_boundary (const t8_element_t *
                                                       elem, int face)
{
  T8_ASSERT (t8_element_is_valid (elem));
  T8_ASSERT (0 <= face && face < T8_DLINE_FACES);

  return t8_dline_is_root_boundary ((const t8_dline_t *) elem, face);
}

int
t8_default_scheme_line_c::t8_element_face_neighbor_inside (const t8_element_t
                                                           * elem,
                                                           t8_element_t *
                                                           neigh, int face,
                                                           int *neigh_face)
{
  T8_ASSERT (t8_element_is_valid (elem));
  T8_ASSERT (t8_element_is_valid (neigh));
  T8_ASSERT (0 <= face && face < T8_DLINE_FACES);

  t8_dline_face_neighbour ((const t8_dline_t *) elem, (t8_dline_t *) neigh,
                           face, neigh_face);
  return t8_dline_is_inside_root ((t8_dline_t *) neigh);
}

void
t8_default_scheme_line_c::t8_element_set_linear_id (t8_element_t * elem,
                                                    int level,
                                                    t8_linearidx_t id)
{
  T8_ASSERT (t8_element_is_valid (elem));
  T8_ASSERT (0 <= level && level <= T8_DLINE_MAXLEVEL);
  T8_ASSERT (0 <= id && id < ((t8_linearidx_t) 1) << level);

  t8_dline_init_linear_id ((t8_default_line_t *) elem, level, id);
}

void
t8_default_scheme_line_c::t8_element_successor (const t8_element_t * elem1,
                                                t8_element_t * elem2,
                                                int level)
{
  T8_ASSERT (t8_element_is_valid (elem1));
  T8_ASSERT (t8_element_is_valid (elem2));
  T8_ASSERT (1 <= level && level <= T8_DLINE_MAXLEVEL);

  t8_dline_successor ((const t8_default_line_t *) elem1,
                      (t8_default_line_t *) elem2, level);
}

void
t8_default_scheme_line_c::t8_element_first_descendant (const t8_element_t *
                                                       elem,
                                                       t8_element_t * desc,
                                                       int level)
{
  T8_ASSERT (t8_element_is_valid (elem));
  T8_ASSERT (t8_element_is_valid (desc));

  T8_ASSERT (0 <= level && level <= T8_DLINE_MAXLEVEL);
  t8_dline_first_descendant ((const t8_dline_t *) elem, (t8_dline_t *) desc,
                             level);
}

void
t8_default_scheme_line_c::t8_element_last_descendant (const t8_element_t *
                                                      elem,
                                                      t8_element_t * desc,
                                                      int level)
{
  T8_ASSERT (t8_element_is_valid (elem));
  T8_ASSERT (t8_element_is_valid (desc));
  T8_ASSERT (0 <= level && level <= T8_DLINE_MAXLEVEL);
  t8_dline_last_descendant ((const t8_dline_t *) elem, (t8_dline_t *) desc,
                            level);
}

void
t8_default_scheme_line_c::t8_element_vertex_coords (const t8_element_t * t,
                                                    int vertex, int coords[])
{
  T8_ASSERT (t8_element_is_valid (t));
  t8_dline_vertex_coords ((const t8_dline_t *) t, vertex, coords);
}

void
t8_default_scheme_line_c::t8_element_vertex_reference_coords (const
                                                              t8_element_t *
                                                              t, int vertex,
                                                              double coords[])
{
  T8_ASSERT (t8_element_is_valid (t));
  t8_dline_vertex_ref_coords ((const t8_dline_t *) t, vertex, coords);
}

int
t8_default_scheme_line_c::t8_element_root_len (const t8_element_t * elem)
{
  T8_ASSERT (t8_element_is_valid (elem));
  return T8_DLINE_ROOT_LEN;
}

t8_linearidx_t
  t8_default_scheme_line_c::t8_element_get_linear_id (const t8_element_t *
                                                      elem, int level)
{
  T8_ASSERT (t8_element_is_valid (elem));
  T8_ASSERT (0 <= level && level <= T8_DLINE_MAXLEVEL);

  return t8_dline_linear_id ((const t8_dline_t *) elem, level);
}

int
t8_default_scheme_line_c::t8_element_num_faces (const t8_element_t * elem)
{
  T8_ASSERT (t8_element_is_valid (elem));
  return T8_DLINE_FACES;
}

int
t8_default_scheme_line_c::t8_element_max_num_faces (const t8_element_t * elem)
{
  T8_ASSERT (t8_element_is_valid (elem));
  return T8_DLINE_FACES;
}

int
t8_default_scheme_line_c::t8_element_num_children (const t8_element_t * elem)
{
  T8_ASSERT (t8_element_is_valid (elem));
  return T8_DLINE_CHILDREN;
}

int
t8_default_scheme_line_c::t8_element_num_face_children (const t8_element_t *
                                                        elem, int face)
{
  T8_ASSERT (t8_element_is_valid (elem));
  T8_ASSERT (0 <= face && face < T8_DLINE_FACES);

  return T8_DLINE_FACE_CHILDREN;
}

int
t8_default_scheme_line_c::t8_element_child_id (const t8_element_t * elem)
{
  T8_ASSERT (t8_element_is_valid (elem));
  return t8_dline_child_id ((const t8_dline_t *) elem);
}

void
t8_default_scheme_line_c::t8_element_children (const t8_element_t * elem,
                                               int length, t8_element_t * c[])
{
  T8_ASSERT (t8_element_is_valid (elem));
  T8_ASSERT (length == T8_DLINE_CHILDREN);

  t8_dline_childrenpv ((const t8_dline_t *) elem, (t8_dline_t **) c);
}

int
t8_default_scheme_line_c::t8_element_ancestor_id (const t8_element_t * elem,
                                                  int level)
{
  T8_ASSERT (t8_element_is_valid (elem));
  return t8_dline_ancestor_id ((const t8_dline_t *) elem, level);
}

int
t8_default_scheme_line_c::t8_element_is_family (t8_element_t ** fam)
{
#ifdef T8_ENABLE_DEBUG
  int                 i;
  for (i = 0; i < T8_DLINE_CHILDREN; i++) {
    T8_ASSERT (t8_element_is_valid (fam[i]));
  }
#endif
  return t8_dline_is_familypv ((const t8_dline_t **) fam);
}

#ifdef T8_ENABLE_DEBUG
/* *INDENT-OFF* */
/* indent bug, indent adds a second "const" modifier */
int
t8_default_scheme_line_c::t8_element_is_valid (const t8_element_t * elem) const
/* *INDENT-ON* */
{
  return t8_dline_is_valid ((const t8_dline_t *) elem);
}
#endif

void
t8_default_scheme_line_c::t8_element_new (int length, t8_element_t ** elem)
{
  /* allocate memory for a line */
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
t8_default_scheme_line_c::t8_element_init (int length, t8_element_t * elem,
                                           int new_called)
{
#ifdef T8_ENABLE_DEBUG
  if (!new_called) {
    int                 i;
    t8_dline_t         *lines = (t8_dline_t *) elem;
    for (i = 0; i < length; i++) {
      t8_dline_init (lines + i);
    }
  }
#endif
}

/* Constructor */
t8_default_scheme_line_c::t8_default_scheme_line_c (void)
{
  eclass = T8_ECLASS_LINE;
  element_size = sizeof (t8_default_line_t);
  ts_context = sc_mempool_new (element_size);
}

t8_default_scheme_line_c::~t8_default_scheme_line_c ()
{
  /* This destructor is empty since the destructor of the
   * default_common scheme is called automatically and it
   * suffices to destroy the quad_scheme.
   * However we need to provide an implementation of the destructor
   * and hence this empty function. */
}

T8_EXTERN_C_END ();
