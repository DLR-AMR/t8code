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
#include "t8_default_tri_cxx.hxx"
#include "t8_dtri_bits.h"
#include "t8_dline_bits.h"
#include "t8_dtet.h"
#include "t8_dtri_connectivity.h"
#include "t8_eclass.h"

/* We want to export the whole implementation to be callable from "C" */
T8_EXTERN_C_BEGIN ();

int
t8_default_scheme_tri_c::t8_element_maxlevel (void)
{
  return T8_DTRI_MAXLEVEL;
}

int
t8_default_scheme_tri_c::t8_element_level (const t8_element_t * elem)
{
  T8_ASSERT (t8_element_is_valid (elem));
  return t8_dtri_get_level ((t8_dtri_t *) elem);
}

void
t8_default_scheme_tri_c::t8_element_copy (const t8_element_t * source,
                                          t8_element_t * dest)
{
  T8_ASSERT (t8_element_is_valid (source));
  T8_ASSERT (t8_element_is_valid (dest));
  t8_dtri_copy ((const t8_dtri_t *) source, (t8_dtri_t *) dest);
}

int
t8_default_scheme_tri_c::t8_element_compare (const t8_element_t * elem1,
                                             const t8_element_t * elem2)
{
  T8_ASSERT (t8_element_is_valid (elem1));
  T8_ASSERT (t8_element_is_valid (elem2));

  return t8_dtri_compare ((const t8_dtri_t *) elem1,
                          (const t8_dtri_t *) elem2);
}

void
t8_default_scheme_tri_c::t8_element_parent (const t8_element_t * elem,
                                            t8_element_t * parent)
{
  const t8_dtri_t    *t = (const t8_dtri_t *) elem;
  t8_dtri_t          *p = (t8_dtri_t *) parent;

  T8_ASSERT (t8_element_is_valid (elem));
  T8_ASSERT (t8_element_is_valid (parent));
  t8_dtri_parent (t, p);
}

void
t8_default_scheme_tri_c::t8_element_sibling (const t8_element_t * elem,
                                             int sibid,
                                             t8_element_t * sibling)
{
  const t8_dtri_t    *t = (const t8_dtri_t *) elem;
  t8_dtri_t          *s = (t8_dtri_t *) sibling;

  T8_ASSERT (t8_element_is_valid (elem));
  t8_dtri_sibling (t, sibid, s);
}

int
t8_default_scheme_tri_c::t8_element_num_faces (const t8_element_t * elem)
{
  T8_ASSERT (t8_element_is_valid (elem));
  return T8_DTRI_FACES;
}

int
t8_default_scheme_tri_c::t8_element_max_num_faces (const t8_element_t * elem)
{
  return T8_DTRI_FACES;
}

int
t8_default_scheme_tri_c::t8_element_num_children (const t8_element_t * elem)
{
  T8_ASSERT (t8_element_is_valid (elem));
  return T8_DTRI_CHILDREN;
}

int
t8_default_scheme_tri_c::t8_element_num_face_children (const t8_element_t *
                                                       elem, int face)
{
  T8_ASSERT (t8_element_is_valid (elem));
  return T8_DTRI_FACE_CHILDREN;
}

int
t8_default_scheme_tri_c::t8_element_get_face_corner (const t8_element_t *
                                                     element, int face,
                                                     int corner)
{
  T8_ASSERT (t8_element_is_valid (element));
  T8_ASSERT (0 <= face && face < T8_DTRI_FACES);
  T8_ASSERT (0 <= corner && corner < 2);
  return t8_dtri_face_corner[face][corner];
}

int
t8_default_scheme_tri_c::t8_element_get_corner_face (const t8_element_t *
                                                     element, int corner,
                                                     int face)
{
  T8_ASSERT (t8_element_is_valid (element));
  T8_ASSERT (0 <= corner && corner < T8_DTRI_CORNERS);
  T8_ASSERT (0 <= face && face < 2);
  return t8_dtri_corner_face[corner][face];
}

void
t8_default_scheme_tri_c::t8_element_child (const t8_element_t * elem,
                                           int childid, t8_element_t * child)
{
  const t8_dtri_t    *t = (const t8_dtri_t *) elem;
  t8_dtri_t          *c = (t8_dtri_t *) child;

  T8_ASSERT (t8_element_is_valid (elem));
  T8_ASSERT (t8_element_is_valid (child));
  t8_dtri_child (t, childid, c);
}

void
t8_default_scheme_tri_c::t8_element_children (const t8_element_t * elem,
                                              int length, t8_element_t * c[])
{
  T8_ASSERT (t8_element_is_valid (elem));
#ifdef T8_ENABLE_DEBUG
  {
    int                 j;
    for (j = 0; j < length; j++) {
      T8_ASSERT (t8_element_is_valid (c[j]));
    }
  }
#endif
  T8_ASSERT (length == T8_DTRI_CHILDREN);

  t8_dtri_childrenpv ((const t8_dtri_t *) elem, (t8_dtri_t **) c);
}

int
t8_default_scheme_tri_c::t8_element_child_id (const t8_element_t * elem)
{
  T8_ASSERT (t8_element_is_valid (elem));
  return t8_dtri_child_id ((t8_dtri_t *) elem);
}

int
t8_default_scheme_tri_c::t8_element_ancestor_id (const t8_element_t * elem,
                                                 int level)
{
  return t8_dtri_ancestor_id ((t8_dtri_t *) elem, level);
}

int
t8_default_scheme_tri_c::t8_element_is_family (t8_element_t ** fam)
{
#ifdef T8_ENABLE_DEBUG
  {
    int                 j;
    for (j = 0; j < T8_DTRI_CHILDREN; j++) {
      T8_ASSERT (t8_element_is_valid (fam[j]));
    }
  }
#endif
  return t8_dtri_is_familypv ((const t8_dtri_t **) fam);
}

void
t8_default_scheme_tri_c::t8_element_nca (const t8_element_t * elem1,
                                         const t8_element_t * elem2,
                                         t8_element_t * nca)
{
  const t8_dtri_t    *t1 = (const t8_dtri_t *) elem1;
  const t8_dtri_t    *t2 = (const t8_dtri_t *) elem2;
  t8_dtri_t          *c = (t8_dtri_t *) nca;

  T8_ASSERT (t8_element_is_valid (elem1));
  T8_ASSERT (t8_element_is_valid (elem2));
  t8_dtri_nearest_common_ancestor (t1, t2, c);
}

t8_element_shape_t
  t8_default_scheme_tri_c::t8_element_face_shape (const t8_element_t * elem,
                                                  int face)
{
  T8_ASSERT (t8_element_is_valid (elem));
  return T8_ECLASS_LINE;
}

void
t8_default_scheme_tri_c::t8_element_children_at_face (const t8_element_t *
                                                      elem, int face,
                                                      t8_element_t *
                                                      children[],
                                                      int num_children,
                                                      int *child_indices)
{
  const t8_dtri_t    *t = (const t8_dtri_t *) elem;
  t8_dtri_t         **c = (t8_dtri_t **) children;

  T8_ASSERT (t8_element_is_valid (elem));
#ifdef T8_ENABLE_DEBUG
  {
    int                 j;
    for (j = 0; j < num_children; j++) {
      T8_ASSERT (t8_element_is_valid (children[j]));
    }
  }
#endif
  T8_ASSERT (0 <= face && face < T8_DTRI_FACES);
  T8_ASSERT (num_children == T8_DTRI_FACE_CHILDREN);

  t8_dtri_children_at_face (t, face, c, num_children, child_indices);
}

int
t8_default_scheme_tri_c::t8_element_face_child_face (const t8_element_t *
                                                     elem, int face,
                                                     int face_child)
{
  T8_ASSERT (t8_element_is_valid (elem));
  T8_ASSERT (0 <= face && face < T8_DTRI_FACES);
  T8_ASSERT (0 <= face_child && face_child < T8_DTRI_FACE_CHILDREN);
  return t8_dtri_face_child_face ((const t8_dtri_t *) elem, face, face_child);
}

int
t8_default_scheme_tri_c::t8_element_face_parent_face (const t8_element_t *
                                                      elem, int face)
{
  T8_ASSERT (0 <= face && face < T8_DTRI_FACES);

  return t8_dtri_face_parent_face ((const t8_dtri_t *) elem, face);
}

int
t8_default_scheme_tri_c::t8_element_tree_face (const t8_element_t * elem,
                                               int face)
{
  T8_ASSERT (t8_element_is_valid (elem));
  T8_ASSERT (0 <= face && face < T8_DTRI_FACES);
  return t8_dtri_tree_face ((t8_dtri_t *) elem, face);
}

void
t8_default_scheme_tri_c::t8_element_transform_face (const t8_element_t *
                                                    elem1,
                                                    t8_element_t * elem2,
                                                    int orientation, int sign,
                                                    int is_smaller_face)
{
  T8_ASSERT (t8_element_is_valid (elem1));
  T8_ASSERT (t8_element_is_valid (elem2));
  t8_dtri_transform_face ((const t8_dtri_t *) elem1, (t8_dtri_t *) elem2,
                          orientation, sign, is_smaller_face);
}

/* Construct the inner element from a boundary element. */
/* This function is defined here instead of in t8_dri_bits.c since
 * the compile logic does not allow for t8_dtri_t and t8_dtet_t to exist
 * both in t8_dtri_bits.c. This would be needed by an implementation, at least
 * for tets. */
int
t8_default_scheme_tri_c::t8_element_extrude_face (const t8_element_t * face,
                                                  const t8_eclass_scheme_c *
                                                  face_scheme,
                                                  t8_element_t * elem,
                                                  int root_face)
{
  const t8_dline_t   *l = (const t8_dline_t *) face;
  t8_dtri_t          *t = (t8_dtri_t *) elem;

  T8_ASSERT (t8_element_is_valid (elem));
  T8_ASSERT (T8_COMMON_IS_TYPE
             (face_scheme, const t8_default_scheme_line_c *));
  T8_ASSERT (face_scheme->eclass == T8_ECLASS_LINE);
  T8_ASSERT (face_scheme->t8_element_is_valid (face));
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
   * Boundary triangles are alway of type 0.
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
t8_default_scheme_tri_c::t8_element_first_descendant_face (const t8_element_t
                                                           * elem, int face,
                                                           t8_element_t *
                                                           first_desc,
                                                           int level)
{
  int                 corner;
  T8_ASSERT (0 <= face && face < T8_DTRI_FACES);
  T8_ASSERT (0 <= level && level <= T8_DTRI_MAXLEVEL);

  /* Compute the first corner of this face */
  corner = t8_dtri_face_corner[face][0];
  /* Compute the descendant in this corner */
  t8_dtri_corner_descendant ((const t8_dtri_t *) elem,
                             (t8_dtri_t *) first_desc, corner, level);
}

void
t8_default_scheme_tri_c::t8_element_last_descendant_face (const t8_element_t *
                                                          elem, int face,
                                                          t8_element_t *
                                                          last_desc,
                                                          int level)
{
  int                 corner;
  T8_ASSERT (0 <= face && face < T8_DTRI_FACES);
  T8_ASSERT (0 <= level && level <= T8_DTRI_MAXLEVEL);
  /* Compute the last corner of this face */
  corner = t8_dtri_face_corner[face][1];
  /* Compute the descendant in this corner */
  t8_dtri_corner_descendant ((const t8_dtri_t *) elem,
                             (t8_dtri_t *) last_desc, corner, level);
}

/* Construct the boundary element at a specific face. */
void
t8_default_scheme_tri_c::t8_element_boundary_face (const t8_element_t * elem,
                                                   int face,
                                                   t8_element_t * boundary,
                                                   const t8_eclass_scheme_c *
                                                   boundary_scheme)
{
  const t8_dtri_t    *t = (const t8_dtri_t *) elem;
  t8_dline_t         *l = (t8_dline_t *) boundary;

  T8_ASSERT (t8_element_is_valid (elem));
  T8_ASSERT (T8_COMMON_IS_TYPE
             (boundary_scheme, const t8_default_scheme_line_c *));
  T8_ASSERT (boundary_scheme->eclass == T8_ECLASS_LINE);
  T8_ASSERT (boundary_scheme->t8_element_is_valid (boundary));
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

void
t8_default_scheme_tri_c::t8_element_boundary (const t8_element_t * elem,
                                              int min_dim, int length,
                                              t8_element_t ** boundary)
{

  SC_ABORT ("Not implemented\n");
#if 0
  int                 iface;

  T8_ASSERT (length == T8_DTRI_FACES);
  for (iface = 0; iface < T8_DTRI_FACES; iface++) {
    t8_element_boundary_face (elem, iface, boundary[iface]);
  }
#endif
}

int
t8_default_scheme_tri_c::t8_element_is_root_boundary (const t8_element_t *
                                                      elem, int face)
{
  const t8_dtri_t    *t = (const t8_dtri_t *) elem;
  T8_ASSERT (t8_element_is_valid (elem));

  return t8_dtri_is_root_boundary (t, face);
}

int
t8_default_scheme_tri_c::t8_element_face_neighbor_inside (const t8_element_t *
                                                          elem,
                                                          t8_element_t *
                                                          neigh, int face,
                                                          int *neigh_face)
{
  const t8_dtri_t    *t = (const t8_dtri_t *) elem;
  t8_dtri_t          *n = (t8_dtri_t *) neigh;

  T8_ASSERT (t8_element_is_valid (elem));
  T8_ASSERT (t8_element_is_valid (neigh));
  T8_ASSERT (0 <= face && face < T8_DTRI_FACES);
  T8_ASSERT (neigh_face != NULL);

  *neigh_face = t8_dtri_face_neighbour (t, face, n);
  /* return true if neigh is inside the root */
  return t8_dtri_is_inside_root (n);
}

void
t8_default_scheme_tri_c::t8_element_set_linear_id (t8_element_t * elem,
                                                   int level,
                                                   t8_linearidx_t id)
{
  T8_ASSERT (t8_element_is_valid (elem));
  T8_ASSERT (0 <= level && level <= T8_DTRI_MAXLEVEL);
  T8_ASSERT (0 <= id && id < ((t8_linearidx_t) 1) << (2 * level));

  t8_dtri_init_linear_id ((t8_dtri_t *) elem, id, level);
}

t8_linearidx_t
  t8_default_scheme_tri_c::t8_element_get_linear_id (const t8_element_t *
                                                     elem, int level)
{
  T8_ASSERT (t8_element_is_valid (elem));
  T8_ASSERT (0 <= level && level <= T8_DTRI_MAXLEVEL);

  return t8_dtri_linear_id ((t8_dtri_t *) elem, level);
}

void
t8_default_scheme_tri_c::t8_element_first_descendant (const t8_element_t *
                                                      elem,
                                                      t8_element_t * desc,
                                                      int level)
{
  T8_ASSERT (t8_element_is_valid (elem));
  T8_ASSERT (t8_element_is_valid (desc));
  T8_ASSERT (0 <= level && level <= T8_DTRI_MAXLEVEL);
  t8_dtri_first_descendant ((t8_dtri_t *) elem, (t8_dtri_t *) desc, level);
}

void
t8_default_scheme_tri_c::t8_element_last_descendant (const t8_element_t *
                                                     elem,
                                                     t8_element_t * desc,
                                                     int level)
{
  T8_ASSERT (t8_element_is_valid (elem));
  T8_ASSERT (t8_element_is_valid (desc));
  T8_ASSERT (0 <= level && level <= T8_DTRI_MAXLEVEL);
  t8_dtri_last_descendant ((t8_dtri_t *) elem, (t8_dtri_t *) desc, level);
}

void
t8_default_scheme_tri_c::t8_element_successor (const t8_element_t *
                                               elem1,
                                               t8_element_t * elem2,
                                               int level)
{
  T8_ASSERT (t8_element_is_valid (elem1));
  T8_ASSERT (t8_element_is_valid (elem2));
  T8_ASSERT (0 <= level && level <= T8_DTRI_MAXLEVEL);

  t8_dtri_successor ((const t8_dtri_t *) elem1, (t8_dtri_t *) elem2, level);
}

void
t8_default_scheme_tri_c::t8_element_anchor (const t8_element_t * elem,
                                            int anchor[3])
{
  t8_dtri_t          *tri = (t8_dtri_t *) elem;

  T8_ASSERT (t8_element_is_valid (elem));
  anchor[0] = tri->x;
  anchor[1] = tri->y;
  anchor[2] = 0;
}

int
t8_default_scheme_tri_c::t8_element_root_len (const t8_element_t * elem)
{
  T8_ASSERT (t8_element_is_valid (elem));
  return T8_DTRI_ROOT_LEN;
}

void
t8_default_scheme_tri_c::t8_element_vertex_coords (const t8_element_t * t,
                                                   int vertex, int coords[])
{
  T8_ASSERT (t8_element_is_valid (t));
  t8_dtri_compute_coords ((const t8_dtri_t *) t, vertex, coords);
}

void
t8_default_scheme_tri_c::t8_element_general_function (const t8_element_t *
                                                      elem,
                                                      const void *indata,
                                                      void *outdata)
{
  T8_ASSERT (t8_element_is_valid (elem));
  T8_ASSERT (outdata != NULL);
  *((int8_t *) outdata) = ((const t8_dtri_t *) elem)->type;
  /* Safety check to catch datatype conversion errors */
  T8_ASSERT (*((int8_t *) outdata) == ((const t8_dtri_t *) elem)->type);
}

#ifdef T8_ENABLE_DEBUG
/* *INDENT-OFF* */
/* indent bug, indent adds a second "const" modifier */
int
t8_default_scheme_tri_c::t8_element_is_valid (const t8_element_t * t) const
/* *INDENT-ON* */
{
  return t8_dtri_is_valid ((const t8_dtri_t *) t);
}
#endif

void
t8_default_scheme_tri_c::t8_element_new (int length, t8_element_t ** elem)
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
t8_default_scheme_tri_c::t8_element_init (int length, t8_element_t * elem,
                                          int new_called)
{
#ifdef T8_ENABLE_DEBUG
  if (!new_called) {
    int                 i;
    t8_dtri_t          *tris = (t8_dtri_t *) elem;
    for (i = 0; i < length; i++) {
      t8_dtri_init (tris + i);
    }
  }
#endif
}

/* Constructor */
t8_default_scheme_tri_c::t8_default_scheme_tri_c (void)
{
  eclass = T8_ECLASS_TRIANGLE;
  element_size = sizeof (t8_dtri_t);
  ts_context = sc_mempool_new (element_size);
}

/* Destructor */
t8_default_scheme_tri_c::~t8_default_scheme_tri_c ()
{
  /* This destructor is empty since the destructor of the
   * default_common scheme is called automatically and it
   * suffices to destroy the quad_scheme.
   * However we need to provide an implementation of the destructor
   * and hence this empty function. */
}

T8_EXTERN_C_END ();
