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

#include <t8_schemes/t8_default/t8_default_common/t8_default_common_cxx.hxx>
#include <t8_schemes/t8_default/t8_default_prism/t8_default_prism_cxx.hxx>
#include <t8_schemes/t8_default/t8_default_prism/t8_dprism_bits.h>
#include <t8_schemes/t8_default/t8_default_prism/t8_dprism.h>

typedef t8_dprism_t t8_default_prism_t;

T8_EXTERN_C_BEGIN ();

void
t8_default_scheme_prism_c::t8_element_new (int length, t8_element_t **elem)
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
t8_default_scheme_prism_c::t8_element_init (int length, t8_element_t *elem,
                                            int new_called)
{
#ifdef T8_ENABLE_DEBUG
  if (!new_called) {
    int                 i;
    t8_dprism_t        *prism = (t8_dprism_t *) elem;
    /* Set all values to 0 */
    for (i = 0; i < length; i++) {
      t8_dprism_init_linear_id (prism + i, 0, 0);
      T8_ASSERT (t8_dprism_is_valid (prism + i));
    }
  }
#endif
}

int
t8_default_scheme_prism_c::t8_element_maxlevel (void)
{
  return T8_DPRISM_MAXLEVEL;
}

int
t8_default_scheme_prism_c::t8_element_level (const t8_element_t *elem)
{
  T8_ASSERT (t8_element_is_valid (elem));
  return t8_dprism_get_level ((const t8_dprism_t *) elem);
}

t8_element_shape_t
t8_default_scheme_prism_c::t8_element_face_shape (const t8_element_t *elem,
                                                  int face)
{
  T8_ASSERT (t8_element_is_valid (elem));
  T8_ASSERT (0 <= face && face < T8_DPRISM_FACES);

  return t8_dprism_face_shape ((const t8_dprism_t *) elem, face);
}

void
t8_default_scheme_prism_c::t8_element_copy (const t8_element_t *source,
                                            t8_element_t *dest)
{
  T8_ASSERT (t8_element_is_valid (source));
  t8_dprism_copy ((const t8_dprism_t *) source, (t8_dprism_t *) dest);
  T8_ASSERT (t8_element_is_valid (dest));
}

int
t8_default_scheme_prism_c::t8_element_compare (const t8_element_t *elem1,
                                               const t8_element_t *elem2)
{
  T8_ASSERT (t8_element_is_valid (elem1));
  T8_ASSERT (t8_element_is_valid (elem2));
  return t8_dprism_compare ((const t8_dprism_t *) elem1,
                            (const t8_dprism_t *) elem2);
}

void
t8_default_scheme_prism_c::t8_element_parent (const t8_element_t *elem,
                                              t8_element_t *parent)
{
  T8_ASSERT (t8_element_is_valid (elem));
  t8_dprism_parent ((const t8_dprism_t *) elem, (t8_dprism_t *) parent);
  T8_ASSERT (t8_element_is_valid (parent));
}

int
t8_default_scheme_prism_c::t8_element_num_children (const t8_element_t *elem)
{
  T8_ASSERT (t8_element_is_valid (elem));
  return T8_DPRISM_CHILDREN;
}

int
t8_default_scheme_prism_c::t8_element_num_face_children (const t8_element_t
                                                         *elem, int face)
{
  T8_ASSERT (t8_element_is_valid (elem));
  return t8_dprism_num_face_children ((const t8_dprism_t *) elem, face);
}

int
t8_default_scheme_prism_c::t8_element_get_face_corner (const t8_element_t
                                                       *element, int face,
                                                       int corner)
{
  T8_ASSERT (t8_element_is_valid (element));
  T8_ASSERT (0 <= face && face < T8_DPRISM_FACES);
  T8_ASSERT (0 <= corner && corner < T8_DPRISM_CORNERS);

  return t8_dprism_get_face_corner ((const t8_dprism_t *) element, face,
                                    corner);
}

int
t8_default_scheme_prism_c::t8_element_num_faces (const t8_element_t *elem)
{
  T8_ASSERT (t8_element_is_valid (elem));
  return T8_DPRISM_FACES;
}

int
t8_default_scheme_prism_c::t8_element_child_id (const t8_element_t *elem)
{
  T8_ASSERT (t8_element_is_valid (elem));
  return t8_dprism_child_id ((const t8_dprism_t *) elem);
}

void
t8_default_scheme_prism_c::t8_element_child (const t8_element_t *elem,
                                             int childid, t8_element_t *child)
{
  T8_ASSERT (t8_element_is_valid (elem));
  t8_dprism_child ((const t8_dprism_t *) elem, childid,
                   (t8_dprism_t *) child);
  T8_ASSERT (t8_element_is_valid (child));
}

int
t8_default_scheme_prism_c::t8_element_max_num_faces (const t8_element_t *elem)
{
  T8_ASSERT (t8_element_is_valid (elem));
  return T8_DPRISM_FACES;
}

void
t8_default_scheme_prism_c::t8_element_children (const t8_element_t *elem,
                                                int length, t8_element_t *c[])
{
  T8_ASSERT (t8_element_is_valid (elem));
  T8_ASSERT (length == T8_DPRISM_CHILDREN);
  t8_dprism_childrenpv ((const t8_dprism_t *) elem, length,
                        (t8_dprism_t **) c);
#if T8_ENABLE_DEBUG
  for (int i = 0; i < length; i++) {
    T8_ASSERT (t8_element_is_valid (c[i]));
  }
#endif
}

int
t8_default_scheme_prism_c::t8_element_ancestor_id (const t8_element_t *elem,
                                                   int level)
{
  T8_ASSERT (t8_element_is_valid (elem));
  return t8_dprism_ancestor_id ((t8_dprism_t *) elem, level);
}

void
t8_default_scheme_prism_c::t8_element_children_at_face (const t8_element_t
                                                        *elem, int face,
                                                        t8_element_t
                                                        *children[],
                                                        int num_children,
                                                        int *child_indices)
{
  T8_ASSERT (t8_element_is_valid (elem));
  T8_ASSERT (0 <= face && face < T8_DPRISM_FACES);
  T8_ASSERT (num_children ==
             t8_dprism_num_face_children ((const t8_dprism_t *) elem, face));
  t8_dprism_children_at_face ((const t8_dprism_t *) elem, face,
                              (t8_dprism_t **) children, num_children);
  /* TODO: Properly implement child_indices */
  T8_ASSERT (child_indices == NULL);
#if T8_ENABLE_DEBUG
  for (int i = 0; i < num_children; i++) {
    T8_ASSERT (t8_element_is_valid (children[i]));
  }
#endif
}

int
t8_default_scheme_prism_c::t8_element_face_child_face (const t8_element_t
                                                       *elem, int face,
                                                       int face_child)
{
  T8_ASSERT (t8_element_is_valid (elem));
  T8_ASSERT (0 <= face && face < T8_DPRISM_FACES);
  T8_ASSERT (face_child <
             t8_dprism_num_face_children ((const t8_dprism_t *) elem, face));
  return t8_dprism_face_child_face ((const t8_dprism_t *) elem, face,
                                    face_child);
}

int
t8_default_scheme_prism_c::t8_element_face_parent_face (const t8_element_t
                                                        *elem, int face)
{
  T8_ASSERT (t8_element_is_valid (elem));
  T8_ASSERT (0 <= face && face < T8_DPRISM_FACES);
  return t8_dprism_face_parent_face ((const t8_dprism_t *) elem, face);
}

int
t8_default_scheme_prism_c::t8_element_tree_face (const t8_element_t *elem,
                                                 int face)
{
  T8_ASSERT (t8_element_is_valid (elem));
  T8_ASSERT (0 <= face && face < T8_DPRISM_FACES);
  return t8_dprism_tree_face ((const t8_dprism_t *) elem, face);
}

int
t8_default_scheme_prism_c::t8_element_extrude_face (const t8_element_t *face,
                                                    const t8_eclass_scheme_c
                                                    *face_scheme,
                                                    t8_element_t *elem,
                                                    int root_face)
{
  T8_ASSERT (face_scheme->t8_element_is_valid (face));
  T8_ASSERT (0 <= root_face && root_face < T8_DPRISM_FACES);
  t8_dprism_extrude_face (face, elem, root_face);
  T8_ASSERT (t8_element_is_valid (elem));
  /* TODO: Fix return value */
  return t8_dprism_root_face_to_face ((const t8_dprism_t *) elem, root_face);
}

int
t8_default_scheme_prism_c::t8_element_is_family (t8_element_t **fam)
{
#ifdef T8_ENABLE_DEBUG
  int                 i;
  for (i = 0; i < T8_DPRISM_CHILDREN; i++) {
    T8_ASSERT (t8_element_is_valid (fam[i]));
  }
#endif
  return t8_dprism_is_familypv ((t8_dprism_t **) fam);
}

void
t8_default_scheme_prism_c::t8_element_nca (const t8_element_t *elem1,
                                           const t8_element_t *elem2,
                                           t8_element_t *nca)
{
  T8_ASSERT (t8_element_is_valid (elem1));
  T8_ASSERT (t8_element_is_valid (elem2));
  const t8_default_prism_t *p1 = (const t8_default_prism_t *) elem1;
  const t8_default_prism_t *p2 = (const t8_default_prism_t *) elem2;
  t8_default_prism_t *c = (t8_default_prism_t *) nca;

  t8_dprism_nearest_common_ancestor (p1, p2, c);
  T8_ASSERT (t8_dprism_is_valid (c));
}

void
t8_default_scheme_prism_c::t8_element_boundary_face (const t8_element_t *elem,
                                                     int face,
                                                     t8_element_t *boundary,
                                                     const t8_eclass_scheme_c
                                                     *boundary_scheme)
{
  T8_ASSERT (t8_element_is_valid (elem));
  T8_ASSERT (0 <= face && face < T8_DPRISM_FACES);
  t8_dprism_boundary_face ((const t8_dprism_t *) elem, face, boundary);
  T8_ASSERT (boundary_scheme->t8_element_is_valid (boundary));
}

const int           t8_dprism_face_corner[5][4] = {
  {1, 2, 4, 5},
  {0, 2, 3, 5},
  {0, 1, 3, 4},
  {0, 1, 2, -1},                /*this face is a triangle -> -1 for the 4th corner */
  {3, 4, 5, -1}                 /*this face is a triangle -> -1 for the 4th corner */
};

void
t8_default_scheme_prism_c::t8_element_first_descendant_face (const
                                                             t8_element_t
                                                             *elem, int face,
                                                             t8_element_t
                                                             *first_desc,
                                                             int level)
{
  int                 corner;
  T8_ASSERT (0 <= face && face < T8_DPRISM_FACES);
  T8_ASSERT (0 <= level && level <= T8_DPRISM_MAXLEVEL);
  T8_ASSERT (t8_element_is_valid (elem));
  corner = t8_dprism_face_corner[face][0];
  t8_dprism_corner_descendant ((const t8_dprism_t *) elem,
                               (t8_dprism_t *) first_desc, corner, level);
  T8_ASSERT (t8_element_is_valid (first_desc));
}

void
t8_default_scheme_prism_c::t8_element_last_descendant_face (const t8_element_t
                                                            *elem, int face,
                                                            t8_element_t
                                                            *last_desc,
                                                            int level)
{
  int                 corner;
  T8_ASSERT (0 <= face && face < T8_DPRISM_FACES);
  T8_ASSERT (0 <= level && level <= T8_DPRISM_MAXLEVEL);
  T8_ASSERT (t8_element_is_valid (elem));
  if (face < 3) {
    corner = t8_dprism_face_corner[face][3];
  }
  else {
    corner = t8_dprism_face_corner[face][2];
  }
  t8_dprism_corner_descendant ((const t8_dprism_t *) elem,
                               (t8_dprism_t *) last_desc, corner, level);
  T8_ASSERT (t8_element_is_valid (last_desc));
}

int
t8_default_scheme_prism_c::t8_element_is_root_boundary (const t8_element_t
                                                        *elem, int face)
{
  T8_ASSERT (t8_element_is_valid (elem));
  return t8_dprism_is_root_boundary ((const t8_dprism_t *) elem, face);
}

int
t8_default_scheme_prism_c::t8_element_face_neighbor_inside (const t8_element_t
                                                            *elem,
                                                            t8_element_t
                                                            *neigh, int face,
                                                            int *neigh_face)
{
  const t8_dprism_t  *p = (const t8_dprism_t *) elem;
  t8_dprism_t        *n = (t8_dprism_t *) neigh;

  T8_ASSERT (0 <= face && face < T8_DPRISM_FACES);
  T8_ASSERT (neigh_face != NULL);
  T8_ASSERT (t8_element_is_valid (elem));
  T8_ASSERT (t8_element_is_valid (neigh));

  *neigh_face = t8_dprism_face_neighbour (p, face, n);
  /* return true if neigh is inside the root */
  return t8_dprism_is_inside_root (n);
}

void
t8_default_scheme_prism_c::t8_element_set_linear_id (t8_element_t *elem,
                                                     int level, uint64_t id)
{
  T8_ASSERT (0 <= level && level <= T8_DPRISM_MAXLEVEL);
  T8_ASSERT (0 <= id && id < ((u_int64_t) 1) << 3 * level);

  t8_dprism_init_linear_id ((t8_default_prism_t *) elem, level, id);

  T8_ASSERT (t8_element_is_valid (elem));
}

void
t8_default_scheme_prism_c::t8_element_successor (const t8_element_t *elem,
                                                 t8_element_t *s, int level)
{
  T8_ASSERT (1 <= level && level <= T8_DPRISM_MAXLEVEL);
  T8_ASSERT (t8_element_is_valid (elem));

  t8_dprism_successor ((const t8_default_prism_t *) elem,
                       (t8_default_prism_t *) s, level);
  T8_ASSERT (t8_element_is_valid (s));
}

void
t8_default_scheme_prism_c::t8_element_first_descendant (const t8_element_t
                                                        *elem,
                                                        t8_element_t *desc,
                                                        int level)
{
  T8_ASSERT (0 <= level && level <= T8_DPRISM_MAXLEVEL);
  T8_ASSERT (t8_element_is_valid (elem));
  t8_dprism_first_descendant ((const t8_default_prism_t *) elem,
                              (t8_default_prism_t *) desc, level);
  T8_ASSERT (t8_element_is_valid (desc));
}

void
t8_default_scheme_prism_c::t8_element_last_descendant (const t8_element_t
                                                       *elem,
                                                       t8_element_t *desc,
                                                       int level)
{
  T8_ASSERT (0 <= level && level <= T8_DPRISM_MAXLEVEL);
  T8_ASSERT (t8_element_is_valid (elem));
  t8_dprism_last_descendant ((const t8_default_prism_t *) elem,
                             (t8_default_prism_t *) desc, level);
  T8_ASSERT (t8_element_is_valid (desc));
}

void
t8_default_scheme_prism_c::t8_element_anchor (const t8_element_t *elem,
                                              int anchor[3])
{
  t8_dprism_t        *prism = (t8_dprism_t *) elem;
  T8_ASSERT (t8_element_is_valid (elem));
  anchor[0] = prism->tri.x / T8_DTRI_ROOT_LEN * T8_DPRISM_ROOT_LEN;
  anchor[1] = prism->tri.y / T8_DTRI_ROOT_LEN * T8_DPRISM_ROOT_LEN;
  anchor[2] = prism->line.x / T8_DLINE_ROOT_LEN * T8_DPRISM_ROOT_LEN;
}

int
t8_default_scheme_prism_c::t8_element_root_len (const t8_element_t *elem)
{
  T8_ASSERT (t8_element_is_valid (elem));
  return T8_DPRISM_ROOT_LEN;
}

void
t8_default_scheme_prism_c::t8_element_vertex_coords (const t8_element_t *elem,
                                                     int vertex, int coords[])
{
  T8_ASSERT (t8_element_is_valid (elem));
  t8_dprism_vertex_coords ((const t8_dprism_t *) elem, vertex, coords);
}

void
t8_default_scheme_prism_c::t8_element_vertex_reference_coords (const
                                                               t8_element_t
                                                               *elem,
                                                               int vertex,
                                                               double
                                                               coords[])
{
  T8_ASSERT (t8_element_is_valid (elem));
  t8_dprism_vertex_ref_coords ((const t8_dprism_t *) elem, vertex, coords);
}

void
t8_default_scheme_prism_c::t8_element_general_function (const t8_element_t
                                                        *elem,
                                                        const void *indata,
                                                        void *outdata)
{
  T8_ASSERT (outdata != NULL);
  T8_ASSERT (t8_element_is_valid (elem));
  *((int8_t *) outdata) = ((const t8_dprism_t *) elem)->tri.type;
  /* Safety check to catch datatype conversion errors */
  T8_ASSERT (*((int8_t *) outdata) == ((const t8_dprism_t *) elem)->tri.type);
}

u_int64_t
  t8_default_scheme_prism_c::t8_element_get_linear_id (const t8_element_t
                                                       *elem, int level)
{
  T8_ASSERT (t8_element_is_valid (elem));
  return t8_dprism_linear_id ((const t8_dprism_t *) elem, level);
}

int
t8_default_scheme_prism_c::t8_element_refines_irregular (void)
{
  /*prisms refine regularly */
  return 0;
}

#ifdef T8_ENABLE_DEBUG
/* *INDENT-OFF* */
/* Indent bug, indent adds an additional const modifier at the end */
int
t8_default_scheme_prism_c::t8_element_is_valid (const t8_element_t * elem) const
{
  T8_ASSERT (elem != NULL);
  return t8_dprism_is_valid ((const t8_dprism_t *) elem);
}
/* *INDENT-ON* */

void
t8_default_scheme_prism_c::t8_element_debug_print (const t8_element_t *elem) const
{
  T8_ASSERT (t8_element_is_valid (elem));
  t8_dprism_debug_print ((const t8_dprism_t *) elem);
}
#endif /* T8_ENABLE_DEBUG */

/* Constructor */
t8_default_scheme_prism_c::t8_default_scheme_prism_c (void)
{
  eclass = T8_ECLASS_PRISM;
  element_size = sizeof (t8_default_prism_t);
  ts_context = sc_mempool_new (element_size);
}

t8_default_scheme_prism_c::~t8_default_scheme_prism_c ()
{
  /* This destructor is empty since the destructor of the
   * default_common scheme is called automatically and it
   * suffices to destroy the quad_scheme.
   * However we need to provide an implementation of the destructor
   * and hence this empty function. */
}

T8_EXTERN_C_END ();
