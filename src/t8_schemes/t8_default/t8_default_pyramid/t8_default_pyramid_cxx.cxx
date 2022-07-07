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
#include <t8_schemes/t8_default/t8_default_pyramid/t8_default_pyramid_cxx.hxx>
#include <t8_schemes/t8_default/t8_default_tet/t8_default_tet_cxx.hxx>
#include <t8_schemes/t8_default/t8_default_pyramid/t8_dpyramid_bits.h>
#include <t8_schemes/t8_default/t8_default_pyramid/t8_dpyramid.h>

typedef t8_dpyramid_t t8_default_pyramid_t;

T8_EXTERN_C_BEGIN ();

int
t8_default_scheme_pyramid_c::t8_element_maxlevel (void)
{
  return T8_DPYRAMID_MAXLEVEL;
}

int
t8_default_scheme_pyramid_c::t8_element_max_num_faces (const t8_element_t
                                                       *elem)
{
  T8_ASSERT (t8_element_is_valid (elem));
  return t8_dpyramid_max_num_faces ((const t8_dpyramid_t *) elem);
}

int
t8_default_scheme_pyramid_c::t8_element_ancestor_id (const t8_element_t *elem,
                                                     int level)
{
  T8_ASSERT (t8_element_is_valid (elem));
  T8_ASSERT (0 <= level && level <= T8_DPYRAMID_MAXLEVEL);
  return t8_dpyramid_ancestor_id ((const t8_dpyramid_t *) elem, level);
}

void
t8_default_scheme_pyramid_c::t8_element_init (int length, t8_element_t *elem,
                                              int called_new)
{
#ifdef T8_ENABLE_DEBUG
  if (!called_new) {
    int                 i;
    t8_dpyramid_t      *pyramid = (t8_dpyramid_t *) elem;
    /* Set all values to 0 */
    for (i = 0; i < length; i++) {
      t8_dpyramid_init_linear_id (pyramid + i, 0, 0);
      T8_ASSERT (t8_dpyramid_is_valid (pyramid + i));
    }
    pyramid->type = 6;
  }
#endif
}

int
t8_default_scheme_pyramid_c::t8_element_num_children (const t8_element_t
                                                      *elem)
{
  T8_ASSERT (t8_element_is_valid (elem));
  return t8_dpyramid_num_children ((const t8_dpyramid_t *) elem);
}

int
t8_default_scheme_pyramid_c::t8_element_num_corners (const t8_element_t *elem)
{
  T8_ASSERT (t8_element_is_valid (elem));
  return t8_dpyramid_num_corners ((const t8_dpyramid_t *) elem);
}

int
t8_default_scheme_pyramid_c::t8_element_num_siblings (const t8_element_t
                                                      *elem) const
{
  T8_ASSERT (t8_element_is_valid (elem));
  return t8_dpyramid_num_siblings ((const t8_dpyramid_t *) elem);
}

int
t8_default_scheme_pyramid_c::t8_element_num_faces (const t8_element_t *elem)
{
  T8_ASSERT (t8_element_is_valid (elem));
  return t8_dpyramid_num_faces ((const t8_dpyramid_t *) elem);
}

int
t8_default_scheme_pyramid_c::t8_element_compare (const t8_element_t *elem1,
                                                 const t8_element_t *elem2)
{
  T8_ASSERT (t8_element_is_valid (elem1));
  T8_ASSERT (t8_element_is_valid (elem2));
  return t8_dpyramid_compare ((const t8_dpyramid_t *) elem1,
                              (const t8_dpyramid_t *) elem2);
}

void
t8_default_scheme_pyramid_c::t8_element_copy (const t8_element_t *source,
                                              t8_element_t *dest)
{
  T8_ASSERT (t8_element_is_valid (source));
  t8_dpyramid_copy ((const t8_dpyramid_t *) source, (t8_dpyramid_t *) dest);
  T8_ASSERT (t8_element_is_valid (dest));
}

void
t8_default_scheme_pyramid_c::t8_element_child (const t8_element_t *elem,
                                               int childid,
                                               t8_element_t *child)
{
  T8_ASSERT (t8_element_is_valid (elem));
  T8_ASSERT (0 <= childid
             && childid <
             t8_dpyramid_num_children ((const t8_dpyramid_t *) elem));
  t8_dpyramid_child ((t8_dpyramid_t *) elem, childid,
                     (t8_dpyramid_t *) child);
  T8_ASSERT (t8_element_is_valid (child));
}

void
t8_default_scheme_pyramid_c::t8_element_children (const t8_element_t *elem,
                                                  int length,
                                                  t8_element_t *c[])
{
  T8_ASSERT (t8_element_is_valid (elem));
  t8_dpyramid_children ((const t8_dpyramid_t *) elem, (t8_dpyramid_t **) c);
#if T8_ENABLE_DEBUG
  for (int i = 0; i < length; i++) {
    T8_ASSERT (t8_element_is_valid (c[i]));
  }
#endif
}

void
t8_default_scheme_pyramid_c::t8_element_children_at_face (const t8_element_t
                                                          *elem, int face,
                                                          t8_element_t
                                                          *children[],
                                                          int num_children,
                                                          int *child_indices)
{
  T8_ASSERT (t8_element_is_valid (elem));
  t8_dpyramid_children_at_face ((const t8_dpyramid_t *) elem, face,
                                (t8_dpyramid_t **) children, num_children,
                                child_indices);
#if T8_ENABLE_DEBUG
  for (int i = 0; i < num_children; i++) {
    T8_ASSERT (t8_element_is_valid (children[i]));
  }
#endif
}

int
t8_default_scheme_pyramid_c::t8_element_face_child_face (const t8_element_t
                                                         *elem, int face,
                                                         int face_child)
{
  T8_ASSERT (t8_element_is_valid (elem));
  return t8_dpyramid_face_child_face ((const t8_dpyramid_t *) elem,
                                      face, face_child);
}

t8_element_shape_t
t8_default_scheme_pyramid_c::t8_element_face_shape (const t8_element_t *elem,
                                                    int face)
{
  T8_ASSERT (t8_element_is_valid (elem));
  return t8_dpyramid_face_shape ((const t8_dpyramid_t *) elem, face);
}

int
t8_default_scheme_pyramid_c::t8_element_child_id (const t8_element_t *p)
{
  T8_ASSERT (t8_element_is_valid (p));
  return t8_dpyramid_child_id ((const t8_dpyramid_t *) p);
}

int
t8_default_scheme_pyramid_c::t8_element_face_neighbor_inside (const
                                                              t8_element_t
                                                              *elem,
                                                              t8_element_t
                                                              *neigh,
                                                              int face,
                                                              int *neigh_face)
{
  T8_ASSERT (t8_element_is_valid (elem));
  return t8_dpyramid_face_neighbor_inside ((const t8_dpyramid_t *) elem,
                                           (t8_dpyramid_t *) neigh,
                                           face, neigh_face);
}

int
t8_default_scheme_pyramid_c::t8_element_face_parent_face (const t8_element_t
                                                          *elem, int face)
{
  T8_ASSERT (t8_element_is_valid (elem));
  return t8_dpyramid_face_parent_face ((const t8_dpyramid_t *) elem, face);
}

void
t8_default_scheme_pyramid_c::t8_element_first_descendant (const t8_element_t
                                                          *elem,
                                                          t8_element_t *desc,
                                                          int level)
{
  T8_ASSERT (t8_element_is_valid (elem));
  t8_dpyramid_first_descendant ((const t8_dpyramid_t *) elem,
                                (t8_dpyramid_t *) desc, level);
  T8_ASSERT (t8_element_is_valid (desc));
}

void
t8_default_scheme_pyramid_c::t8_element_first_descendant_face (const
                                                               t8_element_t
                                                               *elem,
                                                               int face,
                                                               t8_element_t
                                                               *first_desc,
                                                               int level)
{
  T8_ASSERT (t8_element_is_valid (elem));
  t8_dpyramid_first_descendant_face ((const t8_dpyramid_t *) elem, face,
                                     (t8_dpyramid_t *) first_desc, level);
  T8_ASSERT (t8_element_is_valid (first_desc));
}

int
t8_default_scheme_pyramid_c::t8_element_get_face_corner (const
                                                         t8_element_t
                                                         *element, int face,
                                                         int corner)
{
  T8_ASSERT (t8_element_is_valid (element));
  return t8_dpyramid_get_face_corner ((const t8_dpyramid_t *) element,
                                      face, corner);
}

int
t8_default_scheme_pyramid_c::t8_element_level (const t8_element_t *elem)
{
  T8_ASSERT (t8_element_is_valid (elem));
  return t8_dpyramid_get_level ((const t8_dpyramid_t *) elem);
}

int
t8_default_scheme_pyramid_c::t8_element_root_len (const t8_element_t *elem)
{
  T8_ASSERT (t8_element_is_valid (elem));
  return T8_DPYRAMID_ROOT_LEN;
}

void
t8_default_scheme_pyramid_c::t8_element_set_linear_id (t8_element_t *elem,
                                                       int level, uint64_t id)
{
  t8_dpyramid_init_linear_id ((t8_dpyramid_t *) elem, level, id);
  T8_ASSERT (t8_element_is_valid (elem));
}

int
t8_default_scheme_pyramid_c::t8_element_is_family (t8_element_t **fam)
{
#if T8_ENABLE_DEBUG
  int                 num_siblings = t8_element_num_siblings (fam[0]);
  for (int i = 0; i < num_siblings; i++) {
    T8_ASSERT (t8_element_is_valid (fam[i]));
  }
#endif
  return t8_dpyramid_is_family ((const t8_dpyramid_t **) fam);
}

int
t8_default_scheme_pyramid_c::t8_element_is_root_boundary (const t8_element_t
                                                          *elem, int face)
{
  T8_ASSERT (t8_element_is_valid (elem));
  return t8_dpyramid_is_root_boundary ((const t8_dpyramid_t *) elem, face);
}

void
t8_default_scheme_pyramid_c::t8_element_boundary_face (const t8_element_t
                                                       *elem, int face,
                                                       t8_element_t *boundary,
                                                       const
                                                       t8_eclass_scheme_c
                                                       *boundary_scheme)
{
  T8_ASSERT (t8_element_is_valid (elem));
  t8_dpyramid_boundary_face ((const t8_dpyramid_t *) elem, face, boundary);
  T8_ASSERT (boundary_scheme->t8_element_is_valid (boundary));
}

int
t8_default_scheme_pyramid_c::t8_element_extrude_face (const t8_element_t
                                                      *face,
                                                      const t8_eclass_scheme_c
                                                      *face_scheme,
                                                      t8_element_t *elem,
                                                      int root_face)
{
  T8_ASSERT (face_scheme->t8_element_is_valid (face));
  return t8_dpyramid_extrude_face (face, (t8_dpyramid_t *) elem, root_face);
  T8_ASSERT (t8_element_is_valid (elem));
}

t8_element_shape_t
t8_default_scheme_pyramid_c::t8_element_shape (const t8_element_t *elem)
{
  T8_ASSERT (t8_element_is_valid (elem));
  return t8_dpyramid_shape ((const t8_dpyramid_t *) elem);
}

int
t8_default_scheme_pyramid_c::t8_element_tree_face (const t8_element_t *elem,
                                                   int face)
{
  T8_ASSERT (t8_element_is_valid (elem));
  return t8_dpyramid_tree_face ((const t8_dpyramid_t *) elem, face);
}

int
t8_default_scheme_pyramid_c::t8_element_num_face_children (const t8_element_t
                                                           *elem, int face)
{
  T8_ASSERT (t8_element_is_valid (elem));
  return T8_DPYRAMID_FACE_CHILDREN;
}

t8_linearidx_t
  t8_default_scheme_pyramid_c::t8_element_get_linear_id (const t8_element_t
                                                         *elem, int level)
{
  T8_ASSERT (t8_element_is_valid (elem));
  return t8_dpyramid_linear_id ((const t8_dpyramid_t *) elem, level);
}

void
t8_default_scheme_pyramid_c::t8_element_last_descendant (const t8_element_t
                                                         *elem,
                                                         t8_element_t *desc,
                                                         int level)
{
  T8_ASSERT (t8_element_is_valid (elem));
  t8_dpyramid_last_descendant ((const t8_dpyramid_t *) elem,
                               (t8_dpyramid_t *) desc, level);
  T8_ASSERT (t8_element_is_valid (desc));
}

void
t8_default_scheme_pyramid_c::t8_element_last_descendant_face (const
                                                              t8_element_t
                                                              *elem, int face,
                                                              t8_element_t
                                                              *last_desc,
                                                              int level)
{
  T8_ASSERT (t8_element_is_valid (elem));
  t8_dpyramid_last_descendant_face ((const t8_dpyramid_t *) elem,
                                    face, (t8_dpyramid_t *) last_desc, level);
  T8_ASSERT (t8_element_is_valid (last_desc));
}

void
t8_default_scheme_pyramid_c::t8_element_parent (const t8_element_t *elem,
                                                t8_element_t *parent)
{
  T8_ASSERT (t8_element_is_valid (elem));
  t8_dpyramid_parent ((const t8_dpyramid_t *) elem, (t8_dpyramid_t *) parent);
  T8_ASSERT (t8_element_is_valid (parent));
}

void
t8_default_scheme_pyramid_c::t8_element_successor (const t8_element_t *elem,
                                                   t8_element_t *s, int level)
{
  T8_ASSERT (t8_element_is_valid (elem));
  t8_dpyramid_successor ((const t8_dpyramid_t *) elem, (t8_dpyramid_t *) s,
                         level);
  T8_ASSERT (t8_element_is_valid (s));
}

void
t8_default_scheme_pyramid_c::t8_element_vertex_coords (const t8_element_t *t,
                                                       int vertex,
                                                       int coords[])
{
  T8_ASSERT (t8_element_is_valid (t));
  t8_dpyramid_compute_coords ((const t8_dpyramid_t *) t, vertex, coords);
}

void
t8_default_scheme_pyramid_c::t8_element_vertex_reference_coords (const
                                                                 t8_element_t
                                                                 *elem,
                                                                 int vertex,
                                                                 double
                                                                 coords[])
{
  T8_ASSERT (t8_element_is_valid (elem));
  t8_dpyramid_vertex_reference_coords ((const t8_dpyramid_t *) elem,
                                       vertex, coords);
}

int
t8_default_scheme_pyramid_c::t8_element_refines_irregular ()
{
  /*Pyramids do not refine regularly */
  return 1;
}

void
t8_default_scheme_pyramid_c::t8_element_general_function (const t8_element_t
                                                          *elem,
                                                          const void *indata,
                                                          void *outdata)
{
  T8_ASSERT (outdata != NULL);
  T8_ASSERT (t8_element_is_valid (elem));
  *((int8_t *) outdata) = ((const t8_dpyramid_t *) elem)->type;
  /* Safety check to catch datatype conversion errors */
  T8_ASSERT (*((int8_t *) outdata) == ((const t8_dpyramid_t *) elem)->type);
}

#ifdef T8_ENABLE_DEBUG
int
t8_default_scheme_pyramid_c::t8_element_is_valid (const t8_element_t *elem) const
{
  T8_ASSERT (elem != NULL);
  return t8_dpyramid_is_valid ((const t8_dpyramid_t *) elem);
}

void
t8_default_scheme_pyramid_c::t8_element_debug_print (const t8_element_t *elem) const
{
  T8_ASSERT (t8_element_is_valid (elem));
  t8_dpyramid_debug_print ((const t8_dpyramid_t *) elem);
}
#endif

/* Constructor */
t8_default_scheme_pyramid_c::t8_default_scheme_pyramid_c (void)
{
  eclass = T8_ECLASS_PYRAMID;
  element_size = sizeof (t8_default_pyramid_t);
  ts_context = sc_mempool_new (element_size);
}

t8_default_scheme_pyramid_c::~t8_default_scheme_pyramid_c ()
{
  /* This destructor is empty since the destructor of the
   * default_common scheme is called automatically and it
   * suffices to destroy the quad_scheme.
   * However we need to provide an implementation of the destructor
   * and hence this empty function. */
}

T8_EXTERN_C_END ();
