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

#include "t8_standalone_element_cxx.hxx"

template <t8_eclass_t eclass_T>
t8_standalone_scheme_c<eclass_T>::t8_standalone_scheme_c (void)
{
  eclass = eclass_T;
  element_size = sizeof (t8_standalone_element_t<eclass_T>);
  ts_context = sc_mempool_new (element_size);
}

template <t8_eclass_t eclass_T>
t8_standalone_scheme_c<eclass_T>::~t8_standalone_scheme_c ()
{
  T8_ASSERT (ts_context != NULL);
  sc_mempool_destroy ((sc_mempool_t *) ts_context);
}

template <t8_eclass_t eclass_T>
int
t8_standalone_scheme_c<eclass_T>::t8_element_refines_irregular (void) const
{
  if constexpr (eclass_T == T8_ECLASS_PYRAMID) {
    return 1;
  }
  return 0;
}

template <t8_eclass_t eclass_T>
int
t8_standalone_scheme_c<eclass_T>::t8_element_maxlevel (void) const
{
  return T8_ELEMENT_MAXLEVEL[eclass_T];
}

template <t8_eclass_t eclass_T>
t8_eclass_t
t8_standalone_scheme_c<eclass_T>::t8_element_child_eclass (int childid) const
{
  SC_ABORT ("This function should not be needed with the current approach to eclass vs shape.\n");
  return T8_ECLASS_ZERO; /* suppresses compiler warning */
}

template <t8_eclass_t eclass_T>
int
t8_standalone_scheme_c<eclass_T>::t8_element_level (const t8_element_t *elem) const
{
  T8_ASSERT (t8_element_is_valid (elem));
  return ((const t8_standalone_element_t<eclass_T> *) elem)->level;
}

template <t8_eclass_t eclass_T>
void
t8_standalone_scheme_c<eclass_T>::t8_element_copy (const t8_element_t *source, t8_element_t *dest) const
{
  T8_ASSERT (t8_element_is_valid (source));
  t8_sele_copy<eclass_T> ((const t8_standalone_element_t<eclass_T> *) source,
                          (t8_standalone_element_t<eclass_T> *) dest);
  T8_ASSERT (t8_element_is_valid (dest));
}

template <t8_eclass_t eclass_T>
int
t8_standalone_scheme_c<eclass_T>::t8_element_compare (const t8_element_t *elem1, const t8_element_t *elem2) const
{
  T8_ASSERT (t8_element_is_valid (elem1));
  T8_ASSERT (t8_element_is_valid (elem2));
  return t8_sele_compare<eclass_T> ((const t8_standalone_element_t<eclass_T> *) elem1,
                                    (const t8_standalone_element_t<eclass_T> *) elem2);
}

template <t8_eclass_t eclass_T>
int
t8_standalone_scheme_c<eclass_T>::t8_element_equal (const t8_element_t *elem1, const t8_element_t *elem2) const
{
  T8_ASSERT (t8_element_is_valid (elem1));
  T8_ASSERT (t8_element_is_valid (elem2));
  return t8_sele_equal<eclass_T> ((const t8_standalone_element_t<eclass_T> *) elem1,
                                  (const t8_standalone_element_t<eclass_T> *) elem2);
}



template <t8_eclass_t eclass_T>
void
t8_standalone_scheme_c<eclass_T>::t8_element_parent (const t8_element_t *elem, t8_element_t *parent) const
{
  T8_ASSERT (t8_element_is_valid (elem));
  t8_sele_parent<eclass_T> ((const t8_standalone_element_t<eclass_T> *) elem,
                            (t8_standalone_element_t<eclass_T> *) parent);
  T8_ASSERT (t8_element_is_valid (parent));
}

template <t8_eclass_t eclass_T>
int
t8_standalone_scheme_c<eclass_T>::t8_element_num_siblings (const t8_element_t *elem) const
{
  T8_ASSERT (t8_element_is_valid (elem));
  return t8_sele_num_siblings<eclass_T> ((const t8_standalone_element_t<eclass_T> *) elem);
}

template <t8_eclass_t eclass_T>
void
t8_standalone_scheme_c<eclass_T>::t8_element_sibling (const t8_element_t *elem, int sibid, t8_element_t *sibling) const
{
  SC_ABORT ("This function is not implemented yet.\n");
}

template <t8_eclass_t eclass_T>
int
t8_standalone_scheme_c<eclass_T>::t8_element_num_corners (const t8_element_t *elem) const
{
  T8_ASSERT (t8_element_is_valid (elem));
  return t8_sele_num_corners<eclass_T> ((const t8_standalone_element_t<eclass_T> *) elem);
}

template <t8_eclass_t eclass_T>
int
t8_standalone_scheme_c<eclass_T>::t8_element_num_faces (const t8_element_t *elem) const
{
  T8_ASSERT (t8_element_is_valid (elem));
  return t8_sele_num_faces<eclass_T> ((const t8_standalone_element_t<eclass_T> *) elem);
}

template <t8_eclass_t eclass_T>
int
t8_standalone_scheme_c<eclass_T>::t8_element_max_num_faces (const t8_element_t *elem) const
{
  T8_ASSERT (t8_element_is_valid (elem));
  return t8_sele_max_num_faces<eclass_T> ((const t8_standalone_element_t<eclass_T> *) elem);
}

template <t8_eclass_t eclass_T>
int
t8_standalone_scheme_c<eclass_T>::t8_element_num_children (const t8_element_t *elem) const
{
  T8_ASSERT (t8_element_is_valid (elem));
  return t8_sele_num_children<eclass_T> ((const t8_standalone_element_t<eclass_T> *) elem);
}

template <t8_eclass_t eclass_T>
int
t8_standalone_scheme_c<eclass_T>::t8_element_num_face_children (const t8_element_t *elem, int face) const
{
  T8_ASSERT (t8_element_is_valid (elem));
  // Not true anymore for 4D with pyramids as faces
  if constexpr (eclass_T == T8_ECLASS_VERTEX) {
    return 0;
  }
  else {
    return 1 << (T8_ELEMENT_DIM[eclass_T] - 1);
  }
}

template <t8_eclass_t eclass_T>
int
t8_standalone_scheme_c<eclass_T>::t8_element_get_face_corner (const t8_element_t *element, int face, int corner) const
{
  SC_ABORT ("This function is not implemented yet.\n");
  return 0;
}

template <t8_eclass_t eclass_T>
int
t8_standalone_scheme_c<eclass_T>::t8_element_get_corner_face (const t8_element_t *element, int corner, int face) const
{
  SC_ABORT ("This function is not implemented yet.\n");
  return 0;
}

template <t8_eclass_t eclass_T>
void
t8_standalone_scheme_c<eclass_T>::t8_element_child (const t8_element_t *elem, int childid, t8_element_t *child) const
{
  T8_ASSERT (t8_element_is_valid (elem));
  T8_ASSERT (0 <= childid);
  T8_ASSERT (childid < t8_sele_num_children<eclass_T> ((const t8_standalone_element_t<eclass_T> *) elem));
  t8_sele_child<eclass_T> ((t8_standalone_element_t<eclass_T> *) elem, childid,
                           (t8_standalone_element_t<eclass_T> *) child);
  T8_ASSERT (t8_element_is_valid (child));
}

template <t8_eclass_t eclass_T>
void
t8_standalone_scheme_c<eclass_T>::t8_element_children (const t8_element_t *elem, int length, t8_element_t *c[]) const
{
  T8_ASSERT (t8_element_is_valid (elem));
  t8_sele_children<eclass_T> ((const t8_standalone_element_t<eclass_T> *) elem,
                              (t8_standalone_element_t<eclass_T> **) c);
#if T8_ENABLE_DEBUG
  for (int i = 0; i < length; i++) {
    T8_ASSERT (t8_element_is_valid (c[i]));
  }
#endif
}

template <t8_eclass_t eclass_T>
int
t8_standalone_scheme_c<eclass_T>::t8_element_child_id (const t8_element_t *elem) const
{
  T8_ASSERT (t8_element_is_valid (elem));
  return t8_sele_child_id<eclass_T> ((const t8_standalone_element_t<eclass_T> *) elem);
}

template <t8_eclass_t eclass_T>
int
t8_standalone_scheme_c<eclass_T>::t8_element_ancestor_id (const t8_element_t *elem, int level) const
{
  T8_ASSERT (t8_element_is_valid (elem));
  T8_ASSERT (0 <= level && level <= T8_ELEMENT_MAXLEVEL[eclass_T]);
  return t8_sele_ancestor_id<eclass_T> ((const t8_standalone_element_t<eclass_T> *) elem, level);
}

template <t8_eclass_t eclass_T>
int
t8_standalone_scheme_c<eclass_T>::t8_element_is_family (t8_element_t *const*fam) const
{
#if T8_ENABLE_DEBUG
  int num_siblings = t8_element_num_siblings (fam[0]);
  for (int i = 0; i < num_siblings; i++) {
    T8_ASSERT (t8_element_is_valid (fam[i]));
  }
#endif
  return t8_sele_is_family<eclass_T> ((t8_standalone_element_t<eclass_T> **) fam);
}

template <t8_eclass_t eclass_T>
void
t8_standalone_scheme_c<eclass_T>::t8_element_nca (const t8_element_t *elem1, const t8_element_t *elem2,
                                                  t8_element_t *nca) const
{
  T8_ASSERT (t8_element_is_valid (elem1));
  T8_ASSERT (t8_element_is_valid (elem2));

  t8_sele_nearest_common_ancestor<eclass_T> ((const t8_standalone_element_t<eclass_T> *) elem1,
                                             (const t8_standalone_element_t<eclass_T> *) elem2,
                                             (t8_standalone_element_t<eclass_T> *) nca);
  T8_ASSERT (t8_element_is_valid (nca));
}

template <t8_eclass_t eclass_T>
t8_element_shape_t
t8_standalone_scheme_c<eclass_T>::t8_element_face_shape (const t8_element_t *elem, int face) const
{
  SC_ABORT ("This function is not implemented yet.\n");
  return T8_ECLASS_ZERO;
}

template <t8_eclass_t eclass_T>
void
t8_standalone_scheme_c<eclass_T>::t8_element_children_at_face (const t8_element_t *elem, int face,
                                                               t8_element_t *children[], int num_children,
                                                               int *child_indices) const
{
  t8_sele_children_at_face<eclass_T> ((const t8_standalone_element_t<eclass_T> *) elem, face,
                                      (t8_standalone_element_t<eclass_T> **) children, num_children, child_indices);
}

template <t8_eclass_t eclass_T>
int
t8_standalone_scheme_c<eclass_T>::t8_element_face_child_face (const t8_element_t *elem, int face, int face_child) const
{
  return t8_sele_face_child_face<eclass_T> ((const t8_standalone_element_t<eclass_T> *) elem, face, face_child);
}

template <t8_eclass_t eclass_T>
int
t8_standalone_scheme_c<eclass_T>::t8_element_face_parent_face (const t8_element_t *elem, int face) const
{
  return t8_sele_face_parent_face<eclass_T> ((const t8_standalone_element_t<eclass_T> *) elem, face);
}

template <t8_eclass_t eclass_T>
int
t8_standalone_scheme_c<eclass_T>::t8_element_tree_face (const t8_element_t *elem, int face) const
{
  return t8_sele_tree_face ((const t8_standalone_element_t<eclass_T> *) elem, face);
}

template <t8_eclass_t eclass_T>
void
t8_standalone_scheme_c<eclass_T>::t8_element_transform_face (const t8_element_t *elem1, t8_element_t *elem2,
                                                             int orientation, int sign, int is_smaller_face) const
{
  t8_sele_transform_face<eclass_T> ((const t8_standalone_element_t<eclass_T> *) elem1,
                                    (t8_standalone_element_t<eclass_T> *) elem2, orientation, sign, is_smaller_face);
}

template <t8_eclass_t eclass_T>
int
t8_standalone_scheme_c<eclass_T>::t8_element_extrude_face (const t8_element_t *face,
                                                           const t8_eclass_scheme_c *face_scheme, t8_element_t *elem,
                                                           int root_face) const
{
  //  T8_ASSERT(t8_eclass_scheme is correct)
  t8_eclass_t face_eclass = T8_ECLASS_ZERO;
  if constexpr (eclass_T == T8_ECLASS_VERTEX)
    SC_ABORT_NOT_REACHED ();
  if constexpr (eclass_T == T8_ECLASS_LINE)
    face_eclass = T8_ECLASS_VERTEX;
  if constexpr (eclass_T == T8_ECLASS_QUAD)
    face_eclass = T8_ECLASS_LINE;
  if constexpr (eclass_T == T8_ECLASS_HEX)
    face_eclass = T8_ECLASS_QUAD;
  if constexpr (T8_ELEMENT_NUM_EQUATIONS[eclass_T]) {
    face_eclass = t8_sele_lut_rootface_to_eclass<eclass_T>[root_face];
  }
  if (face_eclass == T8_ECLASS_VERTEX) {
    return t8_sele_extrude_face<eclass_T, T8_ECLASS_VERTEX> ((const t8_standalone_element_t<T8_ECLASS_VERTEX> *) face,
                                                             (t8_standalone_element_t<eclass_T> *) elem, root_face);
  }
  if (face_eclass == T8_ECLASS_LINE) {
    return t8_sele_extrude_face<eclass_T, T8_ECLASS_LINE> ((const t8_standalone_element_t<T8_ECLASS_LINE> *) face,
                                                           (t8_standalone_element_t<eclass_T> *) elem, root_face);
  }
  if (face_eclass == T8_ECLASS_TRIANGLE) {
    return t8_sele_extrude_face<eclass_T, T8_ECLASS_TRIANGLE> (
      (const t8_standalone_element_t<T8_ECLASS_TRIANGLE> *) face, (t8_standalone_element_t<eclass_T> *) elem,
      root_face);
  }
  if (face_eclass == T8_ECLASS_QUAD) {
    return t8_sele_extrude_face<eclass_T, T8_ECLASS_QUAD> ((const t8_standalone_element_t<T8_ECLASS_QUAD> *) face,
                                                           (t8_standalone_element_t<eclass_T> *) elem, root_face);
  }
  return 0;
}

template <t8_eclass_t eclass_T>
void
t8_standalone_scheme_c<eclass_T>::t8_element_boundary_face (const t8_element_t *elem, int face, t8_element_t *boundary,
                                                            const t8_eclass_scheme_c *boundary_scheme) const
{
  //  T8_ASSERT(t8_eclass_scheme is correct)
  int root_face = t8_element_tree_face (elem, face);
  t8_eclass_t face_eclass = T8_ECLASS_ZERO;
  if constexpr (eclass_T == T8_ECLASS_VERTEX)
    SC_ABORT_NOT_REACHED ();
  if constexpr (eclass_T == T8_ECLASS_LINE)
    face_eclass = T8_ECLASS_VERTEX;
  if constexpr (eclass_T == T8_ECLASS_QUAD)
    face_eclass = T8_ECLASS_LINE;
  if constexpr (eclass_T == T8_ECLASS_HEX)
    face_eclass = T8_ECLASS_QUAD;
  if constexpr (T8_ELEMENT_NUM_EQUATIONS[eclass_T]) {
    face_eclass = t8_sele_lut_rootface_to_eclass<eclass_T>[root_face];
  }
  if (face_eclass == T8_ECLASS_VERTEX) {
    t8_debugf ("call t8_sele_boundary_face with face_eclass = T8_ECLASS_VERTEX\n");
    t8_sele_boundary_face<eclass_T, T8_ECLASS_VERTEX> ((const t8_standalone_element_t<eclass_T> *) elem, root_face,
                                                       (t8_standalone_element_t<T8_ECLASS_VERTEX> *) boundary);
  }
  if (face_eclass == T8_ECLASS_LINE) {
    t8_debugf ("call t8_sele_boundary_face with face_eclass = T8_ECLASS_LINE\n");
    t8_sele_boundary_face<eclass_T, T8_ECLASS_LINE> ((const t8_standalone_element_t<eclass_T> *) elem, root_face,
                                                     (t8_standalone_element_t<T8_ECLASS_LINE> *) boundary);
  }
  if (face_eclass == T8_ECLASS_TRIANGLE) {
    t8_debugf ("call t8_sele_boundary_face with face_eclass = T8_ECLASS_TRIANGLE\n");
    t8_sele_boundary_face<eclass_T, T8_ECLASS_TRIANGLE> ((const t8_standalone_element_t<eclass_T> *) elem, root_face,
                                                         (t8_standalone_element_t<T8_ECLASS_TRIANGLE> *) boundary);
  }
  if (face_eclass == T8_ECLASS_QUAD) {
    t8_debugf ("call t8_sele_boundary_face with face_eclass = T8_ECLASS_QUAD\n");
    t8_sele_boundary_face<eclass_T, T8_ECLASS_QUAD> ((const t8_standalone_element_t<eclass_T> *) elem, root_face,
                                                     (t8_standalone_element_t<T8_ECLASS_QUAD> *) boundary);
  }
}

template <t8_eclass_t eclass_T>
void
t8_standalone_scheme_c<eclass_T>::t8_element_first_descendant_face (const t8_element_t *elem, int face,
                                                                    t8_element_t *first_desc, int level) const
{
  t8_sele_first_descendant_face ((const t8_standalone_element_t<eclass_T> *) elem, face,
                                 (t8_standalone_element_t<eclass_T> *) first_desc, level);
}

template <t8_eclass_t eclass_T>
void
t8_standalone_scheme_c<eclass_T>::t8_element_last_descendant_face (const t8_element_t *elem, int face,
                                                                   t8_element_t *last_desc, int level) const
{
  t8_sele_last_descendant_face ((const t8_standalone_element_t<eclass_T> *) elem, face,
                                (t8_standalone_element_t<eclass_T> *) last_desc, level);
}

template <t8_eclass_t eclass_T>
void
t8_standalone_scheme_c<eclass_T>::t8_element_boundary (const t8_element_t *elem, int min_dim, int length,
                                                       t8_element_t **boundary) const
{
  SC_ABORT ("This function is not implemented yet.\n");
}

template <t8_eclass_t eclass_T>
int
t8_standalone_scheme_c<eclass_T>::t8_element_is_root_boundary (const t8_element_t *elem, int face) const
{
  return t8_sele_is_root_boundary<eclass_T> ((const t8_standalone_element_t<eclass_T> *) elem, face);
}

template <t8_eclass_t eclass_T>
int
t8_standalone_scheme_c<eclass_T>::t8_element_face_neighbor_inside (const t8_element_t *elem, t8_element_t *neigh,
                                                                   int face, int *neigh_face) const
{
  return t8_sele_face_neighbor<eclass_T> ((const t8_standalone_element_t<eclass_T> *) elem,
                                          (t8_standalone_element_t<eclass_T> *) neigh, face, neigh_face);
}

template <t8_eclass_t eclass_T>
t8_element_shape_t
t8_standalone_scheme_c<eclass_T>::t8_element_shape (const t8_element_t *elem) const
{
  T8_ASSERT (t8_element_is_valid (elem));
  return t8_sele_shape<eclass_T> ((const t8_standalone_element_t<eclass_T> *) elem);
}

template <t8_eclass_t eclass_T>
void
t8_standalone_scheme_c<eclass_T>::t8_element_set_linear_id (t8_element_t *elem, int level, t8_linearidx_t id) const
{
  T8_ASSERT (t8_element_is_valid (elem));
  t8_sele_init_linear_id<eclass_T> ((t8_standalone_element_t<eclass_T> *) elem, level, id);
}

template <t8_eclass_t eclass_T>
t8_linearidx_t
t8_standalone_scheme_c<eclass_T>::t8_element_get_linear_id (const t8_element_t *elem, int level) const
{
  T8_ASSERT (t8_element_is_valid (elem));
  return t8_sele_linear_id<eclass_T> ((const t8_standalone_element_t<eclass_T> *) elem, level);
}

template <t8_eclass_t eclass_T>
void
t8_standalone_scheme_c<eclass_T>::t8_element_first_descendant (const t8_element_t *elem, t8_element_t *desc,
                                                               int level) const
{
  T8_ASSERT (t8_element_is_valid (elem));
  t8_sele_first_descendant<eclass_T> ((const t8_standalone_element_t<eclass_T> *) elem,
                                      (t8_standalone_element_t<eclass_T> *) desc, level);
  T8_ASSERT (t8_element_is_valid (desc));
}

template <t8_eclass_t eclass_T>
void
t8_standalone_scheme_c<eclass_T>::t8_element_last_descendant (const t8_element_t *elem, t8_element_t *desc,
                                                              int level) const
{
  T8_ASSERT (t8_element_is_valid (elem));
  t8_sele_last_descendant<eclass_T> ((const t8_standalone_element_t<eclass_T> *) elem,
                                     (t8_standalone_element_t<eclass_T> *) desc, level);
  T8_ASSERT (t8_element_is_valid (desc));
}

template <t8_eclass_t eclass_T>
void
t8_standalone_scheme_c<eclass_T>::t8_element_successor (const t8_element_t *t, t8_element_t *s) const
{
  T8_ASSERT (t8_element_is_valid (t));
  t8_sele_successor<eclass_T> ((const t8_standalone_element_t<eclass_T> *) t, (t8_standalone_element_t<eclass_T> *) s, t8_element_level(t));
  T8_ASSERT (t8_element_is_valid (s));
}

template <t8_eclass_t eclass_T>
int
t8_standalone_scheme_c<eclass_T>::t8_element_root_len (const t8_element_t *elem) const
{
  T8_ASSERT (t8_element_is_valid (elem));
  return t8_sele_get_root_len<eclass_T> ();
}

template <t8_eclass_t eclass_T>
void
t8_standalone_scheme_c<eclass_T>::t8_element_vertex_reference_coords (const t8_element_t *t, const int vertex,
                                                                      double coords[]) const
{
  T8_ASSERT (t8_element_is_valid (t));
  t8_sele_vertex_reference_coords<eclass_T> ((const t8_standalone_element_t<eclass_T> *) t, vertex, coords);
}

template <t8_eclass_t eclass_T>
void
t8_standalone_scheme_c<eclass_T>::t8_element_reference_coords (const t8_element_t *elem, const double *ref_coords,
                                                               const size_t num_coords, double *out_coords) const
{
  SC_ABORT ("This function is not implemented yet.\n");
}

template <t8_eclass_t eclass_T>
t8_gloidx_t
t8_standalone_scheme_c<eclass_T>::t8_element_count_leaves (const t8_element_t *elem, int level) const
{
  return t8_sele_num_descendants_at_leveldiff<eclass_T> ((const t8_standalone_element_t<eclass_T> *) elem,
                                                         level - t8_element_level (elem));
}

template <t8_eclass_t eclass_T>
t8_gloidx_t
t8_standalone_scheme_c<eclass_T>::t8_element_count_leaves_from_root (int level) const
{
  T8_ASSERT (level >= 0);
  if constexpr (eclass_T == T8_ECLASS_PYRAMID) {
    t8_linearidx_t two_to_l = 1LL << level;
    t8_linearidx_t eight_to_l = 1LL << (3 * level);
    return ((eight_to_l << 2) - two_to_l) / 3;
  }
  return 1LL << (level * T8_ELEMENT_DIM[eclass_T]);
}

template <t8_eclass_t eclass_T>
void
t8_standalone_scheme_c<eclass_T>::t8_element_general_function (const t8_element_t *elem, const void *indata,
                                                               void *outdata) const
{
  SC_ABORT ("This function is not implemented yet.\n");
}

#ifdef T8_ENABLE_DEBUG
template <t8_eclass_t eclass_T>
int
t8_standalone_scheme_c<eclass_T>::t8_element_is_valid (const t8_element_t *elem) const
{
  T8_ASSERT (elem != NULL);
  return t8_sele_is_valid<eclass_T> ((const t8_standalone_element_t<eclass_T> *) elem);
}

template <t8_eclass_t eclass_T>
void
t8_standalone_scheme_c<eclass_T>::t8_element_debug_print (const t8_element_t *elem) const
{
  //  T8_ASSERT(t8_element_is_valid(elem));
  t8_sele_debug_print<eclass_T> ((const t8_standalone_element_t<eclass_T> *) elem);
  T8_ASSERT (t8_element_is_valid (elem));
}

 /**
 * \brief Fill a string with readable information about the element
 * 
 * \param[in] elem The element to translate into human-readable information
 * \param[in, out] debug_string The string to fill. 
 */
template <t8_eclass_t eclass_T> void
t8_standalone_scheme_c<eclass_T>::t8_element_to_string (const t8_element_t *elem, char *debug_string, const int string_size) const
{
  SC_ABORT("Not implemented \n");
}

#endif

template <t8_eclass_t eclass_T>
void
t8_standalone_scheme_c<eclass_T>::t8_element_root (t8_element_t *elem) const {
  t8_sele_root((t8_standalone_element_t<eclass_T> *)elem);
}


template <t8_eclass_t eclass_T>
void
t8_standalone_scheme_c<eclass_T>::t8_element_new (int length, t8_element_t **elem) const
{
  /* allocate memory */
  T8_ASSERT (this->ts_context != NULL);
  T8_ASSERT (0 <= length);
  T8_ASSERT (elem != NULL);

  for (int i = 0; i < length; ++i) {
    elem[i] = (t8_element_t *) sc_mempool_alloc ((sc_mempool_t *) this->ts_context);
  }

  /* in debug mode, set sensible default values. */
#ifdef T8_ENABLE_DEBUG
  {
    int i;
    for (i = 0; i < length; i++) {
      t8_element_root (elem[i]);
    }
  }
#endif
}

template <t8_eclass_t eclass_T>
void
t8_standalone_scheme_c<eclass_T>::t8_element_init (int length, t8_element_t *elem) const
{
  t8_standalone_element_t<eclass_T> *sele_array = (t8_standalone_element_t<eclass_T> *)elem;
  for (int i = 0; i < length; i++) {
    t8_sele_root(sele_array + i);
  }
}

template <t8_eclass_t eclass_T>
void
t8_standalone_scheme_c<eclass_T>::t8_element_deinit (int length, t8_element_t *elem) const
{
}



template <t8_eclass_t eclass_T>
void
t8_standalone_scheme_c<eclass_T>::t8_element_destroy (int length, t8_element_t **elem) const
{
  T8_ASSERT (this->ts_context != NULL);
  T8_ASSERT (0 <= length);
  T8_ASSERT (elem != NULL);
  for (int i = 0; i < length; ++i) {
    sc_mempool_free ((sc_mempool_t *) this->ts_context, elem[i]);
  }
}

template <t8_eclass_t eclass_T>
void
t8_standalone_scheme_c<eclass_T>::t8_element_MPI_Pack (t8_element_t **const elements, const unsigned int count, void *send_buffer, int buffer_size,
                       int *position, sc_MPI_Comm comm) const
{
  SC_ABORT("Not implemented\n");
}

template <t8_eclass_t eclass_T>
void
t8_standalone_scheme_c<eclass_T>::t8_element_MPI_Pack_size (const unsigned int count, sc_MPI_Comm comm, int *pack_size) const
{
  SC_ABORT("Not implemented\n");
}

template <t8_eclass_t eclass_T> 
void
t8_standalone_scheme_c<eclass_T>::t8_element_MPI_Unpack (void *recvbuf, const int buffer_size, int *position, t8_element_t **elements,
                         const unsigned int count, sc_MPI_Comm comm) const
{
  SC_ABORT("Not implemented\n");
}