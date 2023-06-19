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

template<t8_eclass_t eclass_T>
t8_standalone_scheme_c<eclass_T>::t8_standalone_scheme_c(void) {
  eclass = eclass_T;
  element_size = sizeof(t8_standalone_element_t<eclass_T>);
  ts_context = sc_mempool_new(element_size);
}

template<t8_eclass_t eclass_T>
t8_standalone_scheme_c<eclass_T>::~t8_standalone_scheme_c() {
  T8_ASSERT (ts_context != NULL);
  sc_mempool_destroy ((sc_mempool_t *) ts_context);
}

template<t8_eclass_t eclass_T>
int
t8_standalone_scheme_c<eclass_T>::t8_element_refines_irregular (void) const {
  if constexpr (eclass_T == T8_ECLASS_PYRAMID) {
     return 1;
  }
  return 0;
}

template<t8_eclass_t eclass_T>
int
t8_standalone_scheme_c<eclass_T>::t8_element_maxlevel (void) const {
  return T8_ELEMENT_MAXLEVEL[eclass_T];
}

template<t8_eclass_t eclass_T>
t8_eclass_t
t8_standalone_scheme_c<eclass_T>::t8_element_child_eclass (int childid) const {
  SC_ABORT("This function is not implemented yet.\n");
  return T8_ECLASS_ZERO; /* suppresses compiler warning */
}

template<t8_eclass_t eclass_T>
int
t8_standalone_scheme_c<eclass_T>::t8_element_level (const t8_element_t *elem) const {
  T8_ASSERT(t8_element_is_valid(elem));
  return ((const t8_standalone_element_t<eclass_T> *)elem)->level;
}

template<t8_eclass_t eclass_T>
void
t8_standalone_scheme_c<eclass_T>::t8_element_copy (const t8_element_t *source,
                                      t8_element_t *dest) const {
  T8_ASSERT(t8_element_is_valid(source));
  t8_sele_copy<eclass_T>((const t8_standalone_element_t<eclass_T> *)source, (t8_standalone_element_t<eclass_T> *)dest);
  T8_ASSERT(t8_element_is_valid(dest));
}

template<t8_eclass_t eclass_T>
int
t8_standalone_scheme_c<eclass_T>::t8_element_compare (const t8_element_t *elem1,
                                        const t8_element_t *elem2) const {
  T8_ASSERT(t8_element_is_valid(elem1));
  T8_ASSERT(t8_element_is_valid(elem2));
  return t8_sele_compare<eclass_T>((const t8_standalone_element_t<eclass_T> *)elem1,
                                  (const t8_standalone_element_t<eclass_T> *)elem2);
}

template<t8_eclass_t eclass_T>
void
t8_standalone_scheme_c<eclass_T>::t8_element_parent (const t8_element_t *elem,
                                        t8_element_t *parent) const {
  T8_ASSERT(t8_element_is_valid(elem));
  t8_sele_parent<eclass_T>((const t8_standalone_element_t<eclass_T> *)elem, (t8_standalone_element_t<eclass_T> *)parent);
  T8_ASSERT(t8_element_is_valid(parent));
}

template<t8_eclass_t eclass_T>
int
t8_standalone_scheme_c<eclass_T>::t8_element_num_siblings (const t8_element_t *elem) const {
  T8_ASSERT(t8_element_is_valid(elem));
  return t8_sele_num_siblings<eclass_T>((const t8_standalone_element_t<eclass_T> *)elem);
}

template<t8_eclass_t eclass_T>
void
t8_standalone_scheme_c<eclass_T>::t8_element_sibling (const t8_element_t *elem,
                                        int sibid,
                                        t8_element_t *sibling) const {
  SC_ABORT("This function is not implemented yet.\n");
}

template<t8_eclass_t eclass_T>
int
t8_standalone_scheme_c<eclass_T>::t8_element_num_corners (const t8_element_t *elem) const {
  T8_ASSERT(t8_element_is_valid(elem));
  return t8_sele_num_corners<eclass_T>((const t8_standalone_element_t<eclass_T> *)elem);
}

template<t8_eclass_t eclass_T>
int
t8_standalone_scheme_c<eclass_T>::t8_element_num_faces (const t8_element_t *elem) const {
  T8_ASSERT(t8_element_is_valid(elem));
  return t8_sele_num_faces<eclass_T>((const t8_standalone_element_t<eclass_T> *) elem);
}

template<t8_eclass_t eclass_T>
int
t8_standalone_scheme_c<eclass_T>::t8_element_max_num_faces (const t8_element_t *elem)
  const {
  T8_ASSERT(t8_element_is_valid(elem));
  return t8_sele_max_num_faces<eclass_T>((const t8_standalone_element_t<eclass_T> *)elem);
}

template<t8_eclass_t eclass_T>
int
t8_standalone_scheme_c<eclass_T>::t8_element_num_children (const t8_element_t *elem) const {
  T8_ASSERT(t8_element_is_valid(elem));
  return t8_sele_num_children<eclass_T>((const t8_standalone_element_t<eclass_T> *)elem);
}

template<t8_eclass_t eclass_T>
int
t8_standalone_scheme_c<eclass_T>::t8_element_num_face_children (const t8_element_t *elem,
                                                  int face) const {
  SC_ABORT("This function is not implemented yet.\n");
  return 0;
}

template<t8_eclass_t eclass_T>
int
t8_standalone_scheme_c<eclass_T>::t8_element_get_face_corner (const t8_element_t *element,
                                                int face,
                                                int corner) const {
  SC_ABORT("This function is not implemented yet.\n");
  return 0;
}

template<t8_eclass_t eclass_T>
int
t8_standalone_scheme_c<eclass_T>::t8_element_get_corner_face (const t8_element_t *element,
                                                int corner,
                                                int face) const {
  SC_ABORT("This function is not implemented yet.\n");
  return 0;
}

template<t8_eclass_t eclass_T>
void
t8_standalone_scheme_c<eclass_T>::t8_element_child (const t8_element_t *elem,
                                      int childid,
                                      t8_element_t *child) const {
  T8_ASSERT(t8_element_is_valid(elem));
  T8_ASSERT(0 <= childid);
  T8_ASSERT(childid < t8_sele_num_children<eclass_T> ((const t8_standalone_element_t<eclass_T> *) elem));
  t8_sele_child<eclass_T>((t8_standalone_element_t<eclass_T> *)elem, childid,
                                                    (t8_standalone_element_t<eclass_T> *)child);
  T8_ASSERT(t8_element_is_valid(child));
}


template<t8_eclass_t eclass_T>
void
t8_standalone_scheme_c<eclass_T>::t8_element_children (const t8_element_t *elem,
                                          int length,
                                          t8_element_t *c[]) const {
  T8_ASSERT(t8_element_is_valid(elem));
  t8_sele_children<eclass_T>((const t8_standalone_element_t<eclass_T> *)elem, (t8_standalone_element_t<eclass_T> **)c);
#if T8_ENABLE_DEBUG
  for (int i = 0; i < length; i++)
  {
    T8_ASSERT(t8_element_is_valid(c[i]));
  }
#endif
}

template<t8_eclass_t eclass_T>
int
t8_standalone_scheme_c<eclass_T>::t8_element_child_id (const t8_element_t *elem) const {
  T8_ASSERT(t8_element_is_valid(elem));
  return t8_sele_child_id<eclass_T>((const t8_standalone_element_t<eclass_T> *) elem);
}

template<t8_eclass_t eclass_T>
int
t8_standalone_scheme_c<eclass_T>::t8_element_ancestor_id (const t8_element_t *elem,
                                            int level) const {
  T8_ASSERT(t8_element_is_valid(elem));
  T8_ASSERT(0 <= level && level <= T8_ELEMENT_MAXLEVEL[eclass_T]);
  return t8_sele_ancestor_id<eclass_T>((const t8_standalone_element_t<eclass_T> *)elem, level);
}

template<t8_eclass_t eclass_T>
int
t8_standalone_scheme_c<eclass_T>::t8_element_is_family (t8_element_t **fam) const {
#if T8_ENABLE_DEBUG
  int num_siblings = t8_element_num_siblings(fam[0]);
  for (int i = 0; i < num_siblings; i++) {
    T8_ASSERT(t8_element_is_valid(fam[i]));
  }
#endif
  return t8_sele_is_family<eclass_T>((t8_standalone_element_t<eclass_T> **)fam);
}

template<t8_eclass_t eclass_T>
void
t8_standalone_scheme_c<eclass_T>::t8_element_nca (const t8_element_t *elem1,
                                    const t8_element_t *elem2,
                                    t8_element_t *nca) const {
  T8_ASSERT(t8_element_is_valid(elem1));
  T8_ASSERT(t8_element_is_valid(elem2));

  t8_sele_nearest_common_ancestor<eclass_T>((const t8_standalone_element_t<eclass_T> *)elem1,
                                                                      (const t8_standalone_element_t<eclass_T> *)elem2,
                                                                      (t8_standalone_element_t<eclass_T> *)nca);
  T8_ASSERT(t8_element_is_valid(nca));
}

template<t8_eclass_t eclass_T>
t8_element_shape_t
t8_standalone_scheme_c<eclass_T>::t8_element_face_shape (const t8_element_t *elem,
                                                  int face) const {
  SC_ABORT("This function is not implemented yet.\n");
  return T8_ECLASS_ZERO;
}

template<t8_eclass_t eclass_T>
void
t8_standalone_scheme_c<eclass_T>::t8_element_children_at_face (const t8_element_t *elem,
                                                  int face,
                                                  t8_element_t *children[],
                                                  int num_children,
                                                  int *child_indices) const {
  SC_ABORT("This function is not implemented yet.\n");
}

template<t8_eclass_t eclass_T>
int
t8_standalone_scheme_c<eclass_T>::t8_element_face_child_face (const t8_element_t *elem,
                                                int face,
                                                int face_child) const {
  SC_ABORT("This function is not implemented yet.\n");
  return 0;
}

template<t8_eclass_t eclass_T>
int
t8_standalone_scheme_c<eclass_T>::t8_element_face_parent_face (const t8_element_t *elem,
                                                  int face) const {
  SC_ABORT("This function is not implemented yet.\n");
  return 0;
}

template<t8_eclass_t eclass_T>
int
t8_standalone_scheme_c<eclass_T>::t8_element_tree_face (const t8_element_t *elem,
                                          int face) const {
  SC_ABORT("This function is not implemented yet.\n");
  return 0;
}

template<t8_eclass_t eclass_T>
void
t8_standalone_scheme_c<eclass_T>::t8_element_transform_face (const t8_element_t *elem1,
                                                t8_element_t *elem2,
                                                int orientation,
                                                int sign,
                                                int is_smaller_face) const {
  SC_ABORT("This function is not implemented yet.\n");
}

template<t8_eclass_t eclass_T>
int
t8_standalone_scheme_c<eclass_T>::t8_element_extrude_face (const t8_element_t *face,
                                              const t8_eclass_scheme_c
                                              *face_scheme,
                                              t8_element_t *elem,
                                              int root_face) const {
  SC_ABORT("This function is not implemented yet.\n");
  return 0;
}

template<t8_eclass_t eclass_T>
void
t8_standalone_scheme_c<eclass_T>::t8_element_boundary_face (const t8_element_t *elem,
                                              int face,
                                              t8_element_t *boundary,
                                              const t8_eclass_scheme_c
                                              *boundary_scheme) const {
  SC_ABORT("This function is not implemented yet.\n");
}

template<t8_eclass_t eclass_T>
void
t8_standalone_scheme_c<eclass_T>::t8_element_first_descendant_face (const t8_element_t
                                                      *elem, int face,
                                                      t8_element_t
                                                      *first_desc,
                                                      int level) const {
  SC_ABORT("This function is not implemented yet.\n");
}

template<t8_eclass_t eclass_T>
void
t8_standalone_scheme_c<eclass_T>::t8_element_last_descendant_face (const t8_element_t
                                                      *elem, int face,
                                                      t8_element_t
                                                      *last_desc,
                                                      int level) const {
  SC_ABORT("This function is not implemented yet.\n");
}

template<t8_eclass_t eclass_T>
void
t8_standalone_scheme_c<eclass_T>::t8_element_boundary (const t8_element_t *elem,
                                          int min_dim, int length,
                                          t8_element_t **boundary) const {
  SC_ABORT("This function is not implemented yet.\n");
}

template<t8_eclass_t eclass_T>
int
t8_standalone_scheme_c<eclass_T>::t8_element_is_root_boundary (const t8_element_t *elem,
                                                  int face) const {
  SC_ABORT("This function is not implemented yet.\n");
  return 0;
}

template<t8_eclass_t eclass_T>
int
t8_standalone_scheme_c<eclass_T>::t8_element_face_neighbor_inside (const t8_element_t
                                                      *elem,
                                                      t8_element_t *neigh,
                                                      int face,
                                                      int *neigh_face) const {
  SC_ABORT("This function is not implemented yet.\n");
  return 0;
}

template<t8_eclass_t eclass_T>
t8_element_shape_t
t8_standalone_scheme_c<eclass_T>::t8_element_shape (const t8_element_t *elem) const {
  T8_ASSERT(t8_element_is_valid(elem));
  return t8_sele_shape<eclass_T>((const t8_standalone_element_t<eclass_T> *)elem);
}

template<t8_eclass_t eclass_T>
void
t8_standalone_scheme_c<eclass_T>::t8_element_set_linear_id (t8_element_t *elem,
                                              int level,
                                              t8_linearidx_t id) const {
  T8_ASSERT(t8_element_is_valid(elem));
  t8_sele_init_linear_id<eclass_T>((t8_standalone_element_t<eclass_T> *)elem, level, id);
}

template<t8_eclass_t eclass_T>
t8_linearidx_t
t8_standalone_scheme_c<eclass_T>::t8_element_get_linear_id (const
                                                  t8_element_t *elem,
                                                  int level) const {
  T8_ASSERT(t8_element_is_valid(elem));
  return t8_sele_linear_id<eclass_T>((const t8_standalone_element_t<eclass_T> *)elem, level);
}

template<t8_eclass_t eclass_T>
void
t8_standalone_scheme_c<eclass_T>::t8_element_first_descendant (const t8_element_t *elem,
                                                  t8_element_t *desc,
                                                  int level) const {
  T8_ASSERT(t8_element_is_valid(elem));
  t8_sele_first_descendant<eclass_T>((const t8_standalone_element_t<eclass_T> *)elem,
                                                                (t8_standalone_element_t<eclass_T> *)desc, level);
  T8_ASSERT(t8_element_is_valid(desc));
}

template<t8_eclass_t eclass_T>
void
t8_standalone_scheme_c<eclass_T>::t8_element_last_descendant (const t8_element_t *elem,
                                                t8_element_t *desc,
                                                int level) const {
  T8_ASSERT(t8_element_is_valid(elem));
  t8_sele_last_descendant<eclass_T>((const t8_standalone_element_t<eclass_T> *)elem,
                                                              (t8_standalone_element_t<eclass_T> *)desc, level);
  T8_ASSERT(t8_element_is_valid(desc));
}

template<t8_eclass_t eclass_T>
void
t8_standalone_scheme_c<eclass_T>::t8_element_successor (const t8_element_t *t,
                                          t8_element_t *s,
                                          int level) const {
  T8_ASSERT(t8_element_is_valid(t));
  t8_sele_successor<eclass_T>((const t8_standalone_element_t<eclass_T> *)t, (t8_standalone_element_t<eclass_T> *)s, level);
  T8_ASSERT(t8_element_is_valid(s));
}

template<t8_eclass_t eclass_T>
int
t8_standalone_scheme_c<eclass_T>::t8_element_root_len (const t8_element_t *elem) const {
  T8_ASSERT(t8_element_is_valid(elem));
  return t8_sele_get_root_len<eclass_T> ();
}

template<t8_eclass_t eclass_T>
void
t8_standalone_scheme_c<eclass_T>::t8_element_vertex_reference_coords (const t8_element_t
                                                        *t,
                                                        const int vertex,
                                                        double coords[])
  const {
  T8_ASSERT(t8_element_is_valid(t));
  t8_sele_vertex_reference_coords<eclass_T>((const t8_standalone_element_t<eclass_T> *) t, vertex, coords);
}

template<t8_eclass_t eclass_T>
void
t8_standalone_scheme_c<eclass_T>::t8_element_reference_coords (const t8_element_t *elem,
                                                  const double *ref_coords,
                                                  const void *user_data,
                                                  double *out_coords)
  const {
  SC_ABORT("This function is not implemented yet.\n");
}

template<t8_eclass_t eclass_T>
t8_gloidx_t
t8_standalone_scheme_c<eclass_T>::t8_element_count_leafs (const t8_element_t *elem, int level) const {
  return t8_sele_num_descendants_at_leveldiff<eclass_T>((const t8_standalone_element_t<eclass_T> *) elem, level - t8_element_level(elem));
}

template<t8_eclass_t eclass_T>
t8_gloidx_t
t8_standalone_scheme_c<eclass_T>::t8_element_count_leafs_from_root (int level) const {
  T8_ASSERT(level >= 0);
  if constexpr(eclass_T == T8_ECLASS_PYRAMID){
    t8_linearidx_t      two_to_l = 1LL << level;
    t8_linearidx_t      eight_to_l = 1LL << (3 * level);
    return ((eight_to_l << 2) - two_to_l) / 3;
  }
  return 1 << (level * T8_ELEMENT_DIM[eclass_T]);
}

template<t8_eclass_t eclass_T>
void
t8_standalone_scheme_c<eclass_T>::t8_element_general_function (const t8_element_t *elem,
                                                  const void *indata,
                                                  void *outdata) const {
  SC_ABORT("This function is not implemented yet.\n");
}

#ifdef T8_ENABLE_DEBUG
template<t8_eclass_t eclass_T>
int
t8_standalone_scheme_c<eclass_T>::t8_element_is_valid (const t8_element_t *elem) const {
  T8_ASSERT(elem != NULL);
  return t8_sele_is_valid<eclass_T>((const t8_standalone_element_t<eclass_T> *) elem);
}

template<t8_eclass_t eclass_T>
void
t8_standalone_scheme_c<eclass_T>::t8_element_debug_print (const t8_element_t *elem) const {
  T8_ASSERT(t8_element_is_valid(elem));
  t8_sele_debug_print<eclass_T>((const t8_standalone_element_t<eclass_T> *)elem);
}
#endif

template<t8_eclass_t eclass_T>
void
t8_standalone_scheme_c<eclass_T>::t8_element_new (int length, t8_element_t **elem) const {
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
    for (i = 0; i < length; i++)
    {
      t8_element_init(1, elem[i], 0);
    }
  }
#endif
}

template<t8_eclass_t eclass_T>
void
t8_standalone_scheme_c<eclass_T>::t8_element_init (int length, t8_element_t *elem, int called_new) const {
#ifdef T8_ENABLE_DEBUG
  if (!called_new)
  {
    int i;
    t8_standalone_element_t<eclass_T> *el = (t8_standalone_element_t<eclass_T> *)elem;
    /* Set all values to 0 */
    for (i = 0; i < length; i++)
    {
      t8_sele_init_linear_id<eclass_T>(el + i, 0, 0);
      T8_ASSERT(t8_sele_is_valid<eclass_T>(el + i));
    }
  }
#endif
}

template<t8_eclass_t eclass_T>
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
