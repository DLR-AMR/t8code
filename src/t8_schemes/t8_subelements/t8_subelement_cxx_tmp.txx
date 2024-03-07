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

#include "t8_subelements_cxx.hxx"

template <t8_eclass_t eclass_T> //TODO
t8_subelement_scheme_c<eclass_T>::t8_subelement_scheme_c (void)
{
  eclass = eclass_T;
  standalone = t8_standalone_scheme_c<eclass_T>();
  element_size = sizeof (t8_subelement_scheme_t<eclass_T>);
  ts_context = sc_mempool_new (element_size);
}

template <t8_eclass_t eclass_T>
t8_subelement_scheme_c<eclass_T>::~t8_subelement_scheme_c ()
{
  T8_ASSERT (ts_context != NULL);
  sc_mempool_destroy ((sc_mempool_t *) ts_context);
}

template <t8_eclass_t eclass_T>
int
t8_subelement_scheme_c<eclass_T>::t8_element_refines_irregular (void) const
{
  /* In general, subelements do not refine regularly */
  return 1;
}

template <t8_eclass_t eclass_T>
int
t8_subelement_scheme_c<eclass_T>::t8_element_maxlevel (void) const
{
  return standalone->t8_element_maxlevel();
}

template <t8_eclass_t eclass_T>
t8_eclass_t
t8_subelement_scheme_c<eclass_T>::t8_element_child_eclass (int childid) const
{
  return standalone->t8_element_child_eclass( childid );
}

template <t8_eclass_t eclass_T>
int
t8_subelement_scheme_c<eclass_T>::t8_element_level (const t8_element_t *elem) const
{
  return standalone->t8_element_level( ((const t8_element_with_subelements *) elem)->elem );
}

template <t8_eclass_t eclass_T>
void
t8_subelement_scheme_c<eclass_T>::t8_element_copy (const t8_element_t *source, t8_element_t *dest) const
{
  //TODO
}

template <t8_eclass_t eclass_T>
int
t8_subelement_scheme_quad_c::t8_element_is_subelement (const
                                                       t8_element *
                                                       elem) const
{
  const t8_element_with_subelements *sub =
    (const t8_element_with_subelements *) elem;

  T8_ASSERT (sub->transition_type >= 0);

  /* transition_type == 0 => elem is no subelement.
   * transition_type != 0 => elem is subelement 
   */
  return (sub->transition_type == 0 ? false : true);
}

template <t8_eclass_t eclass_T>
void
t8_subelement_scheme_quad_c::t8_element_reset_subelement_values (t8_element *
                                                                 elem) const
{
  t8_element_with_subelements *sub = (t8_element_with_subelements *) elem;

  sub->transition_type = 0;
  sub->subelement_id = 0;
}

template <t8_eclass_t eclass_T>
int
t8_subelement_scheme_c<eclass_T>::t8_element_compare (const t8_element_t *elem1, const t8_element_t *elem2) const
{
  T8_ASSERT (t8_element_is_valid (elem1));
  T8_ASSERT (t8_element_is_valid (elem2));
  int compare =  t8_sele_compare<eclass_T> ((const t8_standalone_element_t<eclass_T> *) elem1,
                                    (const t8_standalone_element_t<eclass_T> *) elem2);
  if (compare == 0 && (t8_element_is_subelement (elem1)
                       || t8_element_is_subelement (elem2))) {
    t8_debugf ("Caution, t8_element_compare is used with subelements.\n");
    if (t8_element_is_subelement (elem1)
        && t8_element_is_subelement (elem2)) {
      /* Caution: The compare function is used for two subelements. */

      if (pquad_w_sub_elem1->transition_type ==
          pquad_w_sub_elem2->transition_type
          && pquad_w_sub_elem1->subelement_id ==
          pquad_w_sub_elem2->subelement_id) {
        /* both subelements are identical */
        return 0;
      }
      /* return != 0 to avoid debug abortion in t8_ghost_add_remote */
      return 1;
    }
    else if (t8_element_is_subelement (elem1)) {
      return -1;                /* elem1 is subelement and therefore smaller */
    }
    else if (t8_element_is_subelement (elem2)) {
      return 1;                 /* elem2 is subelement and therefore smaller */
    }
  }
}

template <t8_eclass_t eclass_T>
void
t8_subelement_scheme_c<eclass_T>::t8_element_parent (const t8_element_t *elem, t8_element_t *parent) const
{
  //TODO
}

template <t8_eclass_t eclass_T>
int
t8_subelement_scheme_c<eclass_T>::t8_element_num_siblings (const t8_element_t *elem) const
{
  return t8_sele_num_corners->t8_element_num_siblings( ( (t8_element_with_subelements *) elem )->elem );
}

template <t8_eclass_t eclass_T>
void
t8_subelement_scheme_c<eclass_T>::t8_element_sibling (const t8_element_t *elem, int sibid, t8_element_t *sibling) const
{
   t8_sele_num_corners->t8_element_num_siblings( ( (t8_element_with_subelements *) elem )->elem, sibid, *sibling );
}

template <t8_eclass_t eclass_T>
int
t8_subelement_scheme_c<eclass_T>::t8_element_num_corners (const t8_element_t *elem) const
{
   return t8_sele_num_corners->t8_element_num_corners( ( (t8_element_with_subelements *) elem )->elem );
}

template <t8_eclass_t eclass_T>
int
t8_subelement_scheme_c<eclass_T>::t8_element_num_faces (const t8_element_t *elem) const
{
   return t8_sele_num_corners->t8_element_num_faces( ( (t8_element_with_subelements *) elem )->elem );
}

template <t8_eclass_t eclass_T>
int
t8_subelement_scheme_c<eclass_T>::t8_element_max_num_faces (const t8_element_t *elem) const
{
   return t8_sele_num_corners->t8_element_num_faces( ( (t8_element_with_subelements *) elem )->elem );
}

template <t8_eclass_t eclass_T>
int
t8_subelement_scheme_c<eclass_T>::t8_element_num_children (const t8_element_t *elem) const
{
   return t8_sele_num_corners->t8_element_num_children( ( (t8_element_with_subelements *) elem )->elem );
}

template <t8_eclass_t eclass_T>
int
t8_subelement_scheme_c<eclass_T>::t8_element_num_face_children (const t8_element_t *elem, int face) const
{
  /* this function is not implemented for subelements */
  T8_ASSERT (!t8_element_is_subelement (elem));
  T8_ASSERT (t8_element_is_valid (elem));

  /* if we use this scheme without set_transition, then we are only balanced and two neighbors are possible */
  return 2;
}

template <t8_eclass_t eclass_T>
int
t8_subelement_scheme_c<eclass_T>::t8_element_get_face_corner (const t8_element_t *element, int face, int corner) const
{
  return t8_sele_num_corners->t8_element_get_face_corner( ( (t8_element_with_subelements *) elem )->elem, face, corner );
}

template <t8_eclass_t eclass_T>
int
t8_subelement_scheme_c<eclass_T>::t8_element_get_corner_face (const t8_element_t *element, int corner, int face) const
{
  return t8_sele_num_corners->t8_element_get_face_corner( ( (t8_element_with_subelements *) elem )->elem, corner, face );
}

template <t8_eclass_t eclass_T>
void
t8_subelement_scheme_c<eclass_T>::t8_element_child (const t8_element_t *elem, int childid, t8_element_t *child) const
{
  //TODO
}

template <t8_eclass_t eclass_T>
void
t8_subelement_scheme_c<eclass_T>::t8_element_children (const t8_element_t *elem, int length, t8_element_t *c[]) const
{
  /* if elem is a subelement, then this function will construct the children of its parent p4est quadrant */
  standalone->t8_element_child_eclass( ( (t8_element_with_subelements *) elem )->elem, length, c );
}

template <t8_eclass_t eclass_T>
int
t8_subelement_scheme_c<eclass_T>::t8_element_child_id (const t8_element_t *elem) const
{
  T8_ASSERT (t8_element_is_valid (elem));
  return (t8_element_is_subelement ( ( (t8_element_with_subelements *) elem )->elem) ? pquad_w_sub->subelement_id :
        standalone->t8_element_child_id( elem ) );
}








template <t8_eclass_t eclass_T>
int
t8_subelement_scheme_c<eclass_T>::t8_element_ancestor_id (const t8_element_t *elem, int level) const
{

}








template <t8_eclass_t eclass_T>
int
t8_subelement_scheme_c<eclass_T>::t8_element_is_family (t8_element_t **fam) const
{

}

template <t8_eclass_t eclass_T>
void
t8_subelement_scheme_c<eclass_T>::t8_element_nca (const t8_element_t *elem1, const t8_element_t *elem2,
                                                  t8_element_t *nca) const
{
 
}

template <t8_eclass_t eclass_T>
t8_element_shape_t
t8_subelement_scheme_c<eclass_T>::t8_element_face_shape (const t8_element_t *elem, int face) const
{

}

template <t8_eclass_t eclass_T>
void
t8_subelement_scheme_c<eclass_T>::t8_element_children_at_face (const t8_element_t *elem, int face,
                                                               t8_element_t *children[], int num_children,
                                                               int *child_indices) const
{

}

template <t8_eclass_t eclass_T>
int
t8_subelement_scheme_c<eclass_T>::t8_element_face_child_face (const t8_element_t *elem, int face, int face_child) const
{

}

template <t8_eclass_t eclass_T>
int
t8_subelement_scheme_c<eclass_T>::t8_element_face_parent_face (const t8_element_t *elem, int face) const
{

}

template <t8_eclass_t eclass_T>
int
t8_subelement_scheme_c<eclass_T>::t8_element_tree_face (const t8_element_t *elem, int face) const
{

}

template <t8_eclass_t eclass_T>
void
t8_subelement_scheme_c<eclass_T>::t8_element_transform_face (const t8_element_t *elem1, t8_element_t *elem2,
                                                             int orientation, int sign, int is_smaller_face) const
{

}

template <t8_eclass_t eclass_T>
int
t8_subelement_scheme_c<eclass_T>::t8_element_extrude_face (const t8_element_t *face,
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
t8_subelement_scheme_c<eclass_T>::t8_element_boundary_face (const t8_element_t *elem, int face, t8_element_t *boundary,
                                                            const t8_eclass_scheme_c *boundary_scheme) const
{

}

template <t8_eclass_t eclass_T>
void
t8_subelement_scheme_c<eclass_T>::t8_element_first_descendant_face (const t8_element_t *elem, int face,
                                                                    t8_element_t *first_desc, int level) const
{

}

template <t8_eclass_t eclass_T>
void
t8_subelement_scheme_c<eclass_T>::t8_element_last_descendant_face (const t8_element_t *elem, int face,
                                                                   t8_element_t *last_desc, int level) const
{

}

template <t8_eclass_t eclass_T>
void
t8_subelement_scheme_c<eclass_T>::t8_element_boundary (const t8_element_t *elem, int min_dim, int length,
                                                       t8_element_t **boundary) const
{

}

template <t8_eclass_t eclass_T>
int
t8_subelement_scheme_c<eclass_T>::t8_element_is_root_boundary (const t8_element_t *elem, int face) const
{

}

template <t8_eclass_t eclass_T>
int
t8_subelement_scheme_c<eclass_T>::t8_element_face_neighbor_inside (const t8_element_t *elem, t8_element_t *neigh,
                                                                   int face, int *neigh_face) const
{

}

template <t8_eclass_t eclass_T>
t8_element_shape_t
t8_subelement_scheme_c<eclass_T>::t8_element_shape (const t8_element_t *elem) const
{

}

template <t8_eclass_t eclass_T>
void
t8_subelement_scheme_c<eclass_T>::t8_element_set_linear_id (t8_element_t *elem, int level, t8_linearidx_t id) const
{

}

template <t8_eclass_t eclass_T>
t8_linearidx_t
t8_subelement_scheme_c<eclass_T>::t8_element_get_linear_id (const t8_element_t *elem, int level) const
{

}

template <t8_eclass_t eclass_T>
void
t8_subelement_scheme_c<eclass_T>::t8_element_first_descendant (const t8_element_t *elem, t8_element_t *desc,
                                                               int level) const
{

}

template <t8_eclass_t eclass_T>
void
t8_subelement_scheme_c<eclass_T>::t8_element_last_descendant (const t8_element_t *elem, t8_element_t *desc,
                                                              int level) const
{

}

template <t8_eclass_t eclass_T>
void
t8_subelement_scheme_c<eclass_T>::t8_element_successor (const t8_element_t *t, t8_element_t *s, int level) const
{

}

template <t8_eclass_t eclass_T>
int
t8_subelement_scheme_c<eclass_T>::t8_element_root_len (const t8_element_t *elem) const
{

}

template <t8_eclass_t eclass_T>
void
t8_subelement_scheme_c<eclass_T>::t8_element_vertex_reference_coords (const t8_element_t *t, const int vertex,
                                                                      double coords[]) const
{

}

template <t8_eclass_t eclass_T>
void
t8_subelement_scheme_c<eclass_T>::t8_element_reference_coords (const t8_element_t *elem, const double *ref_coords,
                                                               const size_t num_coords, double *out_coords) const
{

}

template <t8_eclass_t eclass_T>
t8_gloidx_t
t8_subelement_scheme_c<eclass_T>::t8_element_count_leafs (const t8_element_t *elem, int level) const
{

}

template <t8_eclass_t eclass_T>
t8_gloidx_t
t8_subelement_scheme_c<eclass_T>::t8_element_count_leafs_from_root (int level) const
{

}

template <t8_eclass_t eclass_T>
void
t8_subelement_scheme_c<eclass_T>::t8_element_general_function (const t8_element_t *elem, const void *indata,
                                                               void *outdata) const
{

}

#ifdef T8_ENABLE_DEBUG
template <t8_eclass_t eclass_T>
int
t8_subelement_scheme_c<eclass_T>::t8_element_is_valid (const t8_element_t *elem) const
{

}

template <t8_eclass_t eclass_T>
void
t8_subelement_scheme_c<eclass_T>::t8_element_debug_print (const t8_element_t *elem) const
{

}
#endif

template <t8_eclass_t eclass_T>
void
t8_subelement_scheme_c<eclass_T>::t8_element_new (int length, t8_element_t **elem) const
{

}

template <t8_eclass_t eclass_T>
void
t8_subelement_scheme_c<eclass_T>::t8_element_init (int length, t8_element_t *elem, int called_new) const
{

}

template <t8_eclass_t eclass_T>
void
t8_subelement_scheme_c<eclass_T>::t8_element_destroy (int length, t8_element_t **elem) const
{

}
