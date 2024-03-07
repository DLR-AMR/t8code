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
  standalone = new t8_standalone_scheme_c<eclass_T>();
  element_size = sizeof (t8_subelement_scheme_c<eclass_T>);
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
  T8_ASSERT (t8_element_is_valid (source));

  t8_element_copy_subelement_values (source, dest);
  standalone->t8_element_copy( (((const t8_element_with_subelements *) source)->elem), (((const t8_element_with_subelements *) dest)->elem) );
}

template <t8_eclass_t eclass_T>
int
t8_subelement_scheme_c<eclass_T>::t8_element_is_subelement (const
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
t8_subelement_scheme_c<eclass_T>::t8_element_reset_subelement_values (t8_element *
                                                                 elem) const
{
  t8_element_with_subelements *sub = (t8_element_with_subelements *) elem;

  sub->transition_type = 0;
  sub->subelement_id = 0;
}

template <t8_eclass_t eclass_T>
void
t8_subelement_scheme_c<eclass_T>::t8_element_get_location_of_subelement (const t8_element_t *elem, int location[]) const
{
  //TODO
}

template <t8_eclass_t eclass_T>
void
t8_subelement_scheme_c<eclass_T>::t8_element_copy_subelement_values (const
                                                                t8_element *
                                                                source,
                                                                t8_element *
                                                                dest) const
{
  const t8_element_with_subelements *sub_source =
    (const t8_element_with_subelements *) source;
  t8_element_with_subelements *sub_dest =
    (t8_element_with_subelements *) dest;

  sub_dest->transition_type = sub_source->transition_type;
  sub_dest->transition_type = sub_source->transition_type;
  sub_dest->subelement_id = sub_source->subelement_id;
}

template <t8_eclass_t eclass_T>
int
t8_subelement_scheme_c<eclass_T>::t8_element_subelement_values_are_valid (const
                                                                 t8_element_t *
                                                                 elem) const
{
  //TODO
}

template <t8_eclass_t eclass_T>
int
t8_subelement_scheme_c<eclass_T>::t8_element_compare (const t8_element_t *elem1, const t8_element_t *elem2) const
{
  const t8_element_with_subelements *sub_elem1 =
    (const t8_element_with_subelements *) elem1;
  const t8_element_with_subelements *sub_elem2 =
    (const t8_element_with_subelements *) elem2;

  const t8_element_t *q = sub_elem1->elem;
  const t8_element_t *r = sub_elem2->elem;

  int elem_comp = standalone->t8_element_compare( q, r );

  if( elem_comp == 0 ) {
    if( t8_element_is_subelement( elem1 ) && t8_element_is_subelement( elem2 ) ) {
      /* Caution: The compare function is used for two subelements. */

      if (sub_elem1->transition_type ==
          sub_elem2->transition_type
          && sub_elem1->subelement_id ==
          sub_elem2->subelement_id) {
        /* both subelements are identical */
        return 0;
      }
    }
    else if( t8_element_is_subelement( elem1 ) ) {
      return -1;                /* elem1 is subelement and therefore smaller */
    }
    else if(t8_element_is_subelement( elem2 )) {
      return 1;                 /* elem2 is subelement and therefore smaller */
    }
  }
  return elem_comp;
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
  return standalone->t8_element_num_siblings( ( (t8_element_with_subelements *) elem )->elem );
}

template <t8_eclass_t eclass_T>
void
t8_subelement_scheme_c<eclass_T>::t8_element_sibling (const t8_element_t *elem, int sibid, t8_element_t *sibling) const
{
   standalone->t8_element_sibling( ( (t8_element_with_subelements *) elem )->elem, sibid, sibling );
}

template <t8_eclass_t eclass_T>
int
t8_subelement_scheme_c<eclass_T>::t8_element_num_corners (const t8_element_t *elem) const
{
   return standalone->t8_element_num_corners( ( (t8_element_with_subelements *) elem )->elem );
}

template <t8_eclass_t eclass_T>
int
t8_subelement_scheme_c<eclass_T>::t8_element_num_faces (const t8_element_t *elem) const
{
   return standalone->t8_element_num_faces( ( (t8_element_with_subelements *) elem )->elem );
}

template <t8_eclass_t eclass_T>
int
t8_subelement_scheme_c<eclass_T>::t8_element_max_num_faces (const t8_element_t *elem) const
{
   return standalone->t8_element_max_num_faces( ( (t8_element_with_subelements *) elem )->elem );
}

template <t8_eclass_t eclass_T>
int
t8_subelement_scheme_c<eclass_T>::t8_element_num_children (const t8_element_t *elem) const
{
   return standalone->t8_element_num_children( ( (t8_element_with_subelements *) elem )->elem );
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
t8_subelement_scheme_c<eclass_T>::t8_element_get_face_corner (const t8_element_t *elem, int face, int corner) const
{
  return standalone->t8_element_get_face_corner( ( (t8_element_with_subelements *) elem )->elem, face, corner );
}

template <t8_eclass_t eclass_T>
int
t8_subelement_scheme_c<eclass_T>::t8_element_get_corner_face (const t8_element_t *elem, int corner, int face) const
{
  return standalone->t8_element_get_face_corner( ( (t8_element_with_subelements *) elem )->elem, corner, face );
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
  standalone->t8_element_children( ( (t8_element_with_subelements *) elem )->elem, length, c );
}

template <t8_eclass_t eclass_T>
int
t8_subelement_scheme_c<eclass_T>::t8_element_child_id (const t8_element_t *elem) const
{
  T8_ASSERT (t8_element_is_valid (elem));
  return (t8_element_is_subelement ( ( (t8_element_with_subelements *) elem )->elem) ? ( (t8_element_with_subelements *) elem )->subelement_id :
        standalone->t8_element_child_id( elem ) );
}

template <t8_eclass_t eclass_T>
int
t8_subelement_scheme_c<eclass_T>::t8_element_ancestor_id (const t8_element_t *elem, int level) const
{
  /* this function is not implemented for subelements */
  T8_ASSERT (!t8_element_is_subelement (elem));

  return standalone->t8_element_ancestor_id( ( (t8_element_with_subelements *) elem )->elem, level );
}

template <t8_eclass_t eclass_T>
int
t8_subelement_scheme_c<eclass_T>::t8_element_is_family (t8_element_t **fam) const
{
#ifdef T8_ENABLE_DEBUG
  {
    int                 i;
    int                 num_siblings = t8_element_num_siblings (fam[0]);
    for (i = 0; i < num_siblings; i++) {
      T8_ASSERT (t8_element_is_valid (fam[i]));
    }
  }
#endif

  /* Subelements can not be refined into other elements of a higher level. 
   * So if the first element of fam is a subelement, we assume that the following num_siblings 
   * many elements are its siblings and therefore form a family. */
  /* Note that this test is very rudimentary, especially when there subelements are in fam */
  t8_element_with_subelements **element_w_sub_family =
    (t8_element_with_subelements **) fam;

  if( element_w_sub_family[0]->transition_type != 0 ) {
    return 1;
  }  /* If the first element of fam is no subelement we check the following elements of fam */
  else {
    /* If any of the following elements is a subelement, then they can not form a family */
    for( int i_sib = 0; i_sib < T8_ELEMENT_NUM_CHILDREN[eclass_T]; i_sib++ )
    {
      if( element_w_sub_family[i_sib]->transition_type != 0 )
        return 0;
    }
  }
  return standalone->t8_element_is_family( fam );
}

template <t8_eclass_t eclass_T>
void
t8_subelement_scheme_c<eclass_T>::t8_element_nca (const t8_element_t *elem1, const t8_element_t *elem2,
                                                  t8_element_t *nca) const
{
 //TODO
}

template <t8_eclass_t eclass_T>
t8_element_shape_t
t8_subelement_scheme_c<eclass_T>::t8_element_face_shape (const t8_element_t *elem, int face) const
{
  return standalone->t8_element_face_shape( elem, face );
}

template <t8_eclass_t eclass_T>
void
t8_subelement_scheme_c<eclass_T>::t8_element_children_at_face (const t8_element_t *elem, int face,
                                                               t8_element_t *children[], int num_children,
                                                               int *child_indices) const
{
  standalone->t8_element_children_at_face( elem, face, children, num_children, child_indices );
}

template <t8_eclass_t eclass_T>
int
t8_subelement_scheme_c<eclass_T>::t8_element_face_child_face (const t8_element_t *elem, int face, int face_child) const
{
  T8_ASSERT (t8_element_is_valid (elem));

  if (t8_element_is_subelement (elem)) {
    T8_ASSERT (face == 1);
    return t8_element_face_parent_face (elem, face);
  }
  else {
    return standalone->t8_element_face_child_face( elem, face, face_child );
  }
}

template <t8_eclass_t eclass_T>
int
t8_subelement_scheme_c<eclass_T>::t8_element_face_parent_face (const t8_element_t *elem, int face) const
{
  if (face == -1) {
    return -1;
  }

  int                 child_id;

  /* For subelements we need to adjust the output of this function.
   * A subelements face is a subface of the parent quadrant (the transition cell) if and only if the face number is 1. */
  if (t8_element_is_subelement (elem)) {
    if (face == 1) {
      /* In this case the face is a subface of the parent. We use the location function in order
       * to determine which of the parents faces intersects the subelements face. */
      int                 location[3] = { };
      t8_element_get_location_of_subelement (elem, location);

      /* subelements in location are enumerated clockwise (not as quadrant faces) */
      return 0;//TODO: subelement_location_to_parent_face[location[0]]; -> for T8_ECLASS
    }
    else {
      return -1;
    }
  }

  return standalone->t8_element_face_parent_face( elem, face );
}

template <t8_eclass_t eclass_T>
int
t8_subelement_scheme_c<eclass_T>::t8_element_tree_face (const t8_element_t *elem, int face) const
{
  T8_ASSERT (t8_element_is_valid (elem));

  /* If elem is a subelement, then this function should only be called together with 
   * face = 1 since other faces will never intersect a tree face. */
  if (t8_element_is_subelement (elem)) {
    T8_ASSERT (face == 1);

    return t8_element_face_parent_face (elem, face);
  }
  else {
    return standalone->t8_element_tree_face( elem, face );
  }
}

template <t8_eclass_t eclass_T>
void
t8_subelement_scheme_c<eclass_T>::t8_element_transform_face (const t8_element_t *elem1, t8_element_t *elem2,
                                                             int orientation, int sign, int is_smaller_face) const
{
  /* this function is not implemented for subelements */
  T8_ASSERT (!t8_element_is_subelement (elem1));

  standalone->t8_element_transform_face( elem1, elem2, orientation, sign, is_smaller_face );
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
  //TODO
}

template <t8_eclass_t eclass_T>
void
t8_subelement_scheme_c<eclass_T>::t8_element_first_descendant_face (const t8_element_t *elem, int face,
                                                                    t8_element_t *first_desc, int level) const
{
  /* this function is not implemented for subelements */
  T8_ASSERT (!t8_element_is_subelement (elem));
  standalone->t8_element_first_descendant_face( elem, face, first_desc, level );
}

template <t8_eclass_t eclass_T>
void
t8_subelement_scheme_c<eclass_T>::t8_element_last_descendant_face (const t8_element_t *elem, int face,
                                                                   t8_element_t *last_desc, int level) const
{
  /* this function is not implemented for subelements */
  T8_ASSERT (!t8_element_is_subelement (elem));
  T8_ASSERT (!t8_element_is_subelement (last_desc));
  standalone->t8_element_last_descendant_face( elem, face, last_desc, level );
}

template <t8_eclass_t eclass_T>
void
t8_subelement_scheme_c<eclass_T>::t8_element_boundary (const t8_element_t *elem, int min_dim, int length,
                                                       t8_element_t **boundary) const
{
  if( t8_element_is_subelement( elem ) ) {
    //TODO
  }
  else {
    standalone->t8_element_boundary( elem, min_dim, length, boundary );
  }
}

template <t8_eclass_t eclass_T>
int
t8_subelement_scheme_c<eclass_T>::t8_element_is_root_boundary (const t8_element_t *elem, int face) const
{
  T8_ASSERT (t8_element_is_valid (elem));

  /* In case of a subelement, we need to change its face number to the face number of the parent quad */
  if (t8_element_is_subelement (elem)) {
    if (face == 1) {
      /* adjust face of subelement to face of parent */
      face = t8_element_face_parent_face (elem, face); //really elem? or: ( (t8_element_with_subelements *) elem )->elem?
    }
    else {                      /* in case of a subelement and face 0 or 2 the face is no subface of the root boundary */
      return 0;
    }
  }
  return standalone->t8_element_is_root_boundary( elem, face );
}

template <t8_eclass_t eclass_T>
int
t8_subelement_scheme_c<eclass_T>::t8_element_face_neighbor_inside (const t8_element_t *elem, t8_element_t *neigh,
                                                                   int face, int *neigh_face) const
{
  //TODO
}

template <t8_eclass_t eclass_T>
t8_element_shape_t
t8_subelement_scheme_c<eclass_T>::t8_element_shape (const t8_element_t *elem) const
{
  return standalone->t8_element_shape( ( (t8_element_with_subelements *) elem )->elem );
}

template <t8_eclass_t eclass_T>
void
t8_subelement_scheme_c<eclass_T>::t8_element_set_linear_id (t8_element_t *elem, int level, t8_linearidx_t id) const
{
  /* this function is not implemented for subelements */
  T8_ASSERT (!t8_element_is_subelement (elem));

  standalone->t8_element_set_linear_id( ( (t8_element_with_subelements *) elem )->elem, level, id );
}

template <t8_eclass_t eclass_T>
t8_linearidx_t
t8_subelement_scheme_c<eclass_T>::t8_element_get_linear_id (const t8_element_t *elem, int level) const
{
  return standalone->t8_element_get_linear_id( ( (t8_element_with_subelements *) elem )->elem, level );
}

template <t8_eclass_t eclass_T>
void
t8_subelement_scheme_c<eclass_T>::t8_element_first_descendant (const t8_element_t *elem, t8_element_t *desc,
                                                               int level) const
{
  standalone->t8_element_first_descendant( ( (t8_element_with_subelements *) elem )->elem, ( (t8_element_with_subelements *) desc )->elem, level );
}

template <t8_eclass_t eclass_T>
void
t8_subelement_scheme_c<eclass_T>::t8_element_last_descendant (const t8_element_t *elem, t8_element_t *desc,
                                                              int level) const
{
  standalone->t8_element_last_descendant( ( (t8_element_with_subelements *) elem )->elem, ( (t8_element_with_subelements *) desc )->elem, level );
}

template <t8_eclass_t eclass_T>
void
t8_subelement_scheme_c<eclass_T>::t8_element_successor (const t8_element_t *t, t8_element_t *s, int level) const
{
  /* this function is not implemented for subelements */
  T8_ASSERT (!t8_element_is_subelement (t));

  standalone->t8_element_successor( ( (t8_element_with_subelements *) t )->elem, ( (t8_element_with_subelements *) s )->elem, level );
}

template <t8_eclass_t eclass_T>
int
t8_subelement_scheme_c<eclass_T>::t8_element_root_len (const t8_element_t *elem) const
{
  return standalone->t8_element_root_len( ( (t8_element_with_subelements *) elem )->elem );
}

template <t8_eclass_t eclass_T>
void
t8_subelement_scheme_c<eclass_T>::t8_element_vertex_reference_coords (const t8_element_t *t, const int vertex,
                                                                      double coords[]) const
{
  standalone->t8_element_vertex_reference_coords( ( (t8_element_with_subelements *) t )->elem, vertex, coords );
}

template <t8_eclass_t eclass_T>
void
t8_subelement_scheme_c<eclass_T>::t8_element_reference_coords (const t8_element_t *elem, const double *ref_coords,
                                                               const size_t num_coords, double *out_coords) const
{
  standalone->t8_element_reference_coords( ( (t8_element_with_subelements *) elem )->elem, ref_coords, num_coords, out_coords );
}

template <t8_eclass_t eclass_T>
t8_gloidx_t
t8_subelement_scheme_c<eclass_T>::t8_element_count_leafs (const t8_element_t *elem, int level) const
{
  /* this function is not implemented for subelements */
  T8_ASSERT (!t8_element_is_subelement (elem));

  return standalone->t8_element_count_leafs( ( (t8_element_with_subelements *) elem )->elem, level );
}

template <t8_eclass_t eclass_T>
t8_gloidx_t
t8_subelement_scheme_c<eclass_T>::t8_element_count_leafs_from_root (int level) const
{
  return standalone->t8_element_count_leafs_from_root( level );
}

template <t8_eclass_t eclass_T>
void
t8_subelement_scheme_c<eclass_T>::t8_element_general_function (const t8_element_t *elem, const void *indata,
                                                               void *outdata) const
{
  /* this function is not implemented for subelements */
  T8_ASSERT (!t8_element_is_subelement (elem));

  standalone->t8_element_general_function( ( (t8_element_with_subelements *) elem )->elem, indata, outdata );
}

#ifdef T8_ENABLE_DEBUG
template <t8_eclass_t eclass_T>
int
t8_subelement_scheme_c<eclass_T>::t8_element_is_valid (const t8_element_t *elem) const
{
  return ( standalone->t8_element_is_valid( ( (t8_element_with_subelements *) elem )->elem ) && t8_element_subelement_values_are_valid (elem) );
}

template <t8_eclass_t eclass_T>
void
t8_subelement_scheme_c<eclass_T>::t8_element_debug_print (const t8_element_t *elem) const
{
  const t8_element_with_subelements *sub =
    (const t8_element_with_subelements *) elem;

  t8_productionf ("\n|------------ t8_element_debug_print subelement: ------------|"
                  "\n|    Transition Type:     %i"
                  "\n|    Subelement ID:       %i"
                  "\n|-------------t8_element_debug_print element: ---------------|\n",
                  sub->transition_type, sub->subelement_id);
  standalone->t8_element_debug_print( ( (t8_element_with_subelements *) elem )->elem );

  /* if the element is not valid, abort, but after printing */
  T8_ASSERT (t8_element_is_valid (elem));
}
#endif

template <t8_eclass_t eclass_T>
void
t8_subelement_scheme_c<eclass_T>::t8_element_new (int length, t8_element_t **elem) const
{
  /* allocate memory */
  T8_ASSERT (this->ts_context != NULL);
  T8_ASSERT (0 <= length);
  T8_ASSERT (elem != NULL);

  for (int i = 0; i < length; ++i) {
    elem[i] = (t8_element_t *) sc_mempool_alloc ((sc_mempool_t *) this->ts_context);
  }

  int                 elem_count;
  for (elem_count = 0; elem_count < length; elem_count++) {
    t8_element_with_subelements *sub =
      (t8_element_with_subelements *) elem[elem_count];
    }
  
  /* in debug mode, set sensible default values. */
#ifdef T8_ENABLE_DEBUG
  {
    for (int i = 0; i < length; i++) {
      t8_element_init (1, elem[i], 0);
    }
  }
#endif
}

template <t8_eclass_t eclass_T>
void
t8_subelement_scheme_c<eclass_T>::t8_element_init (int length, t8_element_t *elem, int called_new) const
{
  t8_element_with_subelements *sub = (t8_element_with_subelements *) elem;

  for (int elem_count = 0; elem_count < length; elem_count++) {
    /* initalize subelement parameters */
    sub[elem_count].subelement_type = 0;
    sub[elem_count].transition_type = 0;
    sub[elem_count].subelement_id = 0;
  }

  standalone->t8_element_init( length, ( (t8_element_with_subelements *) elem )->elem, called_new );
}

template <t8_eclass_t eclass_T>
void
t8_subelement_scheme_c<eclass_T>::t8_element_destroy (int length, t8_element_t **elem) const
{
  //TODO
  for( int i= 0; i<length; i++ )
  {
    standalone->t8_element_destroy( 1, &( (t8_element_with_subelements *) elem )->elem ) ; 
  }
}
