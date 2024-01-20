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

/** \file t8_default_common_cxx.hxx
 * We provide some functions that are useful across element classes.
 */

#ifndef T8_CONSECUTIVE_COMMON_CXX_HXX
#define T8_CONSECUTIVE_COMMON_CXX_HXX

#include <t8_element_cxx.hxx>

/* Macro to check whether a pointer (VAR) to a base class, comes from an
 * implementation of a child class (TYPE). */
#define T8_COMMON_IS_TYPE(VAR, TYPE) ((dynamic_cast<TYPE> (VAR)) != NULL)

class t8_consecutive_scheme_common_c: public t8_eclass_scheme_c {
 public:
  /** Destructor for all default schemes */
  virtual ~t8_consecutive_scheme_common_c ();

  /** Compute the number of corners of a given element. */
  virtual int
  t8_element_num_corners (const t8_element_t *elem) const;

  /** Allocate space for a bunch of elements. */
  virtual void
  t8_element_new (int length, t8_element_t **elem) const;

  /** Deallocate space for a bunch of elements. */
  virtual void
  t8_element_destroy (int length, t8_element_t **elem) const;

  /** Return the shape of an element */
  virtual t8_element_shape_t
  t8_element_shape (const t8_element_t *elem) const;

  /** Count how many leaf descendants of a given uniform level an element would produce.
   * \param [in] t     The element to be checked.
   * \param [in] level A refinement level.
   * \return Suppose \a t is uniformly refined up to level \a level. The return value
   * is the resulting number of elements (of the given level).
   * Each default element (except pyramids) refines into 2^{dim * (level - level(t))}
   * children.
   */
  virtual t8_gloidx_t
  t8_element_count_leafs (const t8_element_t *t, int level) const;

  /** Compute the number of siblings of an element. That is the number of 
   * Children of its parent.
   * \param [in] elem The element.
   * \return          The number of siblings of \a element.
   * Note that this number is >= 1, since we count the element itself as a sibling.
   */
  virtual int
  t8_element_num_siblings (const t8_element_t *elem) const;

  /** Count how many leaf descendants of a given uniform level the root element will produce.
   * \param [in] level A refinement level.
   * \return The value of \ref t8_element_count_leafs if the input element
   *      is the root (level 0) element.
   */
  virtual t8_gloidx_t
  t8_element_count_leafs_from_root (int level) const;

  /** The common implementation of the general function for the default scheme
   * has no effect. This function literally does nothing.
   * The tri, tet and prism scheme override this implementation with a function that
   * stores the type of the element in \a outdata.
   *  \param [in] elem A valid element
   *  \param [in] indata Is ignored. Can be NULL.
   *  \param [out] outdata Is ignored. Can be NULL.
   * \note Calling this function has no effect. See the specialized implementations in
   * t8_default_tri_cxx.hxx, t8_default_tet_cxx.hxx and t8_default_prism_cxx.hxx.
   */
  virtual void
  t8_element_general_function (const t8_element_t *elem, const void *indata, void *outdata) const;

  /** Convert points in the reference space of an element to points in the
   *  reference space of the tree.
   * 
   * \param [in] elem         The element.
   * \param [in] coords_input The coordinates \f$ [0,1]^\mathrm{dim} \f$ of the point
   *                          in the reference space of the element.
   * \param [in] num_coords   Number of \f$ dim\f$-sized coordinates to evaluate.
   * \param [out] out_coords  The coordinates of the points in the
   *                          reference space of the tree.
   */
  virtual void
  t8_element_reference_coords (const t8_element_t *elem, const double *ref_coords, const size_t num_coords,
                               double *out_coords) const
    = 0;

  /** Get the integer coordinates of the anchor node of an element.
   * The default scheme implements the Morton type SFCs. In these SFCs the
   * elements are positioned in a cube [0,1]^(dL) with dimension d (=0,1,2,3) and 
   * L the maximum refinement level. 
   * All element vertices have integer coordinates in this cube and the anchor
   * node is the first of all vertices (index 0). It also has the lowest x,y and z
   * coordinates.
   * \param [in] elem   The element.
   * \param [out] anchor The integer coordinates of the anchor node in the cube [0,1]^(dL)
   */
  virtual void
  t8_element_anchor (const t8_element_t *elem, int anchor[3]) const
    = 0;
#if T8_ENABLE_DEBUG
  virtual void
  t8_element_debug_print (const t8_element_t *elem) const;
#endif
  virtual int
  t8_element_is_family (t8_element_t **fam) const;

  t8_linearidx_t
  t8_hilbert_num_descendants_of_child_at_leveldiff (t8_element_t *elem, const int childid, const int leveldiff);

  virtual t8_linearidx_t
  t8_element_linear_id_recursive (t8_element_t *elem, const t8_linearidx_t id, const int level_diff) const;

  virtual void
  t8_element_init_linear_id_recursive (t8_element_t *elem, const int level_diff, t8_linearidx_t id) const;

  virtual void
  t8_element_set_linear_id (t8_element_t *elem, int level, t8_linearidx_t id) const;

  virtual t8_linearidx_t
  t8_element_get_linear_id (const t8_element_t *elem, int level) const;

  virtual void
  t8_element_first_descendant (const t8_element_t *elem, t8_element_t *desc, int level) const;

  virtual void
  t8_element_last_descendant (const t8_element_t *elem, t8_element_t *desc, int level) const;

  virtual void
  t8_element_successor (const t8_element_t *elem1, t8_element_t *elem2, int level) const;

  virtual int
  t8_element_ancestor_id (const t8_element_t *elem, int level) const;

  virtual void
  t8_element_nca (const t8_element_t *elem1, const t8_element_t *elem2, t8_element_t *nca) const;

  virtual void
  t8_element_children (const t8_element_t *elem, int length, t8_element_t *children[]) const;

  virtual void
  t8_element_sibling (const t8_element_t *elem, int sibid, t8_element_t *sibling) const;

  virtual int
  t8_element_compare (const t8_element_t *elem1, const t8_element_t *elem2) const;

  virtual void
  t8_element_root (t8_element_t *elem) const;

  /*******************************************************************
 ************************ facestuff *********************************
 *******************************************************************/
  virtual int
  t8_element_num_faces (const t8_element_t *elem) const;

  /** Compute the maximum number of faces of a given element and all of its
   *  descendants.
   * \param [in] elem The element.
   * \return          The maximum number of faces of \a elem and its descendants.
   */
  virtual int
  t8_element_max_num_faces (const t8_element_t *elem) const;

  /** Return the number of children of an element when it is refined.
   * \param [in] elem   The element whose number of children is returned.
   * \return            The number of children of \a elem if it is to be refined.
   */

  virtual int
  t8_element_num_face_children (const t8_element_t *elem, int face) const;
  /** Return the corner number of an element's face corner.
   * \param [in] element  The element.
   * \param [in] face     A face index for \a element.
   * \param [in] corner   A corner index for the face 0 <= \a corner < num_face_corners.
   * \return              The corner number of the \a corner-th vertex of \a face.
   */
  virtual int
  t8_element_get_face_corner (const t8_element_t *element, int face, int corner) const;

  /** Return the face numbers of the faces sharing an element's corner.
   * \param [in] element  The element.
   * \param [in] corner   A corner index for the face.
   * \param [in] face     A face index for \a corner.
   * \return              The face number of the \a face-th face at \a corner.
   */
  virtual int
  t8_element_get_corner_face (const t8_element_t *element, int corner, int face) const;

  /** Return the type of each child in the ordering of the implementation.
   * \param [in] childid  Must be between 0 and the number of children (exclusive).
   *                      The number of children is defined in \a t8_element_num_children.
   * \return              The type for the given child.
   */

  virtual t8_element_shape_t
  t8_element_face_shape (const t8_element_t *elem, int face) const;

  /** Given an element and a face of the element, compute all children of
   * the element that touch the face.
   * \param [in] elem     The element.
   * \param [in] face     A face of \a elem.
   * \param [in,out] children Allocated elements, in which the children of \a elem
   *                      that share a face with \a face are stored.
   *                      They will be stored in order of their linear id.
   * \param [in] num_children The number of elements in \a children. Must match
   *                      the number of children that touch \a face.
   *                      \ref t8_element_num_face_children
   * \param [in,out] child_indices If not NULL, an array of num_children integers must be given,
   *                      on output its i-th entry is the child_id of the i-th face_child.
   * It is valid to call this function with elem = children[0].
   */
  virtual void
  t8_element_children_at_face (const t8_element_t *elem, int face, t8_element_t *children[], int num_children,
                               int *child_indices) const;

  /** Given a face of an element and a child number of a child of that face, return the face number
   * of the child of the element that matches the child face.
   * \verbatim
      x ---- x   x      x           x ---- x
      |      |   |      |           |   |  | <-- f
      |      |   |      x           |   x--x
      |      |   |                  |      |
      x ---- x   x                  x ---- x
       elem    face  face_child    Returns the face number f
     \endverbatim

   * \param [in]  elem    The element.
   * \param [in]  face    Then number of the face.
   * \param [in]  face_child A number 0 <= \a face_child < num_face_children,
   *                      specifying a child of \a elem that shares a face with \a face.
   *                      These children are counted in linear order. This coincides with
   *                      the order of children from a call to \ref t8_element_children_at_face.
   * \return              The face number of the face of a child of \a elem
   *                      that coincides with \a face_child.
   */
  virtual int
  t8_element_face_child_face (const t8_element_t *elem, int face, int face_child) const;

  /** Given a face of an element return the face number
     * of the parent of the element that matches the element's face. Or return -1 if
     * no face of the parent matches the face.

     * \param [in]  elem    The element.
     * \param [in]  face    Then number of the face.
     * \return              If \a face of \a elem is also a face of \a elem's parent,
     *                      the face number of this face. Otherwise -1.
     * \note For the root element this function always returns \a face.
     */
  virtual int
  t8_element_face_parent_face (const t8_element_t *elem, int face) const;

  /** Given an element and a face of this element. If the face lies on the
   *  tree boundary, return the face number of the tree face.
   *  If not the return value is arbitrary.
   * \param [in] elem     The element.
   * \param [in] face     The index of a face of \a elem.
   * \return The index of the tree face that \a face is a subface of, if
   *         \a face is on a tree boundary.
   *         Any arbitrary integer if \a is not at a tree boundary.
   */
  virtual int
  t8_element_tree_face (const t8_element_t *elem, int face) const;

  /** Suppose we have two trees that share a common face f.
   *  Given an element e that is a subface of f in one of the trees
   *  and given the orientation of the tree connection, construct the face
   *  element of the respective tree neighbor that logically coincides with e
   *  but lies in the coordinate system of the neighbor tree.
   *  \param [in] elem1     The face element.
   *  \param [in,out] elem2 On return the face element \a elem1 with respect
   *                        to the coordinate system of the other tree.
   *  \param [in] orientation The orientation of the tree-tree connection.
   *                        \see t8_cmesh_set_join
   *  \param [in] sign      Depending on the topological orientation of the two tree faces,
   *                        either 0 (both faces have opposite orientation)
   *                        or 1 (both faces have the same top. orientattion).
   *                        \ref t8_eclass_face_orientation
   *  \param [in] is_smaller_face Flag to declare whether \a elem1 belongs to
   *                        the smaller face. A face f of tree T is smaller than
   *                        f' of T' if either the eclass of T is smaller or if
   *                        the classes are equal and f<f'. The orientation is
   *                        defined in relation to the smaller face.
   * \note \a elem1 and \a elem2 may point to the same element.
   */
  virtual void
  t8_element_transform_face (const t8_element_t *elem1, t8_element_t *elem2, int orientation, int sign,
                             int is_smaller_face) const;

  /** Given a boundary face inside a root tree's face construct
   *  the element inside the root tree that has the given face as a
   *  face.
   * \param [in] face     A face element.
   * \param [in] face_scheme The scheme for the face element.
   * \param [in,out] elem An allocated element. The entries will be filled with
   *                      the data of the element that has \a face as a face and
   *                      lies within the root tree.
   * \param [in] root_face The index of the face of the root tree in which \a face
   *                      lies.
   * \return              The face number of the face of \a elem that coincides
   *                      with \a face.
   */
  virtual int
  t8_element_extrude_face (const t8_element_t *face, const t8_eclass_scheme_c *face_scheme, t8_element_t *elem,
                           int root_face) const;

  /** Construct the first descendant of an element at a given level that touches a given face.
   * \param [in] elem      The input element.
   * \param [in] face      A face of \a elem.
   * \param [in, out] first_desc An allocated element. This element's data will be
   *                       filled with the data of the first descendant of \a elem
   *                       that shares a face with \a face.
   * \param [in] level     The level, at which the first descendant is constructed
   */
  virtual void
  t8_element_first_descendant_face (const t8_element_t *elem, int face, t8_element_t *first_desc, int level) const;

  /** Construct the last descendant of an element at a given level that touches a given face.
   * \param [in] elem      The input element.
   * \param [in] face      A face of \a elem.
   * \param [in, out] last_desc An allocated element. This element's data will be
   *                       filled with the data of the last descendant of \a elem
   *                       that shares a face with \a face.
   * \param [in] level     The level, at which the last descendant is constructed
   */
  virtual void
  t8_element_last_descendant_face (const t8_element_t *elem, int face, t8_element_t *last_desc, int level) const;

  /** Construct the boundary element at a specific face.
   * \param [in] elem     The input element.
   * \param [in] face     The index of the face of which to construct the
   *                      boundary element.
   * \param [in,out] boundary An allocated element of dimension of \a element
   *                      minus 1. The entries will be filled with the entries
   *                      of the face of \a element.
   * \param [in] boundary_scheme The scheme for the eclass of the boundary face.
   */
  virtual void
  t8_element_boundary_face (const t8_element_t *elem, int face, t8_element_t *boundary,
                            const t8_eclass_scheme_c *boundary_scheme) const;

  /** Compute whether a given element shares a given face with its root tree.
   * \param [in] elem     The input element.
   * \param [in] face     A face of \a elem.
   * \return              True if \a face is a subface of the element's root element.
   */
  virtual int
  t8_element_is_root_boundary (const t8_element_t *elem, int face) const;

  /** Construct the face neighbor of a given element if this face neighbor
   * is inside the root tree. Return 0 otherwise.
   * \param [in] elem The element to be considered.
   * \param [in,out] neigh If the face neighbor of \a elem along \a face is inside
   *                  the root tree, this element's data is filled with the
   *                  data of the face neighbor. Otherwise the data can be modified
   *                  arbitrarily.
   * \param [in] face The number of the face along which the neighbor should be
   *                  constructed.
   * \param [out] neigh_face The number of \a face as viewed from \a neigh.
   *                  An arbitrary value, if the neighbor is not inside the root tree.
   * \return          True if \a neigh is inside the root tree.
   *                  False if not. In this case \a neigh's data can be arbitrary
   *                  on output.
   */
  virtual int
  t8_element_face_neighbor_inside (const t8_element_t *elem, t8_element_t *neigh, int face, int *neigh_face) const;
};
#endif /* !T8_CONSECUTIVE_COMMON_CXX_HXX */
