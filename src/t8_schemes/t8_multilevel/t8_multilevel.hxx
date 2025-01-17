/*
  This file is part of t8code.
  t8code is a C library to manage a collection (a forest) of multiple
  connected adaptive space-trees of general element classes in parallel.

  Copyright (C) 2025 the developers

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

/** \file t8_multilevel.hxx
 * This file is the point of entry for our multilevel element implementation.
 * This scheme converts any other scheme into a multievel scheme.
 */

#ifndef T8_MULTILEVEL_HXX
#define T8_MULTILEVEL_HXX

#include <t8_element.h>
#include <t8_schemes/t8_crtp.hxx>

/* Forward declaration of the scheme so we can use it as an argument in the eclass schemes function. */
class t8_scheme;

template <typename TUnderlyingElementType>
struct t8_multilevel_element
{
  int hierarchical_level;
  TUnderlyingElementType linear_element;
};

template <class TUnderlyingEclassScheme, typename TUnderlyingElementType>
class t8_multilevel_scheme: private t8_crtp<TUnderlyingEclassScheme> {

 public:
  using multilevel_element = t8_multilevel_element<TUnderlyingElementType>;
  /* Use constructor with modifiable elem size of base class */
  t8_multilevel_scheme ()
  {
    TUnderlyingEclassScheme::TUnderlyingEclassScheme (sizeof (multilevel_element));
  }

  ~t8_multilevel_scheme ()
  {
    TUnderlyingEclassScheme::~TUnderlyingEclassScheme ();
  }

  /** Return the size of a multilevel element.
   * \return  The size of a multilevel element.
   */
  inline size_t
  get_element_size (void) const
  {
    return this->underlying ().get_element_size ();
  }

  /** Allocate memory for an array of elements and initialize them.
   * \param [in] length   The number of elements to be allocated.
   * \param [in,out] elems On input an array of \b length many unallocated
   *                      element pointers.
   *                      On output all these pointers will point to an allocated
   *                      and initialized element.
   * \note Not every element that is created in t8code will be created by a call
   * to this function. However, if an element is not created using \ref element_new,
   * then it is guaranteed that \ref element_init is called on it.
   * \note In debugging mode, an element that was created with \ref element_new
   * must pass \ref element_is_valid.
   * \note If an element was created by \ref element_new then \ref element_init
   * may not be called for it. Thus, \ref element_new should initialize an element
   * in the same way as a call to \ref element_init would.
   * \see element_init
   * \see element_is_valid
   */
  inline void
  element_new (int length, t8_element_t **elem) const
  {
    this->underlying ().element_new (length, elem);
  }

  /** Initialize an array of allocated elements.
   * \param [in] length   The number of elements to be initialized.
   * \param [in,out] elems On input an array of \b length many allocated
   *                       elements.
   * \param [in] called_new True if the elements in \a elem were created by a call
   *                       to \ref element_new. False if no element in \a elem
   *                       was created in this way. The case that only some elements
   *                       were created by \ref element_new should never occur.
   * \note In debugging mode, an element that was passed to \ref element_init
   * must pass \ref element_is_valid.
   * \note If an element was created by \ref element_new then \ref element_init
   * may not be called for it. Thus, \ref element_new should initialize an element
   * in the same way as a call to \ref element_init would.
   * Thus, if \a called_new is true this function should usually do nothing.
   * \see element_new
   * \see element_is_valid
   */
  inline void
  element_init (int length, t8_element_t *elem) const
  {
    this->underlying ().element_init (length, elem);
  }

  /** Return the refinement level of an element.
   * \param [in] elem    The element whose level should be returned.
   * \return             The level of \b elem.
   */
  inline int
  element_get_level (const t8_element_t *elem) const
  {
    return ((multilevel_element *) elem)->hierarchical_level;
  }

  /** Return the maximum allowed level for any element of a given class.
   * \return                      The maximum allowed level for elements of class \b ts.
   */
  inline int
  get_maxlevel (void) const override
  {
    /* The maxlevel is limited by the size of the linear id datatype. 
       Since need to store all parents in the id as well, we can only one level less. */
    return this->underlying ().t8_element_maxlevel () - 1;
  }

  /** Copy all entries of \b source to \b dest. \b dest must be an existing
   *  element. No memory is allocated by this function.
   * \param [in] source The element whose entries will be copied to \b dest.
   * \param [in,out] dest This element's entries will be overwritten with the
   *                    entries of \b source.
   * \note \a source and \a dest may point to the same element.
   */
  inline void
  element_copy (const t8_element_t *source, t8_element_t *dest) const
  {
    T8_ASSERT (element_is_valid (source));
    T8_ASSERT (element_is_valid (dest));
    ((multilevel_element *) dest)->hierarchical_level = ((multilevel_element *) source)->hierarchical_level;
    this->underlying ().element_copy (&((multilevel_element *) source)->linear_element,
                                      &((multilevel_element *) dest)->linear_element);
  }

  /** Compare two elements.
   * \param [in] elem1  The first element.
   * \param [in] elem2  The second element.
   * \return       negative if elem1 < elem2, zero if elem1 equals elem2
   *               and positive if elem1 > elem2.
   *  If elem2 is a copy of elem1 then the elements are equal.
   */
  inline int
  element_compare (const t8_element_t *elem1, const t8_element_t *elem2) const
  {
    T8_ASSERT (element_is_valid (elem1));
    T8_ASSERT (element_is_valid (elem2));
    const int maxlvl = get_maxlevel ();
    const t8_linearidx_t id1 = element_get_linear_id (elem1, maxlvl);
    const t8_linearidx_t id2 = element_get_linear_id (elem2, maxlvl);
    if (id1 < id2)
      return -1;
    if (id1 > id2)
      return 1;
    else
      return 0;
  }

  /** Check if two elements are equal.
   * \param [in] scheme     Implementation of a class scheme.
   * \param [in] elem1  The first element.
   * \param [in] elem2  The second element.
   * \return            1 if the elements are equal, 0 if they are not equal
   */
  inline int
  element_is_equal (const t8_element_t *elem1, const t8_element_t *elem2) const
  {
    T8_ASSERT (element_is_valid (elem1));
    T8_ASSERT (element_is_valid (elem2));
    if (element_get_level (elem1) != element_get_level (elem2))
      return 0;
    return this->underlying ().element_is_equal (&(static_cast<multilevel_element *> (elem1))->linear_element,
                                                 &(static_cast<multilevel_element *> (elem2))->linear_element);
  }

  /** Compute the parent of a given element \b elem and store it in \b parent.
   *  \b parent needs to be an existing element. No memory is allocated by this function.
   *  \b elem and \b parent can point to the same element, then the entries of
   *  \b elem are overwritten by the ones of its parent.
   * \param [in] elem   The element whose parent will be computed.
   * \param [in,out] parent This element's entries will be overwritten by those
   *                    of \b elem's parent.
   *                    The storage for this element must exist
   *                    and match the element class of the parent.
   *                    For a pyramid, for example, it may be either a
   *                    tetrahedron or a pyramid depending on \b elem's childid.
   */
  inline void
  element_get_parent (const t8_element_t *elem, t8_element_t *parent) const
  {
    T8_ASSERT (element_is_valid (elem));
    T8_ASSERT (element_is_valid (parent));
    const int level = element_get_level (elem);
    multilevel_element *elem_m = static_cast<multilevel_element *> (elem);
    multilevel_element *parent_m = static_cast<multilevel_element *> (parent);
    /* Parent is always one hierarchical level lower. */
    parent_m->hierarchical_level = elem_m->hierarchical_level - 1;
    if (elem_m->hierarchical_level == elem_m->linear_element.level) {
      /* Hierarchical and geometrical level of elem are the same,
         so count the hierarchical one up and compute normal parent. */
      parent_m->hierarchical_level = elem_m->hierarchical_level - 1;
      this->underlying ().element_get_parent (&elem_m->linear_element, &parent_m->linear_element);
    }
    else {
      /* Hierarchical and geometrical level of elem are different,
         so the parent is the same as the element but on a lower hierarchical level. */
      parent_m->hierarchical_level = elem_m->hierarchical_level;
      this->underlying ().element_copy (&elem_m->linear_element, &parent_m->linear_element);
    }
  }

  /** Compute a specific sibling of a given element \b elem and store it in \b sibling.
   *  \b sibling needs to be an existing element. No memory is allocated by this function.
   *  \b elem and \b sibling can point to the same element, then the entries of
   *  \b elem are overwritten by the ones of its \b sibid -th sibling.
   * \param [in] elem   The element whose sibling will be computed.
   * \param [in] sibid  The id of the sibling computed.
   * \param [in,out] sibling This element's entries will be overwritten by those
   *                    of \b elem's sibid-th sibling.
   *                    The storage for this element must exist
   *                    and match the element class of the sibling.
   */
  inline void
  element_get_sibling (const t8_element_t *elem, const int sibid, t8_element_t *sibling) const
  {
    T8_ASSERT (element_is_valid (elem));
    T8_ASSERT (element_is_valid (parent));
    multilevel_element *elem_m = static_cast<multilevel_element *> (elem);
    multilevel_element *sibling_m = static_cast<multilevel_element *> (sibling);
    /* Siblings are always on the same hierarchical level. */
    sibling_m->hierarchical_level = elem_m->hierarchical_level;
    if (sibid == 0) {
      /* The first sibling is the parent drawn to the same hierarchical level. */
      this->underlying ().element_get_parent (&elem_m->linear_element, &sibling_m->linear_element);
    }
    else {
      /* All other siblings are shiftet up one id. */
      this->underlying ().element_get_sibling (&elem_m->linear_element, sibid - 1, &sibling_m->linear_element);
    }
  }

  /** Compute the number of faces of a given element.
   * \param [in] elem The element.
   * \return          The number of faces of \a elem.
   */
  inline int
  element_get_num_faces (const t8_element_t *elem) const
  {
    T8_ASSERT (element_is_valid (elem));
    return this->underlying ().element_get_num_faces (&static_cast<multilevel_element *> (elem)->linear_element);
  }

  /** Compute the maximum number of faces of a given element and all of its
   *  descendants.
   * \param [in] elem The element.
   * \return          The maximum number of faces of \a elem and its descendants.
   */
  inline int
  element_get_max_num_faces (const t8_element_t *elem) const
  {
    T8_ASSERT (element_is_valid (elem));
    return this->underlying ().element_get_max_num_faces (&static_cast<multilevel_element *> (elem)->linear_element);
  }

  /** Return the number of children of an element when it is refined.
   * \param [in] elem   The element whose number of children is returned.
   * \return            The number of children of \a elem if it is to be refined.
   */
  inline int
  element_get_num_children (const t8_element_t *elem) const
  {
    T8_ASSERT (element_is_valid (elem));
    /* Increase the number of children by one so that an element becomes child of itself. */
    return 1 + this->underlying ().element_get_num_children (&static_cast<multilevel_element *> (elem)->linear_element);
  }

  /** Return the number of children of an element's face when the element is refined.
   * \param [in] elem   The element whose face is considered.
   * \param [in] face   A face of \a elem.
   * \return            The number of children of \a face if \a elem is to be refined.
   */
  inline int
  element_get_num_face_children (const t8_element_t *elem, const int face) const
  {
    T8_ASSERT (element_is_valid (elem));
    return this->underlying ().element_get_num_face_children (
      &static_cast<multilevel_element *> (elem)->linear_element);
  }

  /** Return the corner number of an element's face corner.
   * \param [in] element  The element.
   * \param [in] face     A face index for \a element.
   * \param [in] corner   A corner index for the face 0 <= \a corner < num_face_corners.
   * \return              The corner number of the \a corner-th vertex of \a face.
   */
  inline int
  element_get_face_corner (const t8_element_t *elem, const int face, const int corner) const
  {
    T8_ASSERT (element_is_valid (elem));
    return this->underlying ().element_get_face_corner (&static_cast<multilevel_element *> (elem)->linear_element);
  }

  /** Return the face numbers of the faces sharing an element's corner.
   * \param [in] element  The element.
   * \param [in] corner   A corner index for the face.
   * \param [in] face     A face index for \a corner.
   * \return              The face number of the \a face-th face at \a corner.
   */
  inline int
  element_get_corner_face (const t8_element_t *elem, int corner, int face) const
  {
    T8_ASSERT (element_is_valid (elem));
    return this->underlying ().element_get_corner_face (&static_cast<multilevel_element *> (elem)->linear_element);
  }

  /** Construct the child element of a given number.
   * \param [in] elem     This must be a valid element, bigger than maxlevel.
   * \param [in] childid  The number of the child to construct.
   * \param [in,out] child        The storage for this element must exist
   *                              and match the element class of the child.
   *                              On output, a valid element.
   * It is valid to call this function with elem = child.
     */
  inline void
  element_get_child (const t8_element_t *elem, const int childid, t8_element_t *child) const
  {
    T8_ASSERT (element_is_valid (elem));
    T8_ASSERT (element_is_valid (child));
    multilevel_element *elem_m = static_cast<multilevel_element *> (elem);
    multilevel_element *child_m = static_cast<multilevel_element *> (child);
    /* Children are always one hierarchical level higher. */
    child_m->hierarchical_level = elem_m->hierarchical_level + 1;
    if (childid == 0) {
      /* The first child is the element itself. */
      element_copy (elem, child);
    }
    else {
      /* The other children are the normal children shifted by one. */
      this->underlying ().element_get_child (&elem_m->linear_element, childid - 1, &child_m->linear_element);
    }
  }

  /** Construct all children of a given element.
   * \param [in] elem     This must be a valid element, bigger than maxlevel.
   * \param [in] length   The length of the output array \a c must match
   *                      the number of children.
   * \param [in,out] c    The storage for these \a length elements must exist
   *                      and match the element class in the children's ordering.
   *                      On output, all children are valid.
   * It is valid to call this function with elem = c[0].
   * \see element_get_num_children
     */
  inline void
  element_get_children (const t8_element_t *elem, const int length, t8_element_t *children[]) const
  {
    T8_ASSERT (element_is_valid (elem));
    T8_ASSERT (length == element_get_num_children (elem));
    /* Children are always one level higher. */
    const int child_level = element_get_level (elem) + 1;
    multilevel_element *elem_m = static_cast<multilevel_element *> (elem);
    multilevel_element **children_m = static_cast<multilevel_element **> (children);

    /* The first child is the element itself. */
    T8_ASSERT (element_is_valid (*children_m));
    element_copy (elem_m, *children_m);
    *(children_m)->hierarchical_level = child_level;

    /* The rest are the normal children. */
    T8_ASSERT (this->underlying ().element_get_num_children (elem_m->linear_element) == length - 1);
    for (size_t child_id = 0; child_id < length - 1; ++child_id) {
      this->underlying ().element_get_child (elem_m->linear_element, children_m[child_id]);
      children_m[child_id]->hierarchical_level = child_level;
    }
  }

  /** Compute the child id of an element.
   * \param [in] elem     This must be a valid element.
   * \return              The child id of elem.
   */
  inline int
  element_get_child_id (const t8_element_t *elem) const
  {
    T8_ASSERT (element_is_valid (elem));
    multilevel_element *elem_m = static_cast<multilevel_element *> (elem);
    /* If the hierarchical level is one higher than the geometrical level,
    the element is child of itself and has id 0. */
    if (elem_m->hierarchical_level != elem_m->linear_element.level) {
      T8_ASSERT (elem_m->hierarchical_level + 1 == elem_m->linear_element.level);
      return 0;
    }
    /* All other children are shifted by one to make space for the first child. */
    return 1 + this->underlying ().element_get_child_id (&elem_m->linear_element);
  }

  /** Compute the ancestor id of an element, that is the child id
   * at a given level.
   * \param [in] elem     This must be a valid element.
   * \param [in] level    A refinement level. Must satisfy \a level < elem.level
   * \return              The child_id of \a elem in regard to its \a level ancestor.
   */
  inline int
  element_get_ancestor_id (const t8_element_t *elem, int level) const;

  /** Query whether a given set of elements is a family or not.
   * \param [in] fam      An array of as many elements as an element of class
   *                      \b scheme has siblings.
   * \return              Zero if \b fam is not a family, nonzero if it is.
   * \note level 0 elements do not form a family.
   */
  inline int
  elements_are_family (t8_element_t *const *fam) const;

  /** Compute the nearest common ancestor of two elements. That is,
   * the element with highest level that still has both given elements as
   * descendants.
   * \param [in] elem1    The first of the two input elements.
   * \param [in] elem2    The second of the two input elements.
   * \param [in,out] nca  The storage for this element must exist
   *                      and match the element class of the child.
   *                      On output the unique nearest common ancestor of
   *                      \b elem1 and \b elem2.
   */
  inline void
  element_get_nca (const t8_element_t *elem1, const t8_element_t *elem2, t8_element_t *nca) const;

  /** Compute the shape of the face of an element.
   * \param [in] elem     The element.
   * \param [in] face     A face of \a elem.
   * \return              The element shape of the face.
   */
  inline t8_element_shape_t
  element_get_face_shape (const t8_element_t *elem, int face) const;

  /** Given an element and a face of the element, compute all children of
   * the element that touch the face.
   * \param [in] elem     The element.
   * \param [in] face     A face of \a elem.
   * \param [in,out] children Allocated elements, in which the children of \a elem
   *                      that share a face with \a face are stored.
   *                      They will be stored in order of their linear id.
   * \param [in] num_children The number of elements in \a children. Must match
   *                      the number of children that touch \a face.
   *                      \ref element_get_num_face_children
   * \param [in,out] child_indices If not NULL, an array of num_children integers must be given,
   *                      on output its i-th entry is the child_id of the i-th face_child.
   * It is valid to call this function with elem = children[0].
   */
  inline void
  element_get_children_at_face (const t8_element_t *elem, int face, t8_element_t *children[], int num_children,
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
   *                      the order of children from a call to \ref element_get_children_at_face.
   * \return              The face number of the face of a child of \a elem
   *                      that coincides with \a face_child.
   */
  inline int
  element_face_get_child_face (const t8_element_t *elem, int face, int face_child) const;

  /** Given a face of an element return the face number
     * of the parent of the element that matches the element's face. Or return -1 if
     * no face of the parent matches the face.

     * \param [in]  elem    The element.
     * \param [in]  face    Then number of the face.
     * \return              If \a face of \a elem is also a face of \a elem's parent,
     *                      the face number of this face. Otherwise -1.
     * \note For the root element this function always returns \a face.
     */
  inline int
  element_face_get_parent_face (const t8_element_t *elem, int face) const;

  /** Given an element and a face of this element. If the face lies on the
   *  tree boundary, return the face number of the tree face.
   *  If not the return value is arbitrary.
   * \param [in] elem     The element.
   * \param [in] face     The index of a face of \a elem.
   * \return The index of the tree face that \a face is a subface of, if
   *         \a face is on a tree boundary.
   *         Any arbitrary integer if \a is not at a tree boundary.
   */
  inline int
  element_get_tree_face (const t8_element_t *elem, int face) const;

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
  inline void
  element_transform_face (const t8_element_t *elem1, t8_element_t *elem2, int orientation, int sign,
                          int is_smaller_face) const;

  /** Given a boundary face inside a root tree's face construct
   *  the element inside the root tree that has the given face as a
   *  face.
   * \param [in] face     A face element.
   * \param [in,out] elem An allocated element. The entries will be filled with
   *                      the data of the element that has \a face as a face and
   *                      lies within the root tree.
   * \param [in] root_face The index of the face of the root tree in which \a face
   *                      lies.
   * \param [in] scheme   The scheme collection with a scheme for the eclass of the face.
   * \return              The face number of the face of \a elem that coincides
   *                      with \a face.
   */
  inline int
  element_extrude_face (const t8_element_t *face, t8_element_t *elem, int root_face, const t8_scheme *scheme) const;

  /** Construct the first descendant of an element at a given level that touches a given face.
   * \param [in] elem      The input element.
   * \param [in] face      A face of \a elem.
   * \param [in, out] first_desc An allocated element. This element's data will be
   *                       filled with the data of the first descendant of \a elem
   *                       that shares a face with \a face.
   * \param [in] level     The level, at which the first descendant is constructed
   */
  inline void
  element_get_first_descendant_face (const t8_element_t *elem, int face, t8_element_t *first_desc, int level) const;

  /** Construct the last descendant of an element at a given level that touches a given face.
   * \param [in] elem      The input element.
   * \param [in] face      A face of \a elem.
   * \param [in, out] last_desc An allocated element. This element's data will be
   *                       filled with the data of the last descendant of \a elem
   *                       that shares a face with \a face.
   * \param [in] level     The level, at which the last descendant is constructed
   */
  inline void
  element_get_last_descendant_face (const t8_element_t *elem, int face, t8_element_t *last_desc, int level) const;

  /** Construct the boundary element at a specific face.
   * \param [in] elem     The input element.
   * \param [in] face     The index of the face of which to construct the
   *                      boundary element.
   * \param [in,out] boundary An allocated element of dimension of \a element
   *                      minus 1. The entries will be filled with the entries
   *                      of the face of \a element.
   * \param [in] scheme   The scheme containing an eclass scheme for the boundary face.
   */
  inline void
  element_get_boundary_face (const t8_element_t *elem, int face, t8_element_t *boundary, const t8_scheme *scheme) const;

  /** Compute whether a given element shares a given face with its root tree.
   * \param [in] elem     The input element.
   * \param [in] face     A face of \a elem.
   * \return              True if \a face is a subface of the element's root element.
   */
  inline int
  element_is_root_boundary (const t8_element_t *elem, int face) const;

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
  inline int
  element_get_face_neighbor_inside (const t8_element_t *elem, t8_element_t *neigh, int face, int *neigh_face) const;

  /** Initialize the entries of an allocated element according to a
   *  given linear id in a uniform refinement.
   * \param [in,out] elem The element whose entries will be set.
   * \param [in] level    The level of the uniform refinement to consider.
   * \param [in] id       The linear id.
   *                      id must fulfil 0 <= id < 'number of leaves in the uniform refinement'
   */
  inline void
  element_set_linear_id (t8_element_t *elem, int level, t8_linearidx_t id) const;

  /** Compute the linear id of a given element in a hypothetical uniform
   * refinement of a given level.
   * \param [in] elem     The element whose id we compute.
   * \param [in] level    The level of the uniform refinement to consider.
   * \return              The linear id of the element.
   */
  inline t8_linearidx_t
  element_get_linear_id (const t8_element_t *elem, int level) const;

  /** Compute the first descendant of a given element.
   * \param [in] elem     The element whose descendant is computed.
   * \param [out] desc    The first element in a uniform refinement of \a elem
   *                      of the given level.
   * \param [in] level    The level, at which the descendant is computed.
   */
  inline void
  element_get_first_descendant (const t8_element_t *elem, t8_element_t *desc, int level) const;

  /** Compute the last descendant of a given element.
   * \param [in] elem     The element whose descendant is computed.
   * \param [out] desc    The last element in a uniform refinement of \a elem
   *                      of the given level.
   * \param [in] level    The level, at which the descendant is computed.
   */
  inline void
  element_get_last_descendant (const t8_element_t *elem, t8_element_t *desc, int level) const;

  /** Construct the successor in a uniform refinement of a given element.
   * \param [in] elem1    The element whose successor should be constructed.
   * \param [in,out] elem2  The element whose entries will be set.
   * \param [in] level    The level of the uniform refinement to consider.
   */
  inline void
  element_construct_successor (const t8_element_t *elem, t8_element_t *succ) const;

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
  inline void
  element_get_anchor (const t8_element_t *elem, int anchor[3]) const;

  /** Compute the integer coordinates of a given element vertex.
   * The default scheme implements the Morton type SFCs. In these SFCs the
   * elements are positioned in a cube [0,1]^(dL) with dimension d (=0,1,2,3) and 
   * L the maximum refinement level. 
   * All element vertices have integer coordinates in this cube.
   *   \param [in] elem    The element to be considered.
   *   \param [in] vertex  The id of the vertex whose coordinates shall be computed.
   *   \param [out] coords An array of at least as many integers as the element's dimension
   *                      whose entries will be filled with the coordinates of \a vertex.
   */
  inline void
  element_get_vertex_integer_coords (const t8_element_t *elem, int vertex, int coords[]) const;

  /** Compute the coordinates of a given element vertex inside a reference tree
   *  that is embedded into [0,1]^d (d = dimension).
   *   \param [in] elem   The element to be considered.
   *   \param [in] vertex The id of the vertex whose coordinates shall be computed.
   *   \param [out] coords An array of at least as many doubles as the element's dimension
   *                      whose entries will be filled with the coordinates of \a vertex.
   *   \warning           coords should be zero-initialized, as only the first d coords will be set, but when used elsewhere
   *                      all coords might be used. 
   */
  inline void
  element_get_vertex_reference_coords (const t8_element_t *elem, const int vertex, double coords[]) const;

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
  inline void
  element_get_reference_coords (const t8_element_t *elem, const double *ref_coords, const size_t num_coords,
                                double *out_coords) const;

  /** Returns true, if there is one element in the tree, that does not refine into 2^dim children.
   * Returns false otherwise.
   * * \return           0 if the tree is regular, 1 if it is irregular.
   */
  inline int
  refines_irregular (void) const;

#ifdef T8_ENABLE_DEBUG
  /** Query whether a given element can be considered as 'valid' and it is
   *  safe to perform any of the above algorithms on it.
   * \param [in]      elem  The element to be checked.
   * \return          True if \a elem is safe to use. False otherwise.
   * \note            An element that is constructed with \ref element_new
   *                  must pass this test.
   * \note            An element for which \ref element_init was called must pass
   *                  this test.
   * \note            This function is used for debugging to catch certain errors.
   *                  These can for example occur when an element points to a region
   *                  of memory which should not be interpreted as an element.
   * \note            We recommend to use the assertion T8_ASSERT (element_is_valid (elem))
   *                  in the implementation of each of the functions in this file.
   */
  inline int
  element_is_valid (const t8_element_t *t) const;

  /**
  * Print a given element. For a example for a triangle print the coordinates
  * and the level of the triangle. This function is only available in the
  * debugging configuration. 
  * 
  * \param [in]        elem  The element to print
  */
  inline void
  element_to_string (const t8_element_t *elem, char *debug_string, const int string_size) const;
#endif

  /** Fills an element with the root element.
 * \param [in,out] elem   The element to be filled with root.
 */
  inline void
  get_root (t8_element_t *elem) const;

  /** Pack multiple elements into contiguous memory, so they can be sent via MPI.
   * \param [in] elements Array of elements that are to be packed
   * \param [in] count Number of elements to pack
   * \param [in,out] send_buffer Buffer in which to pack the elements
   * \param [in] buffer_size size of the buffer (in order to check that we don't access out of range)
   * \param [in, out] position the position of the first byte that is not already packed
   * \param [in] comm MPI Communicator
  */
  inline void
  element_MPI_Pack (t8_element_t **const elements, const unsigned int count, void *send_buffer, int buffer_size,
                    int *position, sc_MPI_Comm comm) const;

  /** Determine an upper bound for the size of the packed message of \b count elements
   * \param [in] count Number of elements to pack
   * \param [in] comm MPI Communicator
   * \param [out] pack_size upper bound on the message size
  */
  inline void
  element_MPI_Pack_size (const unsigned int count, sc_MPI_Comm comm, int *pack_size) const;

  /** Unpack multiple elements from contiguous memory that was received via MPI.
   * \param [in] recvbuf Buffer from which to unpack the elements
   * \param [in] buffer_size size of the buffer (in order to check that we don't access out of range)
   * \param [in, out] position the position of the first byte that is not already packed
   * \param [in] elements Array of initialised elements that is to be filled from the message
   * \param [in] count Number of elements to unpack
   * \param [in] comm MPI Communicator
  */
  inline void
  element_MPI_Unpack (void *recvbuf, const int buffer_size, int *position, t8_element_t **elements,
                      const unsigned int count, sc_MPI_Comm comm) const;
};

#endif /* !T8_MULTILEVEL_HXX */
