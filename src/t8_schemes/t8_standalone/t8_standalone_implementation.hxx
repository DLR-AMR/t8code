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

#ifndef T8_STANDALONE_IMPLEMENTATION_HXX
#define T8_STANDALONE_IMPLEMENTATION_HXX

#include <t8_schemes/t8_scheme.hxx>
#include <t8_eclass.h>
#include <sc_functions.h>
#include <t8_schemes/t8_standalone/t8_standalone_elements.hxx>

template <t8_eclass TEclass>
struct t8_standalone_scheme
{
 public:
  /** Constructor
   * \param [in] elem_size  The size of the elements this scheme holds.
  */
  t8_standalone_scheme ()
    : element_size (sizeof (t8_standalone_element<TEclass>)), scheme_context (sc_mempool_new (element_size)) {};

 protected:
  size_t element_size;  /**< The size in bytes of an element of class \a eclass */
  void *scheme_context; /**< Anonymous implementation context. */

 public:
  /** Destructor for all default schemes */
  ~t8_standalone_scheme ()
  {
    T8_ASSERT (scheme_context != NULL);
    SC_ASSERT (((sc_mempool_t *) scheme_context)->elem_count == 0);
    sc_mempool_destroy ((sc_mempool_t *) scheme_context);
  }

  /** Move constructor */
  t8_standalone_scheme (t8_standalone_scheme &&other) noexcept
    : element_size (other.element_size), scheme_context (other.scheme_context)
  {
    other.scheme_context = nullptr;
  }

  /** Move assignment operator */
  t8_standalone_scheme &
  operator= (t8_standalone_scheme &&other) noexcept
  {
    if (this != &other) {
      // Free existing resources of moved-to object
      if (scheme_context) {
        sc_mempool_destroy ((sc_mempool_t *) scheme_context);
      }

      // Transfer ownership of resources
      element_size = other.element_size;
      scheme_context = other.scheme_context;

      // Leave the source object in a valid state
      other.scheme_context = nullptr;
    }
    return *this;
  }

  /** Copy constructor */
  t8_standalone_scheme (const t8_standalone_scheme &other)
    : element_size (other.element_size), scheme_context (sc_mempool_new (other.element_size)) {};

  /** Copy assignment operator */
  t8_standalone_scheme &
  operator= (const t8_standalone_scheme &other)
  {
    if (this != &other) {
      // Free existing resources of assigned-to object
      if (scheme_context) {
        sc_mempool_destroy ((sc_mempool_t *) scheme_context);
      }

      // Copy the values from the source object
      element_size = other.element_size;
      scheme_context = sc_mempool_new (other.element_size);
    }
    return *this;
  }

  // ################################################____GENERAL INFO____################################################

  /** Return the tree class of this scheme.
   * \return The tree class of this scheme.
   */
  constexpr t8_eclass_t
  get_eclass (void) const
  {
    return TEclass;
  }

  constexpr size_t
  get_element_size (void) const
  {
    return element_size;
  }

  /** Returns true, if there is one element in the tree, that does not refine into 2^dim children.
   * Returns false otherwise.
   * \return                    non-zero if there is one element in the tree that does not refine into 2^dim children.
   */
  static constexpr int
  refines_irregular (void)
  {
    if constexpr (TEclass == T8_ECLASS_PYRAMID) {
      return 1;
    }
    return 0;
  }

  /** Return the maximum allowed level for any element of a given class.
   * \return                      The maximum allowed level for elements of class \b ts.
   */
  static constexpr int
  get_maxlevel (void)
  {
    return T8_ELEMENT_MAXLEVEL[TEclass];
  }

  // ################################################____SHAPE INFORMATION____################################################

  /** Compute the number of corners of a given element.
   * \param [in] elem The element.
   * \return          The number of corners of \a elem.
   */
  static constexpr int
  element_get_num_corners (const t8_element_t *elem)
  {
    T8_ASSERT (element_is_valid (elem));

    return T8_ELEMENT_NUM_CORNERS[TEclass];
  }

  /** Compute the number of faces of a given element.
   * \param [in] elem The element.
   * \return          The number of faces of \a elem.
   */
  static constexpr int
  element_get_num_faces (const t8_element_t *elem)
  {
    SC_ABORT ("This function is not implemented in this scheme yet.\n");
    return 0;
  }

  /** Compute the maximum number of faces of a given element and all of its
   *  descendants.
   * \param [in] elem The element.
   * \return          The maximum number of faces of \a elem and its descendants.
   */
  static constexpr int
  element_get_max_num_faces (const t8_element_t *elem)
  {
    SC_ABORT ("This function is not implemented in this scheme yet.\n");
    return 0;
  }

  static constexpr t8_element_shape_t
  element_get_shape (const t8_element_t *elem)
  {
    T8_ASSERT (element_is_valid (elem));
    return TEclass;
  }

  /** Return the corner number of an element's face corner.
   * Example quad: 2 x --- x 3
   *                 |     |
   *                 |     |   face 1
   *               0 x --- x 1
   *      Thus for face = 1 the output is: corner=0 : 1, corner=1: 3
   *
   * \param [in] element  The element.
   * \param [in] face     A face index for \a element.
   * \param [in] corner   A corner index for the face 0 <= \a corner < num_face_corners.
   * \return              The corner number of the \a corner-th vertex of \a face.
   *
   * The order in which the corners must be given is determined by the eclass of \a element:
   * LINE/QUAD/TRIANGLE:  No specific order.
   * HEX               :  In Z-order of the face starting with the lowest corner number.
   * TET               :  Starting with the lowest corner number counterclockwise as seen from
   *                      'outside' of the element.
   */
  static constexpr int
  element_get_face_corner (const t8_element_t *element, const int face, const int corner)
  {
    SC_ABORT ("This function is not implemented yet.\n");
    return 0;
  }

  /** Return the face numbers of the faces sharing an element's corner.
   * Example quad: 2 x --- x 3
   *                 |     |
   *                 |     |   face 1
   *               0 x --- x 1
   *                  face 2
   *      Thus for corner = 1 the output is: face=0 : 2, face=1: 1
   * \param [in] element  The element.
   * \param [in] corner   A corner index for the face.
   * \param [in] face     A face index for \a corner.
   * \return              The face number of the \a face-th face at \a corner.
   */
  static constexpr int
  element_get_corner_face (const t8_element_t *element, const int corner, const int face)
  {
    SC_ABORT ("This function is not implemented yet.\n");
    return 0;
  }

  static constexpr t8_element_shape_t
  element_get_face_shape (const t8_element_t *elem, const int face)
  {
    SC_ABORT ("This function is not implemented yet.\n");
    return T8_ECLASS_ZERO;
  }

  // ################################################____GENERAL HELPER____################################################

  /** Copy all entries of \b source to \b dest. \b dest must be an existing
   *  element. No memory is allocated by this function.
   * \param [in] source The element whose entries will be copied to \b dest.
   * \param [in,out] dest This element's entries will be overwrite with the
   *                    entries of \b source.
   * \note \a source and \a dest may point to the same element.
   */
  static constexpr void
  element_copy (const t8_element_t *source, t8_element_t *dest)
  {
    T8_ASSERT (element_is_valid (source));
    if (source == dest)
      return;
    memcpy ((t8_standalone_element<TEclass> *) dest, (const t8_standalone_element<TEclass> *) source,
            sizeof (t8_standalone_element<TEclass>));
    T8_ASSERT (element_is_valid (dest));
  }

  /** Check if two elements are equal.
  * \param [in] elem1  The first element.
  * \param [in] elem2  The second element.
  * \return            true if the elements are equal, false if they are not equal
  */
  static constexpr int
  element_is_equal (const t8_element_t *elem1, const t8_element_t *elem2) noexcept
  {
    T8_ASSERT (element_is_valid (elem1));
    T8_ASSERT (element_is_valid (elem2));

    const t8_standalone_element<TEclass> *el1 = (const t8_standalone_element<TEclass> *) elem1;
    const t8_standalone_element<TEclass> *el2 = (const t8_standalone_element<TEclass> *) elem2;
    if (el1->level != el2->level)
      return 0;
    for (int idim = 0; idim < T8_ELEMENT_DIM[TEclass]; idim++) {
      if (el1->coords[idim] != el2->coords[idim])
        return 0;
    }
    /* return el1->type == el2->type;
    ToDo-Type */
    return 1;
  }

  // ################################################____ACCESSOR____################################################

  /** Return the level of a particular element.
   * \param [in] elem    The element whose level should be returned.
   * \return             The level of \b elem.
   */
  static constexpr int
  element_get_level (const t8_element_t *elem)
  {
    T8_ASSERT (element_is_valid (elem));
    return ((const t8_standalone_element<TEclass> *) elem)->level;
  }

  // ################################################____REFINEMENT____################################################

  /** create the root element
   * \param [in,out] elem The element that is filled with the root
   */
  static constexpr void
  get_root (t8_element_t *elem)
  {
    t8_standalone_element<TEclass> *el = (t8_standalone_element<TEclass> *) elem;
    el->level = 0;
    for (int idim = 0; idim < T8_ELEMENT_DIM[TEclass]; idim++) {
      el->coords[idim] = 0;
    }
    /* el->type = 0;
    ToDo-Type */
    return;
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
  static constexpr void
  element_get_parent (const t8_element_t *elem, t8_element_t *parent)
  {
    T8_ASSERT (element_is_valid (elem));

    const t8_standalone_element<TEclass> *el = (const t8_standalone_element<TEclass> *) elem;
    t8_standalone_element<TEclass> *parent_elem = (t8_standalone_element<TEclass> *) parent;

    T8_ASSERT (el->level > 0);

    if constexpr (T8_ELEMENT_NUM_EQUATIONS[TEclass]) {
      SC_ABORT ("Only implemented for hypercubes.\n");
    }

    const t8_element_coord length = element_get_len ((el->level));
    set_coords_at_level_to_zero (el, parent_elem, length);

    parent_elem->level = el->level - 1;
    T8_ASSERT (parent_elem->level >= 0);

    T8_ASSERT (element_is_valid (parent));
  }

  /** Compute the number of siblings of an element. That is the number of 
   * elements with the same parent (if available).
   * \param [in] elem The element.
   * \return          The number of siblings of \a element.
   * Note that this number is >= 1, since we count the element itself as a sibling.
   * Note that the number of siblings is 1 for the root element.
   */
  static constexpr int
  element_get_num_siblings (const t8_element_t *elem)
  {
    T8_ASSERT (element_is_valid (elem));

    const t8_standalone_element<TEclass> *el = (const t8_standalone_element<TEclass> *) elem;
    if (el->level == 0)
      return 1;
    T8_ASSERT (0 < el->level && el->level <= T8_ELEMENT_MAXLEVEL[TEclass]);
    /* To get the number siblings, we first get the parent and then get the number of children of that parent*/
    if constexpr (refines_irregular ()) {
      SC_ABORT ("This function is not implemented yet.\n");
    }
    else {
      return T8_ELEMENT_NUM_CHILDREN[TEclass];
    }
  }

  /** Compute a specific sibling of a given element \b elem and store it in \b sibling.
   *  \b sibling needs to be an existing element. No memory is allocated by this function.
   *  \b elem and \b sibling can point to the same element, then the entries of
   *  \b elem are overwritten by the ones of its sibid-th sibling.
   * \param [in] elem   The element whose sibling will be computed.
   * \param [in] sibid  The id of the sibling computed.
   * \param [in,out] sibling This element's entries will be overwritten by those
   *                    of \b elem's sibid-th sibling.
   *                    The storage for this element must exist
   *                    and match the element class of the sibling.
   */
  static constexpr void
  element_get_sibling (const t8_element_t *elem, const int sibid, t8_element_t *sibling)
  {
    SC_ABORT ("This function is not implemented yet.\n");
  }

  /** Construct the child element of a given number.
   * \param [in] elem     This must be a valid element, bigger than maxlevel.
   * \param [in] childid  The number of the child to construct.
   * \param [in,out] child        The storage for this element must exist
   *                              and match the element class of the child.
   *                              For a pyramid, for example, it may be either a
   *                              tetrahedron or a pyramid depending on \a childid.
   *                              This can be checked by \a t8_element_child_eclass.
   *                              On output, a valid element.
   * It is valid to call this function with elem = child.
   * \see t8_element_child_eclass
   */
  static constexpr void
  element_get_child (const t8_element_t *elem, const int childid, t8_element_t *child)
  {
    T8_ASSERT (element_is_valid (elem));
    T8_ASSERT (0 <= childid);
    T8_ASSERT (childid < element_get_num_children (elem));

    const t8_standalone_element<TEclass> *el = (const t8_standalone_element<TEclass> *) elem;
    t8_standalone_element<TEclass> *c = (t8_standalone_element<TEclass> *) child;

    T8_ASSERT (0 <= childid && childid < T8_ELEMENT_NUM_CHILDREN[TEclass]);
    T8_ASSERT (0 <= el->level && el->level <= T8_ELEMENT_MAXLEVEL[TEclass]);

    /* Compute the cube id and shift the coordinates accordingly */
    t8_cube_id cube_id;
    if constexpr (T8_ELEMENT_NUM_EQUATIONS[TEclass]) {
      SC_ABORT ("Only implemented for hypercubes.\n");
    }
    else {
      cube_id = childid;
    }

    const t8_element_coord length = element_get_len (el->level + 1);

    put_cube_id_at_level (el, c, length, cube_id);

    c->level = el->level + 1;

    T8_ASSERT (element_is_valid (child));
  }

  /** Return the number of children of an element when it is refined.
   * \param [in] elem   The element whose number of children is returned.
   * \return            The number of children of \a elem if it is to be refined.
   */
  static constexpr int
  element_get_num_children (const t8_element_t *elem)
  {
    T8_ASSERT (element_is_valid (elem));

    return T8_ELEMENT_NUM_CHILDREN[TEclass];
  }

  /** Construct all children of a given element.
   * \param [in] elem     This must be a valid element, bigger than maxlevel.
   * \param [in] length   The length of the output array \a c must match
   *                      the number of children.
   * \param [in,out] c    The storage for these \a length elements must exist
   *                      and match the element class in the children's ordering.
   *                      On output, all children are valid.
   * It is valid to call this function with elem = c[0].
   * \see t8_element_num_children
   * \see t8_element_child_eclass
   */
  static constexpr void
  element_get_children (const t8_element_t *elem, const int length, t8_element_t *c[])
  {
    T8_ASSERT (element_is_valid (elem));

    const t8_standalone_element<TEclass> *el = (const t8_standalone_element<TEclass> *) elem;
    T8_ASSERT (0 <= el->level && el->level <= T8_ELEMENT_MAXLEVEL[TEclass]);

    const int num_children = length;
    T8_ASSERT (length == element_get_num_children ((const t8_element_t *) el));

    for (int ichild = num_children - 1; ichild >= 0; ichild--) {
      element_get_child ((const t8_element_t *) el, ichild, c[ichild]);
      T8_ASSERT (element_is_valid (c[ichild]));
    }
  }

  /** Compute the child id of an element.
   * \param [in] elem     This must be a valid element.
   * \return              The child id of elem.
   */
  static constexpr int
  element_get_child_id (const t8_element_t *elem)
  {
    T8_ASSERT (element_is_valid (elem));

    const t8_standalone_element<TEclass> *el = (const t8_standalone_element<TEclass> *) elem;
    T8_ASSERT (el->level >= 0);
    if (el->level == 0) {
      return -1;
    }
    const t8_cube_id cube_id = compute_cubeid (el, el->level);
    t8_child_id child_id;
    if constexpr (T8_ELEMENT_NUM_EQUATIONS[TEclass]) {
      SC_ABORT ("Only implemented for hypercubes.\n");
    }
    else {
      child_id = cube_id;
    }
    return child_id;
  }

  /** Compute the ancestor id of an element, that is the child id
   * at a given level.
   * \param [in] elem     This must be a valid element.
   * \param [in] level    A refinement level. Must satisfy \a level < elem.level
   * \return              The child_id of \a elem in regard to its \a level ancestor.
   */
  static constexpr int
  element_get_ancestor_id (const t8_element_t *elem, const int level)
  {
    T8_ASSERT (element_is_valid (elem));
    T8_ASSERT (0 <= level && level <= T8_ELEMENT_MAXLEVEL[TEclass]);

    const t8_standalone_element<TEclass> *el = (const t8_standalone_element<TEclass> *) elem;
    t8_standalone_element<TEclass> ancestor;
    T8_ASSERT (0 <= el->level && el->level <= T8_ELEMENT_MAXLEVEL[TEclass]);

    element_get_ancestor (el, level, &ancestor);
    return element_get_child_id ((const t8_element_t *) &ancestor);
  }

  /** Query whether a given set of elements is a family or not.
   * \param [in] fam      An array of as many elements as an element of class
   *                      \b ts has siblings.
   * \return              Zero if \b fam is not a family, nonzero if it is.
   * \note level 0 elements do not form a family.
   */
  static constexpr int
  elements_are_family (t8_element_t *const *fam)
  {
#if T8_ENABLE_DEBUG
    const int num_siblings = element_get_num_siblings (fam[0]);
    for (int isib = 0; isib < num_siblings; isib++) {
      T8_ASSERT (element_is_valid (fam[isib]));
    }
#endif

    t8_standalone_element<TEclass> parent, compare;
    /* Take the parent of the first element as baseline to compare against */
    element_get_parent ((const t8_element_t *) fam[0], (t8_element_t *) &parent);
    const int num_children = element_get_num_children ((const t8_element_t *) &parent);
    for (int childid = 0; childid < num_children; childid++) {
      /* check whether each element has the same parent */
      element_get_parent ((const t8_element_t *) fam[childid], (t8_element_t *) &compare);
      if (element_compare ((const t8_element_t *) &parent, (const t8_element_t *) &compare)) {
        return 0;
      }

      /* check whether each element is the correct child of the collective parent */
      /* Could be replaced by type comparison as level is already checked in parent comparison */
      element_get_child ((const t8_element_t *) &parent, childid, (t8_element_t *) &compare);

      if (element_compare ((const t8_element_t *) fam[childid], (const t8_element_t *) &compare)) {
        return 0;
      }
    }
    return 1;
  }

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
  static constexpr void
  element_get_nca (const t8_element_t *elem1, const t8_element_t *elem2, t8_element_t *nca)
  {
    T8_ASSERT (element_is_valid (elem1));
    T8_ASSERT (element_is_valid (elem2));

    const t8_standalone_element<TEclass> *el1 = (const t8_standalone_element<TEclass> *) elem1;
    const t8_standalone_element<TEclass> *el2 = (const t8_standalone_element<TEclass> *) elem2;
    /* get the first possible level of the nca*/
    int cube_ancestor_level = element_get_cube_nca_level (el1, el2);
    int real_level;
    if constexpr (T8_ELEMENT_NUM_EQUATIONS[TEclass]) {
      SC_ABORT ("Only implemented for hypercubes.\n");
    }
    else {
      real_level = cube_ancestor_level;
    }
    /* get the ancestor at the calculated level*/
    element_get_ancestor (el1, real_level, (t8_standalone_element<TEclass> *) nca);
    T8_ASSERT (element_is_valid (nca));
  }

  /** Compute the first descendant of a given element.
   * \param [in] elem     The element whose descendant is computed.
   * \param [out] desc    The first element in a uniform refinement of \a elem
   *                      of the given level.
   * \param [in] level    The level, at which the descendant is computed.
   */
  static constexpr void
  element_get_first_descendant (const t8_element_t *elem, t8_element_t *desc, const int level)
  {
    T8_ASSERT (element_is_valid (elem));

    const t8_standalone_element<TEclass> *el = (const t8_standalone_element<TEclass> *) elem;
    t8_standalone_element<TEclass> *d = (t8_standalone_element<TEclass> *) desc;

    T8_ASSERT (level >= el->level);
    T8_ASSERT (0 <= level && level <= T8_ELEMENT_MAXLEVEL[TEclass]);

    /* The first descendant of an element has the same anchor coords and type, but another level */
    element_copy ((const t8_element_t *) el, (t8_element_t *) d);
    d->level = level;

    T8_ASSERT (element_is_valid ((t8_element_t *) d));
  }

  /** Compute the last descendant of a given element.
   * \param [in] elem     The element whose descendant is computed.
   * \param [out] desc    The last element in a uniform refinement of \a elem
   *                      of the given level.
   * \param [in] level    The level, at which the descendant is computed.
   */
  static constexpr void
  element_get_last_descendant (const t8_element_t *elem, t8_element_t *desc, const int level)
  {
    T8_ASSERT (element_is_valid (elem));

    const t8_standalone_element<TEclass> *el = (const t8_standalone_element<TEclass> *) elem;
    t8_standalone_element<TEclass> *d = (t8_standalone_element<TEclass> *) desc;

    T8_ASSERT (level >= el->level);
    T8_ASSERT (0 <= level && level <= T8_ELEMENT_MAXLEVEL[TEclass]);

    element_copy ((const t8_element_t *) el, (t8_element_t *) d);
    d->level = level;

    /* Shift the coords to the eighth cube. The type of the last descendant
    * is the type of the input element */
    t8_element_coord coord_offset = element_get_len (el->level) - element_get_len (level);
    for (int idim = 0; idim < T8_ELEMENT_DIM[TEclass]; idim++) {
      d->coords[idim] |= coord_offset;
    }

    T8_ASSERT (element_is_valid (desc));
  }

  // ################################################____FACE REFINEMENT____################################################

  /** Return the number of children of an element's face when the element is refined.
   * \param [in] elem   The element whose face is considered.
   * \param [in] face   A face of \a elem.
   * \return            The number of children of \a face if \a elem is to be refined.
   */
  static constexpr int
  element_get_num_face_children (const t8_element_t *elem, const int face)
  {
    SC_ABORT ("This function is not implemented in this scheme yet.\n");
    return 0;
  }

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
  static constexpr void
  element_get_children_at_face (const t8_element_t *elem, int face, t8_element_t *children[], const int num_children,
                                int *child_indices)
  {
    SC_ABORT ("This function is not implemented in this scheme yet.\n");
  }

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
  static constexpr int
  element_face_get_child_face (const t8_element_t *elem, const int face, const int face_child)
  {
    SC_ABORT ("This function is not implemented in this scheme yet.\n");
    return 0;
  }

  /** Given a face of an element return the face number
   * of the parent of the element that matches the element's face. Or return -1 if
   * no face of the parent matches the face.

    * \param [in]  elem    The element.
    * \param [in]  face    Then number of the face.
    * \return              If \a face of \a elem is also a face of \a elem's parent,
    *                      the face number of this face. Otherwise -1.
    * \note For the root element this function always returns \a face.
    */
  static constexpr int
  element_face_get_parent_face (const t8_element_t *elem, const int face)
  {
    SC_ABORT ("This function is not implemented in this scheme yet.\n");
    return 0;
  }

  /** Construct the first descendant of an element at a given level that touches a given face.
   * \param [in] elem      The input element.
   * \param [in] face      A face of \a elem.
   * \param [in, out] first_desc An allocated element. This element's data will be
   *                       filled with the data of the first descendant of \a elem
   *                       that shares a face with \a face.
   * \param [in] level     The level, at which the first descendant is constructed
   */
  static constexpr void
  element_get_first_descendant_face (const t8_element_t *elem, const int face, t8_element_t *first_desc,
                                     const int level)
  {
    SC_ABORT ("This function is not implemented in this scheme yet.\n");
  }

  /** Construct the last descendant of an element at a given level that touches a given face.
   * \param [in] elem      The input element.
   * \param [in] face      A face of \a elem.
   * \param [in, out] last_desc An allocated element. This element's data will be
   *                       filled with the data of the last descendant of \a elem
   *                       that shares a face with \a face.
   * \param [in] level     The level, at which the last descendant is constructed
   */
  static constexpr void
  element_get_last_descendant_face (const t8_element_t *elem, const int face, t8_element_t *last_desc, const int level)
  {
    SC_ABORT ("This function is not implemented in this scheme yet.\n");
  }

  // ################################################____FACE NEIGHBOR____################################################

  /** Compute whether a given element shares a given face with its root tree.
   * \param [in] elem     The input element.
   * \param [in] face     A face of \a elem.
   * \return              True if \a face is a subface of the element's root element.
   * \note You can compute the corresponding face number of the tree via \ref t8_element_tree_face.
   */
  static constexpr int
  element_is_root_boundary (const t8_element_t *elem, const int face)
  {
    SC_ABORT ("This function is not implemented in this scheme yet.\n");
    return 0;
  }

  /** Given an element and a face of this element. If the face lies on the
   *  tree boundary, return the face number of the tree face.
   *  If not the return value is arbitrary.
   *  You can call \ref t8_element_is_root_boundary to query whether the face is
   *  at the tree boundary.
   * \param [in] elem     The element.
   * \param [in] face     The index of a face of \a elem.
   * \return The index of the tree face that \a face is a subface of, if
   *         \a face is on a tree boundary.
   *         Any arbitrary integer if \a is not at a tree boundary.
   * \warning The return value may look like a valid face of the tree even if 
   *   the element does not lie on the root boundary.
   */
  static constexpr int
  element_get_tree_face (const t8_element_t *elem, const int face)
  {
    SC_ABORT ("This function is not implemented in this scheme yet.\n");
    return 0;
  }

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
  static constexpr int
  element_get_face_neighbor_inside (const t8_element_t *elem, t8_element_t *neigh, const int face, int *neigh_face)
  {
    SC_ABORT ("This function is not implemented in this scheme yet.\n");
    return 0;
  }

  // ################################################____TREE FACE TRANSFORMATION____################################################  */

  /** Suppose we have two trees that share a common face f.
   *  Given an element e that is a subface of f in one of the trees
   *  and given the orientation of the tree connection, construct the face
   *  element of the respective tree neighbor that logically coincides with e
   *  but lies in the coordinate system of the neighbor tree.
   *  \param [in] elem1     The face element.
   *  \param [in,out] elem2 On return the face element \a elem1 with respective
   *                        to the coordinate system of the other tree.
   *  \param [in] orientation The orientation of the tree-tree connection.
   *                        \see t8_cmesh_set_join
   *  \param [in] sign      Depending on the topological orientation of the two tree faces,
   *                        either 0 (both faces have opposite orientation)
   *                        or 1 (both faces have the same top. orientation).
   *                        \ref t8_eclass_face_orientation
   *  \param [in] is_smaller_face Flag to declare whether \a elem1 belongs to
   *                        the smaller face. A face f of tree T is smaller than
   *                        f' of T' if either the eclass of T is smaller or if
   *                        the classes are equal and f<f'. The orientation is
   *                        defined in relation to the smaller face.
   * \note \a elem1 and \a elem2 may point to the same element.
   */
  static constexpr void
  element_transform_face (const t8_element_t *elem1, t8_element_t *elem2, const int orientation, const int sign,
                          const int is_smaller_face)
  {
    SC_ABORT ("This function is not implemented in this scheme yet.\n");
  }

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
  static constexpr int
  element_extrude_face (const t8_element_t *face, t8_element_t *elem, const int root_face, const t8_scheme *face_scheme)
  {
    SC_ABORT ("This function is not implemented in this scheme yet.\n");
    return 0;
  }

  /** Construct the boundary element at a specific face.
   * \param [in] elem     The input element.
   * \param [in] face     The index of the face of which to construct the
   *                      boundary element.
   * \param [in,out] boundary An allocated element of dimension of \a element
   *                      minus 1. The entries will be filled with the entries
   *                      of the face of \a element.
   * \param [in] boundary_scheme The scheme for the eclass of the boundary face.
   * If \a elem is of class T8_ECLASS_VERTEX, then \a boundary must be NULL
   * and will not be modified.
   */
  static constexpr void
  element_get_boundary_face (const t8_element_t *elem, const int face, t8_element_t *boundary,
                             const t8_scheme *boundary_scheme)
  {
    SC_ABORT ("This function is not implemented in this scheme yet.\n");
  }

  // ################################################____LINEAR ID____################################################

  /** Initialize the entries of an allocated element according to a
   *  given linear id in a uniform refinement.
   * \param [in,out] elem The element whose entries will be set.
   * \param [in] level    The level of the uniform refinement to consider.
   * \param [in] id       The linear id.
   *                      id must fulfil 0 <= id < 'number of leaves in the uniform refinement'
   */
  static constexpr void
  element_set_linear_id (t8_element_t *elem, const int level, t8_linearidx_t id)
  {

    t8_standalone_element<TEclass> *el = (t8_standalone_element<TEclass> *) elem;

    get_root ((t8_element_t *) el);

    /* There is only one element at level 0, so it must be root */
    if (level == 0) {
      T8_ASSERT (id == 0);
      return;
    }

    T8_ASSERT (0 <= id);
    T8_ASSERT (1 <= level && level <= T8_ELEMENT_MAXLEVEL[TEclass]);
    t8_standalone_element<TEclass> child;

    while (el->level < level) {
      /* Shortcut if we need the first descendant of the subtree*/
      if (id == 0) {
        element_get_first_descendant ((const t8_element_t *) el, (t8_element_t *) el, level);
        return;
      }

      t8_linearidx_t sum_descendants_of_children_before;
      t8_linearidx_t sum_descendants_of_children_until_current = 0;
      int childindex = -1;

      /* Find the first child id so that the sum of descendants of previous child and the own number of descendants is greater than id */
      do {
        /* Go to the next child */
        sum_descendants_of_children_before = sum_descendants_of_children_until_current;
        childindex++;
        T8_ASSERT (childindex < element_get_num_children ((const t8_element_t *) el));

        element_get_num_children ((const t8_element_t *) el);

        element_get_child ((const t8_element_t *) el, childindex, (t8_element_t *) &child);
        const t8_linearidx_t num_descendants_of_child = element_count_leaves ((t8_element_t *) &child, level);

        /* Add number of descendant of current child to cumulative sum */
        sum_descendants_of_children_until_current = sum_descendants_of_children_before + num_descendants_of_child;

      } while (sum_descendants_of_children_until_current <= id);

      /* Replace el by child to go into next iteration at finer level*/
      element_get_child ((const t8_element_t *) el, childindex, (t8_element_t *) el);
      /* get id in subtree of child */
      id -= sum_descendants_of_children_before;
    }
    T8_ASSERT (id < T8_ELEMENT_NUM_CHILDREN[TEclass]);
    element_get_child ((const t8_element_t *) el, id, (t8_element_t *) el);
    return;
  }

  /** Compute the linear id of a given element in a hypothetical uniform
   * refinement of a given level.
   * \param [in] elem     The element whose id we compute.
   * \param [in] level    The level of the uniform refinement to consider.
   * \return              The linear id of the element.
   */
  static constexpr t8_linearidx_t
  element_get_linear_id (const t8_element_t *elem, const int level)
  {
    T8_ASSERT (element_is_valid (elem));
    const t8_standalone_element<TEclass> *el = (const t8_standalone_element<TEclass> *) elem;
    t8_standalone_element<TEclass> ancestor;

    /* Determine the starting element for the iterative linear id computation. */
    if (level < el->level) {
      /* Throw away child ids up to the coarser level */
      element_get_ancestor (el, level, &ancestor);
    }
    else {
      /* Start with the input element. 
       Copy to have a mutable element. */
      element_copy ((const t8_element_t *) el, (t8_element_t *) &ancestor);
    }

    t8_linearidx_t id = 0;
    t8_standalone_element<TEclass> child;

    while (ancestor.level != 0) {
      const t8_child_id childid = element_get_child_id ((t8_element_t *) &ancestor);
      element_get_parent ((t8_element_t *) &ancestor, (t8_element_t *) &ancestor);
      t8_linearidx_t parent_id = 0;

      for (int ichild = 0; ichild < childid; ichild++) {
        /* el is now parent, so compute child to get sibling of previous el */

        element_get_child ((const t8_element_t *) &ancestor, ichild, (t8_element_t *) &child);
        const t8_linearidx_t num_child_descendants = element_count_leaves ((t8_element_t *) &child, level);
        parent_id += num_child_descendants;
      }
      id += parent_id;
    }
    T8_ASSERT (id >= 0);
    return id;
  }

  /** Construct the successor in a uniform refinement of a given element.
   * \param [in] elem1    The element whose successor should be constructed.
   * \param [in,out] elem2  The element whose entries will be set.
   */
  static constexpr void
  element_construct_successor (const t8_element_t *elem1, t8_element_t *elem2)
  {
    T8_ASSERT (element_is_valid (elem1));

    const t8_standalone_element<TEclass> *elem = (const t8_standalone_element<TEclass> *) elem1;
    t8_standalone_element<TEclass> *succ = (t8_standalone_element<TEclass> *) elem2;

    element_copy ((const t8_element_t *) elem, (t8_element_t *) succ);

    const t8_child_id child_id = element_get_child_id ((const t8_element_t *) elem);
    const int num_siblings = element_get_num_siblings ((const t8_element_t *) elem);
    T8_ASSERT (0 <= child_id && child_id < num_siblings);
    /* If the element is the last child of the parent, we need to go to the parent's successor (go to a coarser level)*/
    if (child_id == num_siblings - 1) {
      element_get_parent ((const t8_element_t *) succ, (t8_element_t *) succ);
      element_construct_successor ((const t8_element_t *) succ, (t8_element_t *) succ);
      element_get_child ((const t8_element_t *) succ, 0, (t8_element_t *) succ);
    }
    else {
      element_get_parent ((const t8_element_t *) succ, (t8_element_t *) succ);
      element_get_child ((const t8_element_t *) succ, child_id + 1, (t8_element_t *) succ);
    }

    T8_ASSERT (element_is_valid (elem2));
  }

  /** Count how many leaf descendants of a given uniform level an element would produce.
   * \param [in] t     The element to be checked.
   * \param [in] level A refinement level. 
   * \return Suppose \a t is uniformly refined up to level \a level. The return value
   * is the resulting number of elements (of the given level).
   * If \a level < t8_element_level(t), the return value should be 0.
   *
   * Example: If \a t is a line element that refines into 2 line elements on each level,
   *  then the return value is max(0, 2^{\a level - level(\a t)}).
   *  Thus, if \a t's level is 0, and \a level = 3, the return value is 2^3 = 8.
   */
  static constexpr t8_gloidx_t
  element_count_leaves (const t8_element_t *elem, const t8_element_level level)
  {
    T8_ASSERT (element_is_valid (elem));
    T8_ASSERT (0 <= level && level <= T8_ELEMENT_MAXLEVEL[TEclass]);
    if (level < element_get_level (elem)) {
      return 0;
    }
    else {
      return num_descendants_at_leveldiff (elem, level - element_get_level (elem));
    }
  }

  /** Count how many leaf descendants of a given uniform level the root element will produce.
   * \param [in] level A refinement level.
   * \return The value of \ref t8_element_count_leaves if the input element
   *      is the root (level 0) element.
   *
   * This is a convenience function, and can be implemented via
   * \ref t8_element_count_leaves.
   */
  static constexpr t8_gloidx_t
  count_leaves_from_root (const int level)
  {
    T8_ASSERT (level <= T8_ELEMENT_MAXLEVEL[TEclass]);
    T8_ASSERT (level >= 0);
    if constexpr (TEclass == T8_ECLASS_PYRAMID) {
      SC_ABORT ("Not implemented yet.\n");
    }
    return 1LL << (level * T8_ELEMENT_DIM[TEclass]);
  }

  /** Compare two elements.
   * \param [in] elem1  The first element.
   * \param [in] elem2  The second element.
   * \return       negative if elem1 < elem2, zero if elem1 equals elem2
   *               and positive if elem1 > elem2.
   *  If elem2 is a copy of elem1 then the elements are equal.
   */
  static constexpr int
  element_compare (const t8_element_t *elem1, const t8_element_t *elem2)
  {
    T8_ASSERT (element_is_valid (elem1));
    T8_ASSERT (element_is_valid (elem2));

    const t8_standalone_element<TEclass> *e1 = (const t8_standalone_element<TEclass> *) elem1;
    const t8_standalone_element<TEclass> *e2 = (const t8_standalone_element<TEclass> *) elem2;

    const int maxlvl = SC_MAX (e1->level, e2->level);

    const t8_linearidx_t id1 = element_get_linear_id ((const t8_element_t *) e1, maxlvl);
    const t8_linearidx_t id2 = element_get_linear_id ((const t8_element_t *) e2, maxlvl);
    if (id1 == id2) {
      if (e1->level == e2->level) {
        return 0;
      }
      else {
        return e1->level - e2->level;
      }
    }
    return id1 < id2 ? -1 : 1;
  }

  // ################################################____VISUALIZATION____################################################

  /** Compute the coordinates of a given element vertex inside a reference tree
   *  that is embedded into [0,1]^d (d = dimension).
   *   \param [in] elem      The element to be considered.
   *   \param [in] vertex The id of the vertex whose coordinates shall be computed.
   *   \param [out] coords An array of at least as many doubles as the element's dimension
   *                      whose entries will be filled with the coordinates of \a vertex.
   */
  static constexpr void
  element_get_vertex_reference_coords (const t8_element_t *elem, const int vertex, double coords[])
  {
    T8_ASSERT (element_is_valid (elem));

    const t8_standalone_element<TEclass> *el = (const t8_standalone_element<TEclass> *) elem;
    if constexpr (TEclass == T8_ECLASS_VERTEX) {
      return;
    }
    else {
      int coords_int[T8_ELEMENT_DIM[TEclass]];
      T8_ASSERT (0 <= vertex && vertex < T8_ELEMENT_NUM_CORNERS[TEclass]);
      element_compute_coords (el, vertex, coords_int);
      for (int idim = 0; idim < T8_ELEMENT_DIM[TEclass]; idim++) {
        coords[idim] = coords_int[idim] / (double) get_root_len ();
      }
    }
  }

  /** Convert a point in the reference space of an element to a point in the
   *  reference space of the tree.
   * 
   * \param [in] elem         The element.
   * \param [in] coords_input The coordinates of the point in the reference space of the element.
   * \param [in] user_data    User data.
   * \param [out] out_coords  The coordinates of the point in the reference space of the tree.
   */
  static constexpr void
  element_get_reference_coords (const t8_element_t *elem, const double *ref_coords, const size_t num_coords,
                                double *out_coords)
  {
    double *current_ref_coords = (double *) ref_coords;
    double *current_out_coords = out_coords;
    t8_element_coord length = element_get_len (element_get_level (elem));

    for (size_t coord = 0; coord < num_coords; ++coord) {
      for (int dim = 0; dim < T8_ELEMENT_DIM[TEclass]; ++dim) {
        current_out_coords[dim]
          = ((t8_standalone_element<TEclass> *) elem)->coords[dim] + current_ref_coords[dim] * length;

        current_out_coords[dim] /= (double) get_root_len ();
      }

      current_ref_coords += T8_ECLASS_MAX_DIM;
      current_out_coords += T8_ELEMENT_DIM[TEclass];
    }
  }

  // ################################################____MEMORY____################################################

  /** Allocate memory for an array of elements of a given class and initialize them.
   * \param [in] length   The number of elements to be allocated.
   * \param [in,out] elems On input an array of \b length many unallocated
   *                      element pointers.
   *                      On output all these pointers will point to an allocated
   *                      and initialized element.
   * \note Not every element that is created in t8code will be created by a call
   * to this function. However, if an element is not created using \ref t8_element_new,
   * then it is guaranteed that \ref t8_element_init is called on it.
   * \note In debugging mode, an element that was created with \ref t8_element_new
   * must pass \ref t8_element_is_valid.
   * \note If an element was created by \ref t8_element_new then \ref t8_element_init
   * may not be called for it. Thus, \ref t8_element_new should initialize an element
   * in the same way as a call to \ref t8_element_init would.
   * \see t8_element_init
   * \see t8_element_is_valid
   */
  /* TODO: would it be better to directly allocate an array of elements,
   *       not element pointers? */
  void
  element_new (const int length, t8_element_t **elem) const
  {
    /* allocate memory */
    T8_ASSERT (this->scheme_context != NULL);
    T8_ASSERT (0 <= length);
    T8_ASSERT (elem != NULL);

    for (int i = 0; i < length; ++i) {
      elem[i] = (t8_element_t *) sc_mempool_alloc ((sc_mempool_t *) this->scheme_context);
    }

/* in debug mode, set sensible default values. */
#ifdef T8_ENABLE_DEBUG
    {
      int i;
      for (i = 0; i < length; i++) {
        element_init (1, elem[i]);
      }
    }
#endif
  }

  /** Initialize an array of allocated elements.
   * \param [in] length   The number of elements to be initialized.
   * \param [in,out] elems On input an array of \b length many allocated
   *                       elements.
   * \note In debugging mode, an element that was passed to \ref t8_element_init
   * must pass \ref t8_element_is_valid.
   * \note If an element was created by \ref t8_element_new then \ref t8_element_init
   * may not be called for it. Thus, \ref t8_element_new should initialize an element
   * in the same way as a call to \ref t8_element_init would.
   * Thus, if \a called_new is true this function should usually do nothing.
   * \see t8_element_new
   * \see t8_element_is_valid
   */
  static inline void
  element_init (const int length, t8_element_t *elem)
  {
#ifdef T8_ENABLE_DEBUG
    int ielem;
    t8_standalone_element<TEclass> *el = (t8_standalone_element<TEclass> *) elem;
    /* Set all values to 0 */
    for (ielem = 0; ielem < length; ielem++) {
      element_set_linear_id ((t8_element_t *) (el + ielem), 0, 0);
      T8_ASSERT (element_is_valid ((t8_element_t *) (el + ielem)));
    }
#endif
  }

  static constexpr void
  element_deinit (const int length, t8_element_t *elem)
  {
  }

  /** Deallocate an array of elements.
   * \param [in] length   The number of elements in the array.
   * \param [in,out] elems On input an array of \b length many allocated
   *                      element pointers.
   *                      On output all these pointers will be freed.
   *                      \b elem itself will not be freed by this function.
   */
  void
  element_destroy (const int length, t8_element_t **elem) const
  {
    T8_ASSERT (this->scheme_context != NULL);
    T8_ASSERT (0 <= length);
    T8_ASSERT (elem != NULL);
    for (int i = 0; i < length; ++i) {
      sc_mempool_free ((sc_mempool_t *) scheme_context, elem[i]);
    }
  }

  // ################################################____DEBUG____################################################

#ifdef T8_ENABLE_DEBUG
  /** Query whether a given element can be considered as 'valid' and it is
   *  safe to perform any of the above algorithms on it.
   *  For example this could mean that all coordinates are in valid ranges
   *  and other membervariables do have meaningful values.
   * \param [in]      elem  The element to be checked.
   * \return          True if \a elem is safe to use. False otherwise.
   * \note            An element that is constructed with \ref t8_element_new
   *                  must pass this test.
   * \note            An element for which \ref t8_element_init was called must pass
   *                  this test.
   * \note            This function is used for debugging to catch certain errors.
   *                  These can for example occur when an element points to a region
   *                  of memory which should not be interpreted as an element.
   * \note            We recommend to use the assertion T8_ASSERT (element_is_valid (elem))
   *                  in the implementation of each of the functions in this file.
   */
  static constexpr int
  element_is_valid (const t8_element_t *elem)
  {
    T8_ASSERT (elem != NULL);

    const t8_standalone_element<TEclass> *el = (const t8_standalone_element<TEclass> *) elem;
    const t8_element_coord max_coord = get_root_len () - 1;

    /* Check the level */
    int is_valid = 0 <= el->level && el->level <= T8_ELEMENT_MAXLEVEL[TEclass];
    /* Check coordinates, we allow a boundary layer around the root-element */
    for (int i = 0; i < T8_ELEMENT_DIM[TEclass]; i++) {
      is_valid = is_valid && -(int64_t) get_root_len () <= el->coords[i] && el->coords[i] <= max_coord;
    }

    return is_valid;
  }

  /**
 * Print a given element. For a example for a triangle print the coordinates
 * and the level of the triangle. This function is only available in the
 * debugging configuration. 
 * 
 * \param [in]        elem  The element to print
 */
  static constexpr void
  element_debug_print (const t8_element_t *elem)
  {

    const t8_standalone_element<TEclass> *el = (const t8_standalone_element<TEclass> *) elem;

    t8_debugf ("level: %i\n", el->level);
    for (int i = 0; i < T8_ELEMENT_DIM[TEclass]; i++) {
      t8_debugf ("x_%i: %i \n", i, el->coords[i]);
    }
    /**  for (int e = 0; e < T8_ELEMENT_NUM_EQUATIONS[TEclass]; e++) {
    *  t8_debugf ("t_%i: %i \n", e, el->type[e]);
    *}
    * ToDo-Type */
  }

  static constexpr void
  element_to_string (const t8_element_t *elem, char *debug_string, const int string_size)
  {
    const t8_standalone_element<TEclass> *el = (const t8_standalone_element<TEclass> *) elem;
    int offset = 0;
    offset += snprintf (debug_string + offset, string_size - offset, "level: %i\n", el->level);
    for (int idim = 0; idim < T8_ELEMENT_DIM[TEclass]; idim++) {
      offset += snprintf (debug_string + offset, string_size - offset, "x_%i: %i \n", idim, el->coords[idim]);
    }
  }

#endif

  // ################################################____MPI____################################################

  /** Pack multiple elements into contiguous memory, so they can be sent via MPI.
     * \param [in] elements Array of elements that are to be packed
     * \param [in] count Number of elements to pack
     * \param [in,out] send_buffer Buffer in which to pack the elements
     * \param [in] buffer_size size of the buffer (in order to check that we don't access out of range)
     * \param [in, out] position the position of the first byte that is not already packed
     * \param [in] comm MPI Communicator
    */
  constexpr void
  element_MPI_Pack (t8_element_t **const elements, const unsigned int count, void *send_buffer, const int buffer_size,
                    int *position, sc_MPI_Comm comm) const

  {
    int mpiret;
    t8_standalone_element<TEclass> **els = (t8_standalone_element<TEclass> **) elements;

    for (unsigned int ielem = 0; ielem < count; ielem++) {
      for (int idim = 0; idim < T8_ELEMENT_DIM[TEclass]; idim++) {
        mpiret = sc_MPI_Pack (&(els[ielem]->coords[idim]), 1, sc_MPI_INT, send_buffer, buffer_size, position, comm);
        SC_CHECK_MPI (mpiret);
      }
      mpiret = sc_MPI_Pack (&els[ielem]->level, 1, sc_MPI_INT8_T, send_buffer, buffer_size, position, comm);
      SC_CHECK_MPI (mpiret);
    }
  }

  /** Determine an upper bound for the size of the packed message of \a count elements
     * \param [in] count Number of elements to pack
     * \param [in] comm MPI Communicator
     * \param [out] pack_size upper bound on the message size
    */
  constexpr void
  element_MPI_Pack_size (const unsigned int count, sc_MPI_Comm comm, int *pack_size) const
  {
    int singlesize = 0;
    int datasize = 0;
    int mpiret;

    /* x,y,z */
    mpiret = sc_MPI_Pack_size (1, sc_MPI_INT, comm, &datasize);
    SC_CHECK_MPI (mpiret);
    singlesize += T8_ELEMENT_DIM[TEclass] * datasize;

    /* level */
    mpiret = sc_MPI_Pack_size (1, sc_MPI_INT8_T, comm, &datasize);
    SC_CHECK_MPI (mpiret);
    singlesize += datasize;

    *pack_size = count * singlesize;
  }

  /** Unpack multiple elements from contiguous memory that was received via MPI.
     * \param [in] recvbuf Buffer from which to unpack the elements
     * \param [in] buffer_size size of the buffer (in order to check that we don't access out of range)
     * \param [in, out] position the position of the first byte that is not already packed
     * \param [in] elements Array of initialised elements that is to be filled from the message
     * \param [in] count Number of elements to unpack
     * \param [in] comm MPI Communicator
    */
  constexpr void
  element_MPI_Unpack (void *recvbuf, const int buffer_size, int *position, t8_element_t **elements,
                      const unsigned int count, sc_MPI_Comm comm) const
  {
    int mpiret;
    t8_standalone_element<TEclass> **els = (t8_standalone_element<TEclass> **) elements;

    for (unsigned int ielem = 0; ielem < count; ielem++) {
      for (int idim = 0; idim < T8_ELEMENT_DIM[TEclass]; idim++) {
        mpiret = sc_MPI_Unpack (recvbuf, buffer_size, position, &(els[ielem]->coords[idim]), 1, sc_MPI_INT, comm);
        SC_CHECK_MPI (mpiret);
      }
      mpiret = sc_MPI_Unpack (recvbuf, buffer_size, position, &(els[ielem]->level), 1, sc_MPI_INT8_T, comm);
      SC_CHECK_MPI (mpiret);
    }
  }

 private:
  // ################################################____HELPER____################################################
  static constexpr t8_element_coord
  element_get_len (const t8_element_level level)
  {
    return 1 << (T8_ELEMENT_MAXLEVEL[TEclass] - (level));
  }

  static t8_cube_id
  compute_cubeid (const t8_standalone_element<TEclass> *elem, const int level)
  {
    t8_cube_id cube_id = 0;

    T8_ASSERT (0 <= elem->level && elem->level <= T8_ELEMENT_MAXLEVEL[TEclass]);
    const t8_element_coord h = element_get_len (level);

    /* The cube id of the root element is 0.*/
    if (level == 0) {
      return 0;
    }
    for (int i = 0; i < T8_ELEMENT_DIM[TEclass]; i++) {
      cube_id |= ((elem->coords[i] & h) ? 1 << i : 0);
    }
    return cube_id;
  }

  /**
 * Compute the ancestor of \a el at a given level via the equation properties
 * 
 * \param[in] elem      Input element
 * \param[in] level     Level of the ancestor to compute
 * \param[in, out] and  Allocated element that will be filled with the data of the ancestor.
 */
  static constexpr void
  element_get_ancestor (const t8_standalone_element<TEclass> *elem, const int level,
                        t8_standalone_element<TEclass> *ancestor)
  {
    T8_ASSERT (0 <= level && level <= elem->level);
    if (elem != ancestor) {
      element_copy ((const t8_element_t *) elem, (t8_element_t *) ancestor);
    }
    if (elem->level == level) {
      return;
    }

    /* Set type */
    if constexpr (T8_ELEMENT_NUM_EQUATIONS[TEclass]) {
      SC_ABORT ("Only implemented for hypercubes.\n");
    }

    /* The coordinates and the type of the ancestor are defined by the level. */
    element_cut_coordinates (ancestor, T8_ELEMENT_MAXLEVEL[TEclass] - level);

    ancestor->level = level;
  }

  static constexpr int
  element_get_cube_nca_level (const t8_standalone_element<TEclass> *elem1, const t8_standalone_element<TEclass> *elem2)
  {
    /* XOR all coordinates. The number of zeros on the left determines the level needed, so that the coordinates equal. 
    OR over all these bit representations. The number of zeros on the left in this new number equals the coarses of all of these levels. 
    Therefore this is the level needed so that all coordinates equal.*/
    t8_element_coord maxexclor = 0;

    for (int idim = 0; idim < T8_ELEMENT_DIM[TEclass]; idim++) {
      maxexclor |= (elem1->coords[idim] ^ elem2->coords[idim]);
    }

    const int num_zeros = number_of_leading_zeros (maxexclor);
    /* If one element already is the ancestor of the other element num_zeros evaluates to maxlevel, in that case return the coarser of both levels*/
    return SC_MIN (num_zeros, (int) SC_MIN (elem1->level, elem2->level));
  }

  static constexpr int
  number_of_leading_zeros (const t8_element_coord maxexclor)
  {
    const int num_of_active_bits_used = SC_LOG2_32 (maxexclor) + 1;
    T8_ASSERT (num_of_active_bits_used <= T8_ELEMENT_MAXLEVEL[TEclass]);

    return T8_ELEMENT_MAXLEVEL[TEclass] - num_of_active_bits_used;
  }

  /**
 * Set the \a shift last bits of every coordinate to zero. 
 * 
 * \param[in, out]  elem     Input element
 * \param[in]       shift Number of bits to set to zero
 */
  static constexpr void
  element_cut_coordinates (t8_standalone_element<TEclass> *elem, const int shift)
  {
    T8_ASSERT (0 <= shift && shift <= T8_ELEMENT_MAXLEVEL[TEclass]);
    for (int idim = 0; idim < T8_ELEMENT_DIM[TEclass]; idim++) {
      elem->coords[idim] = (elem->coords[idim] >> shift) << shift;
    }
  }

  /**
 * Set the least significant coordinates bits to zero. 
 * 
 * \param[in]  elem        Input element
 * \param[in, out]       parent_elem Parent element
 * \param[in]       length      int that is 1 at the level of the input element
 * Note length is used as additional input to avoid recomputation. 
 */
  static constexpr void
  set_coords_at_level_to_zero (const t8_standalone_element<TEclass> *elem, t8_standalone_element<TEclass> *parent_elem,
                               const t8_element_coord length)
  {
    for (int idim = 0; idim < T8_ELEMENT_DIM[TEclass]; idim++) {
      parent_elem->coords[idim] = elem->coords[idim] & ~length;
    }
  }

  /**
 * Adjust the coordinates based on the cube ID.
 * 
 * \param[in]           parent       Input element
 * \param[in, out]      child     Output element
 * \param[in]           length   int that is 1 at the level of the child element 
 * \param[in]           cube_id  Cube ID for bitwise operation
 * Note length is used as additional input to avoid recomputation. 
 */
  static constexpr void
  put_cube_id_at_level (const t8_standalone_element<TEclass> *parent, t8_standalone_element<TEclass> *child,
                        const t8_element_coord length, const t8_cube_id cube_id)
  {
    for (int idim = 0; idim < T8_ELEMENT_DIM[TEclass]; idim++) {
      child->coords[idim] = parent->coords[idim] + ((cube_id & (1 << idim)) ? length : 0);
    }
  }

  static constexpr t8_element_coord
  get_root_len ()
  {
    if constexpr (TEclass == T8_ECLASS_VERTEX) {
      return 0;
    }
    else {
      return 1 << T8_ELEMENT_MAXLEVEL[TEclass];
    }
  }

  /**Caller is responsible for taking the absolute value of leveldiff */
  static constexpr t8_linearidx_t
  num_descendants_at_leveldiff (const t8_element_t *elem, const t8_element_level leveldiff)
  {
    if (leveldiff < 0)
      return 0;
    if constexpr (TEclass == T8_ECLASS_PYRAMID) {
      SC_ABORT ("Not implemented yet.\n");
    }
    return 1LL << (T8_ELEMENT_DIM[TEclass] * leveldiff);
  }

  /** Compute the coordinates of a vertex of an element.
 * \param [in] elem    Input element.
 * \param [in] vertex The number of the vertex.
 * \param [out] coords An array of 3 t8_element_coord that
 * 		     will be filled with the coordinates of the vertex.
 */
  static constexpr void
  element_compute_coords (const t8_standalone_element<TEclass> *elem, const int vertex, int coords[])
  {
    T8_ASSERT (0 <= vertex && vertex < element_get_num_corners ((const t8_element_t *) elem));

    if constexpr (T8_ELEMENT_NUM_EQUATIONS[TEclass]) {
      SC_ABORT ("Only implemented for hypercubes.\n");
    }
    else {
      //Hypercubes
      for (int idim = 0; idim < T8_ELEMENT_DIM[TEclass]; idim++) {
        coords[idim] = elem->coords[idim] + ((vertex & (1 << idim)) >> idim) * element_get_len (elem->level);
      }
    }
  }
};

#endif /* T8_STANDALONE_IMPLEMENTATION_HXX */
