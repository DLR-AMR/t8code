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

/** \file t8_scheme_implementation.hxx
 *  An implementation for the class \ref t8_scheme in \ref t8_scheme.hxx.
 */
#pragma once

#include <t8.h>
#include <t8_schemes/t8_scheme.hxx>
#include <t8_eclass/t8_eclass.h>
#include <sc_functions.h>
#include <t8_schemes/t8_standalone/t8_standalone_implementation.hxx>
#include <t8_schemes/t8_standalone/t8_standalone_elements.hxx>
#include <t8_schemes/t8_subelement/t8_subelement_type.hxx>
#include <t8_schemes/t8_scheme_helpers.hxx>
#include <utility>
#include <algorithm>

/** TODO. */
struct t8_subelementquad_scheme: public t8_scheme_helpers<T8_ECLASS_QUAD, t8_subelementquad_scheme>
{
 public:
  using standalone_scheme = t8_standalone_scheme<T8_ECLASS_QUAD>;
  /** Constructor
  */
  t8_subelementquad_scheme () noexcept
    : element_size (sizeof (t8_subelement_element)), scheme_context (sc_mempool_new (element_size)) {};

 protected:
  // What do i need this for?
  size_t element_size;  /**< The size in bytes of an element of class \a eclass */
  void *scheme_context; /**< Anonymous implementation context. */

 public:
  /** Destructor for all default schemes */
  ~t8_subelementquad_scheme ()
  {
    T8_ASSERT (scheme_context != NULL);
    SC_ASSERT (((sc_mempool_t *) scheme_context)->elem_count == 0);
    sc_mempool_destroy ((sc_mempool_t *) scheme_context);
  }

  /** Move constructor */
  t8_subelementquad_scheme (t8_subelementquad_scheme &&other) noexcept
    : element_size (other.element_size), scheme_context (std::exchange (other.scheme_context, nullptr))
  {
  }

  /** Move assignment operator */
  t8_subelementquad_scheme &
  operator= (t8_subelementquad_scheme &&other) noexcept
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
  t8_subelementquad_scheme (const t8_subelementquad_scheme &other)
    : element_size (other.element_size), scheme_context (sc_mempool_new (other.element_size)) {};

  /** Copy assignment operator */
  t8_subelementquad_scheme &
  operator= (const t8_subelementquad_scheme &other)
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

  /** Return the size of any element of a given class.
   * \return                      The size of an element.
   */
  static constexpr size_t
  get_element_size (void) noexcept
  {
    return sizeof (t8_subelement_element);
  }

  /** Returns true, if there is one element in the tree, that does not refine into 2^dim children.
   * Returns false otherwise.
   * \return                    non-zero if there is one element in the tree that does not refine into 2^dim children.
   */
  static constexpr int
  refines_irregular (void) noexcept
  {
    return true;  // Potentially there are subelements.
  }

  /** Return the maximum allowed level for any element of a given class.
   * \return                      The maximum allowed level for elements of class \b ts.
   */
  static constexpr int
  get_maxlevel (void) noexcept
  {
    return standalone_scheme::get_maxlevel () - 1;  // We need to reserve one level for the subelements.
  }

  // ################################################____SHAPE INFORMATION____################################################

  /** Compute the number of corners of a given element.
   * \param [in] elem The element.
   * \return          The number of corners of \a elem.
   */
  static int
  element_get_num_corners (const t8_element_t *elem) noexcept
  {
    T8_ASSERT (element_is_valid (elem));
    const t8_subelement_element *subelement = (const t8_subelement_element *) elem;
    if (subelement->subelement_type == 0) {
      return standalone_scheme::element_get_num_corners (subelement_to_element (subelement));
    }

    return T8_ELEMENT_NUM_CORNERS[T8_ECLASS_TRIANGLE];
  }

  /** Compute the number of faces of a given element.
   * \param [in] elem The element.
   * \return          The number of faces of \a elem.
   */
  static int
  element_get_num_faces (const t8_element_t *elem) noexcept
  {
    T8_ASSERT (element_is_valid (elem));
    const t8_subelement_element *subelement = (const t8_subelement_element *) elem;
    if (subelement->subelement_type == 0) {
      return standalone_scheme::element_get_num_faces (subelement_to_element (subelement));
    }

    return T8_ELEMENT_NUM_FACES[T8_ECLASS_TRIANGLE];
  }

  /** Compute the maximum number of faces of a given element and all of its
   *  descendants.
   * \param [in] elem The element.
   * \return          The maximum number of faces of \a elem and its descendants.
   */
  static int
  element_get_max_num_faces (const t8_element_t *elem) noexcept
  {
    T8_ASSERT (element_is_valid (elem));
    const t8_subelement_element *subelement = (const t8_subelement_element *) elem;
    // As subelements are discarded for the next adaptation cycle, descendants are also quads.
    return standalone_scheme::element_get_max_num_faces (subelement_to_element (subelement));
  }

  /** Return the shape of an allocated element according its type.
   * For example, a child of an element can be an element of a different shape
   * and has to be handled differently - according to its shape.
   * \param [in] elem     The element to be considered
   * \return              The shape of the element as an eclass
   */
  static t8_element_shape_t
  element_get_shape (const t8_element_t *elem) noexcept
  {
    T8_ASSERT (element_is_valid (elem));
    if (!element_is_subelement (elem)) {
      return T8_ECLASS_QUAD;
    }

    return T8_ECLASS_TRIANGLE;
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
  static int
  element_get_face_corner ([[maybe_unused]] const t8_element_t *element, [[maybe_unused]] const int face,
                           [[maybe_unused]] const int corner) noexcept
  {
    SC_ABORT ("element_get_face_corner is not implemented for subelements yet.\n");
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
  static int
  element_get_corner_face ([[maybe_unused]] const t8_element_t *element, [[maybe_unused]] const int corner,
                           [[maybe_unused]] const int face) noexcept
  {
    SC_ABORT ("element_get_face_corner is not implemented for subelements yet.\n");
  }

  /** Compute the shape of the face of an element.
   * \param [in] elem     The element.
   * \param [in] face     A face of \a elem.
   * \return              The element shape of the face.
   * I.e. T8_ECLASS_LINE for quads, T8_ECLASS_TRIANGLE for tets
   *      and depending on the face number either T8_ECLASS_QUAD or
   *      T8_ECLASS_TRIANGLE for prisms.
   */
  static t8_element_shape_t
  element_get_face_shape ([[maybe_unused]] const t8_element_t *elem, [[maybe_unused]] const int face) noexcept
  {
    T8_ASSERT (element_is_valid (elem));
    if (element_is_subelement (elem)) {
      T8_ASSERT (0 <= face && face < T8_ELEMENT_NUM_FACES[T8_ECLASS_TRIANGLE]);
    }
    else {
      T8_ASSERT (0 <= face && face < T8_ELEMENT_NUM_FACES[T8_ECLASS_QUAD]);
    }
    return T8_ECLASS_LINE;
  }

  // ################################################____GENERAL HELPER____################################################

  /** Copy all entries of \b source to \b dest. \b dest must be an existing
   *  element. No memory is allocated by this function.
   * \param [in] source The element whose entries will be copied to \b dest.
   * \param [in,out] dest This element's entries will be overwrite with the
   *                    entries of \b source.
   * \note \a source and \a dest may point to the same element.
   */
  static void
  element_copy (const t8_element_t *source, t8_element_t *dest) noexcept
  {
    T8_ASSERT (element_is_valid (source));
    if (source == dest)
      return;
    memcpy ((t8_subelement_element *) dest, (const t8_subelement_element *) source, sizeof (t8_subelement_element));
    T8_ASSERT (element_is_valid (dest));
  }

  /** Check if two elements are equal.
  * \param [in] elem1  The first element.
  * \param [in] elem2  The second element.
  * \return            true if the elements are equal, false if they are not equal
  */
  static int
  element_is_equal (const t8_element_t *elem1, const t8_element_t *elem2) noexcept
  {
    T8_ASSERT (element_is_valid (elem1));
    T8_ASSERT (element_is_valid (elem2));

    const t8_subelement_element *el1 = (const t8_subelement_element *) elem1;
    const t8_subelement_element *el2 = (const t8_subelement_element *) elem2;
    if (el1->subelement_type != el2->subelement_type) {
      return 0;
    }
    if (el1->subelement_id != el2->subelement_id) {
      return 0;
    }
    return standalone_scheme::element_is_equal (subelement_to_element (el1), subelement_to_element (el2));
  }

  // ################################################____ACCESSOR____################################################

  /** Return the level of a particular element.
   * \param [in] elem    The element whose level should be returned.
   * \return             The level of \b elem.
   */
  static int
  element_get_level (const t8_element_t *elem) noexcept
  {
    T8_ASSERT (element_is_valid (elem));
    return standalone_scheme::element_get_level (element_to_element (elem));
  }

  // ################################################____REFINEMENT____################################################

  /** create the root element
   * \param [in,out] elem The element that is filled with the root
   */
  static void
  set_to_root (t8_element_t *elem) noexcept
  {
    t8_subelement_element *subelement = (t8_subelement_element *) elem;
    reset_subelement_values (subelement);
    standalone_scheme::set_to_root (subelement_to_element (subelement));
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
  static void
  element_get_parent (const t8_element_t *elem, t8_element_t *parent) noexcept
  {
    T8_ASSERT (element_is_valid (elem));

    const t8_subelement_element *el = (const t8_subelement_element *) elem;
    t8_subelement_element *parent_elem = (t8_subelement_element *) parent;
    reset_subelement_values (parent_elem);
    if (element_is_subelement (elem)) {
      // For subelements, the parent is the element from which they are refined.
      standalone_scheme::element_copy (subelement_to_element (el), subelement_to_element (parent_elem));
      return;
    }
    standalone_scheme::element_get_parent (subelement_to_element (el), subelement_to_element (parent_elem));
  }

  /** Compute the number of siblings of an element. That is the number of
   * elements with the same parent (if available).
   * \param [in] elem The element.
   * \return          The number of siblings of \a element.
   * Note that this number is >= 1, since we count the element itself as a sibling.
   * Note that the number of siblings is 1 for the root element.
   */
  static int
  element_get_num_siblings (const t8_element_t *elem) noexcept
  {
    T8_ASSERT (element_is_valid (elem));
    const t8_subelement_element *subelement = (const t8_subelement_element *) elem;
    if (!element_is_subelement (elem)) {
      return standalone_scheme::element_get_num_siblings (subelement_to_element (subelement));
    }
    return element_get_number_of_subelements (subelement->subelement_type);
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
  static void
  element_get_sibling ([[maybe_unused]] const t8_element_t *elem, [[maybe_unused]] const int sibid,
                       [[maybe_unused]] t8_element_t *sibling) noexcept
  {
    SC_ABORT ("This function is not implemented yet.\n");
  }

  /** ONLY for non subelements!
   * Construct the child element of a given number.
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
  static void
  element_get_child (const t8_element_t *elem, const int childid, t8_element_t *child) noexcept
  {
    T8_ASSERT (!element_is_subelement (elem));
    T8_ASSERT (element_is_refinable (elem));
    standalone_scheme::element_get_child (element_to_element (elem), childid, element_to_element (child));
    T8_ASSERT (element_is_valid (child));
  }

  /** Return the number of children of an element when it is refined.
   * \param [in] elem   The element whose number of children is returned.
   * \return            The number of children of \a elem if it is to be refined.
   */
  static int
  element_get_num_children ([[maybe_unused]] const t8_element_t *elem) noexcept
  {
    /* Note that children of subelements equal the children of the parent quadrant. 
     * Therefore, the number of children of a subelement equals T8_ECLASS_QUAD. */
    T8_ASSERT (element_is_valid (elem));
    return T8_ELEMENT_NUM_CHILDREN[T8_ECLASS_QUAD];
  }

  /** Return the max number of children of an eclass.
   * \return            The max number of children of \a element.
   */
  static int
  get_max_num_children () noexcept
  {
    return T8_SUB_QUAD_MAX_SUBELEMENT_ID + 1;
  }

  /**
   * Indicates if an element is refinable. Possible reasons for being not refinable could be
   * that the element has reached its max level.
   * \param [in] elem   The element to check.
   * \return            True if the element is refinable.
   */
  static bool
  element_is_refinable (const t8_element_t *elem) noexcept
  {
    T8_ASSERT (element_is_valid (elem));
    if (element_is_subelement (elem)) {
      // Subelements are not refinable, as they are discarded for the next adaptation cycle.
      return false;
    }
    return standalone_scheme::element_get_level (element_to_element (elem)) < get_maxlevel ();
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
  static void
  element_get_children (const t8_element_t *elem, const int length, t8_element_t *c[]) noexcept
  {
    /* if elem is a subelement, then this function will construct the children of its parent. */
    T8_ASSERT (element_is_valid (elem));
    const t8_subelement_element *subelement = (const t8_subelement_element *) elem;

    t8_element_t *standalone_children_ptrs[T8_ELEMENT_NUM_CHILDREN[T8_ECLASS_QUAD]];
    for (int ichild = 0; ichild < length; ++ichild) {
      standalone_children_ptrs[ichild] = subelement_to_element ((t8_subelement_element *) c[ichild]);
    }

    if (element_is_subelement (elem)) {
      t8_subelement_element parent_storage;
      element_get_parent (elem, (t8_element_t *) &parent_storage);
      standalone_scheme::element_get_children (subelement_to_element (&parent_storage), length,
                                               standalone_children_ptrs);
    }
    else {
      standalone_scheme::element_get_children (subelement_to_element (subelement), length, standalone_children_ptrs);
    }
    for (int ichild = 0; ichild < length; ++ichild) {
      reset_subelement_values ((t8_subelement_element *) c[ichild]);
    }
  }

  /** Compute the child id of an element.
   * \param [in] elem     This must be a valid element.
   * \return              The child id of elem.
   */
  static int
  element_get_child_id (const t8_element_t *elem) noexcept
  {
    T8_ASSERT (element_is_valid (elem));
    const t8_subelement_element *subelement = (const t8_subelement_element *) elem;
    if (element_is_subelement (elem)) {
      // For subelements, the child id is the subelement id.
      return subelement->subelement_id;
    }
    return standalone_scheme::element_get_child_id (subelement_to_element (subelement));
  }

  /** Compute the ancestor id of an element, that is the child id
   * at a given level.
   * \param [in] elem     This must be a valid element.
   * \param [in] level    A refinement level. Must satisfy \a level < elem.level
   * \return              The child_id of \a elem in regard to its \a level ancestor.
   */
  static int
  element_get_ancestor_id (const t8_element_t *elem, const t8_element_level level) noexcept
  {
    SC_CHECK_ABORT (!element_is_subelement (elem), "element_get_ancestor_id is not implemented for subelements yet.\n");
    return standalone_scheme::element_get_ancestor_id (element_to_element (elem), level);
  }

  /** Query whether a given set of elements is a family or not.
   * \param [in] fam      An array of as many elements as an element of class
   *                      \b ts has siblings.
   * \return              Zero if \b fam is not a family, nonzero if it is.
   * \note level 0 elements do not form a family.
   */
  static int
  elements_are_family (t8_element_t *const *fam) noexcept
  {
#if T8_ENABLE_DEBUG
    const int num_siblings = element_get_num_siblings (fam[0]);
    for (int isib = 0; isib < num_siblings; isib++) {
      T8_ASSERT (element_is_valid (fam[isib]));
    }
#endif
    /* If the first element is a subelement, the remaining elements also have to be subelements and the elements must be equal. */
    if (element_is_subelement (fam[0])) {
      auto element_0 = element_to_element (fam[0]);
      for (int isib = 1; isib < element_get_num_siblings (fam[0]); ++isib) {
        if (!element_is_subelement (fam[isib])
            || !standalone_scheme::element_is_equal (element_0, element_to_element (fam[isib]))) {
          return 0;
        }
      }
      return 1;
    }
    /* If the first element is no subelement, the remaining elements also have to be no subelements and they must form a family. */
    t8_element_t *standalone_children[T8_ELEMENT_NUM_CHILDREN[T8_ECLASS_QUAD]];
    for (int isib = 0; isib < element_get_num_siblings (fam[0]); ++isib) {
      if (element_is_subelement (fam[isib])) {
        return 0;
      }
      standalone_children[isib] = element_to_element (fam[isib]);
    }
    return standalone_scheme::elements_are_family (standalone_children);
  }

  /** Query whether element A is an ancestor of the element B.
   * An element A is ancestor of an element B if A == B or if B can 
   * be obtained from A via successive refinement.
   * \param [in] element_A An element of class \a eclass in scheme \a scheme.
   * \param [in] element_B An element of class \a eclass in scheme \a scheme.
   * \return     True if and only if \a element_A is an ancestor of \a element_B.
  */
  static bool
  element_is_ancestor (const t8_element_t *element_A, const t8_element_t *element_B) noexcept
  {
    T8_ASSERT (element_is_valid (element_A));
    T8_ASSERT (element_is_valid (element_B));
    if (element_is_equal (element_A, element_B)) {
      return true;
    }
    if (element_is_subelement (element_A)) {
      // Subelements are not ancestors of any element, as they are discarded for the next adaptation cycle.
      // B could be a subelement if the underlying element is an ancestor of A.
      return false;
    }
    standalone_scheme tmp {};
    return tmp.element_is_ancestor (element_to_element (element_A), element_to_element (element_B));
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
  static void
  element_get_nca ([[maybe_unused]] const t8_element_t *elem1, [[maybe_unused]] const t8_element_t *elem2,
                   [[maybe_unused]] t8_element_t *nca) noexcept
  {
    SC_ABORT ("This function is not implemented yet.\n");
  }

  /** Compute the first descendant of a given element.
   * \param [in] elem     The element whose descendant is computed.
   * \param [out] desc    The first element in a uniform refinement of \a elem
   *                      of the given level.
   * \param [in] level    The level, at which the descendant is computed.
   */
  static void
  element_get_first_descendant (const t8_element_t *elem, t8_element_t *desc, const t8_element_level level) noexcept
  {
    SC_CHECK_ABORT (!element_is_subelement (elem),
                    "element_get_first_descendant is not implemented for subelements yet.\n");
    return standalone_scheme::element_get_first_descendant (element_to_element (elem), element_to_element (desc),
                                                            level);
  }

  /** Compute the last descendant of a given element.
   * \param [in] elem     The element whose descendant is computed.
   * \param [out] desc    The last element in a uniform refinement of \a elem
   *                      of the given level.
   * \param [in] level    The level, at which the descendant is computed.
   */
  static void
  element_get_last_descendant (const t8_element_t *elem, t8_element_t *desc, const t8_element_level level) noexcept
  {
    // Info: For subelements, the last descendant is the same as the first descendant, as they are discarded for the next adaptation cycle.
    standalone_scheme::element_get_last_descendant (element_to_element (elem), element_to_element (desc), level);
    reset_subelement_values ((t8_subelement_element *) desc);
  }

  // ################################################____FACE REFINEMENT____################################################

  /** Return the number of children of an element's face when the element is refined.
   * \param [in] elem   The element whose face is considered.
   * \param [in] face   A face of \a elem.
   * \return            The number of children of \a face if \a elem is to be refined.
   */
  static int
  element_get_num_face_children ([[maybe_unused]] const t8_element_t *elem, [[maybe_unused]] const int face) noexcept
  {
    SC_CHECK_ABORT (!element_is_subelement (elem),
                    "element_get_num_face_children is not implemented for subelements yet.\n");
    return standalone_scheme::element_get_num_face_children (element_to_element (elem), face);
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
   *                      \ref element_get_num_face_children
   * \param [in,out] child_indices If not NULL, an array of num_children integers must be given,
   *                      on output its i-th entry is the child_id of the i-th face_child.
   * It is valid to call this function with elem = children[0].
   */
  static void
  element_get_children_at_face ([[maybe_unused]] const t8_element_t *elem, const int face, t8_element_t *children[],
                                const int num_children, [[maybe_unused]] int *child_indices) noexcept
  {
    SC_CHECK_ABORT (!element_is_subelement (elem),
                    "element_get_children_at_face is not implemented for subelements yet.\n");
    t8_element_t *standalone_children_ptrs[T8_ELEMENT_NUM_CHILDREN[T8_ECLASS_QUAD]];
    for (int ichild = 0; ichild < num_children; ++ichild) {
      standalone_children_ptrs[ichild] = subelement_to_element ((t8_subelement_element *) children[ichild]);
    }
    standalone_scheme::element_get_children_at_face (element_to_element (elem), face, standalone_children_ptrs,
                                                     num_children, child_indices);
    for (int ichild = 0; ichild < num_children; ++ichild) {
      reset_subelement_values ((t8_subelement_element *) children[ichild]);
    }
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
    *                      the order of children from a call to \ref element_get_children_at_face.
    * \return              The face number of the face of a child of \a elem
    *                      that coincides with \a face_child.
    */
  static int
  element_face_get_child_face ([[maybe_unused]] const t8_element_t *elem, const int face,
                               [[maybe_unused]] const int face_child) noexcept
  {
    SC_CHECK_ABORT (!element_is_subelement (elem),
                    "element_face_get_child_face is not implemented for subelements yet.\n");
    return standalone_scheme::element_face_get_child_face (element_to_element (elem), face, face_child);
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
  static int
  element_face_get_parent_face ([[maybe_unused]] const t8_element_t *elem, const int face) noexcept
  {
    SC_CHECK_ABORT (!element_is_subelement (elem),
                    "element_face_get_parent_face is not implemented for subelements yet.\n");
    return standalone_scheme::element_face_get_parent_face (element_to_element (elem), face);
  }

  /** Construct the first descendant of an element at a given level that touches a given face.
   * \param [in] elem      The input element.
   * \param [in] face      A face of \a elem.
   * \param [in, out] first_desc An allocated element. This element's data will be
   *                       filled with the data of the first descendant of \a elem
   *                       that shares a face with \a face.
   * \param [in] level     The level, at which the first descendant is constructed
   */
  static void
  element_get_first_descendant_face (const t8_element_t *elem, const int face, t8_element_t *first_desc,
                                     const t8_element_level level) noexcept
  {
    SC_CHECK_ABORT (!element_is_subelement (elem),
                    "element_get_first_descendant_face is not implemented for subelements yet.\n");
    return standalone_scheme::element_get_first_descendant_face (element_to_element (elem), face,
                                                                 element_to_element (first_desc), level);
  }

  /** Construct the last descendant of an element at a given level that touches a given face.
   * \param [in] elem      The input element.
   * \param [in] face      A face of \a elem.
   * \param [in, out] last_desc An allocated element. This element's data will be
   *                       filled with the data of the last descendant of \a elem
   *                       that shares a face with \a face.
   * \param [in] level     The level, at which the last descendant is constructed
   */
  static void
  element_get_last_descendant_face ([[maybe_unused]] const t8_element_t *elem, const int face, t8_element_t *last_desc,
                                    const t8_element_level level) noexcept
  {
    SC_CHECK_ABORT (!element_is_subelement (elem),
                    "element_get_last_descendant_face is not implemented for subelements yet.\n");
    return standalone_scheme::element_get_last_descendant_face (element_to_element (elem), face,
                                                                element_to_element (last_desc), level);
  }

  // ################################################____FACE NEIGHBOR____################################################

  /** Compute whether a given element shares a given face with its root tree.
   * \param [in] elem     The input element.
   * \param [in] face     A face of \a elem.
   * \return              True if \a face is a subface of the element's root element.
   * \note You can compute the corresponding face number of the tree via \ref element_get_tree_face.
   */
  static int
  element_is_root_boundary (const t8_element_t *elem, const int face) noexcept
  {
    SC_CHECK_ABORT (!element_is_subelement (elem),
                    "element_is_root_boundary is not implemented for subelements yet.\n");
    return standalone_scheme::element_is_root_boundary (element_to_element (elem), face);
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
  static int
  element_get_tree_face (const t8_element_t *elem, const int face) noexcept
  {
    SC_CHECK_ABORT (!element_is_subelement (elem), "element_get_tree_face is not implemented for subelements yet.\n");
    return standalone_scheme::element_get_tree_face (element_to_element (elem), face);
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
  static int
  element_get_face_neighbor_inside (const t8_element_t *elem, t8_element_t *neigh, const int face,
                                    int *neigh_face) noexcept
  {
    SC_CHECK_ABORT (!element_is_subelement (elem),
                    "element_get_face_neighbor_inside is not implemented for subelements yet.\n");
    return standalone_scheme::element_get_face_neighbor_inside (element_to_element (elem), element_to_element (neigh),
                                                                face, neigh_face);
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
  static void
  element_transform_face ([[maybe_unused]] const t8_element_t *elem1, [[maybe_unused]] t8_element_t *elem2,
                          [[maybe_unused]] const int orientation, [[maybe_unused]] const int sign,
                          [[maybe_unused]] const int is_smaller_face) noexcept
  {
    /* This function has an explicit template specialization outside of t8_subelementquad_scheme*/
    SC_ABORT ("Not implemented for this eclass.\n");
  }

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
  static int
  element_extrude_face ([[maybe_unused]] const t8_element_t *face, [[maybe_unused]] t8_element_t *elem,
                        [[maybe_unused]] const int root_face, [[maybe_unused]] const t8_scheme *scheme) noexcept
  {
    SC_ABORT ("This function is nt implemented yet.");
  }

  /** Construct the boundary element at a specific face.
   * \param [in] elem     The input element.
   * \param [in] face     The index of the face of which to construct the
   *                      boundary element.
   * \param [in,out] boundary An allocated element of dimension of \a element
   *                      minus 1. The entries will be filled with the entries
   *                      of the face of \a element.
   * \param [in] scheme   The scheme containing an eclass scheme for the boundary face.
   * If \a elem is of class T8_ECLASS_VERTEX, then \a boundary must be NULL
   * and will not be modified.
   */
  static void
  element_get_boundary_face ([[maybe_unused]] const t8_element_t *elem, [[maybe_unused]] const int face,
                             [[maybe_unused]] t8_element_t *boundary, [[maybe_unused]] const t8_scheme *scheme) noexcept
  {
    SC_ABORT ("This function is not implemented yet.\n");
  }

  // ################################################____LINEAR ID____################################################

  /** Initialize the entries of an allocated element according to a
   *  given linear id in a uniform refinement.
   * \param [in,out] elem The element whose entries will be set.
   * \param [in] level    The level of the uniform refinement to consider.
   * \param [in] id       The linear id.
   *                      id must fulfil 0 <= id < 'number of leaves in the uniform refinement'
   */
  static void
  element_set_linear_id (t8_element_t *elem, const t8_element_level level, t8_linearidx_t id) noexcept
  {
    SC_CHECK_ABORT (!element_is_subelement (elem), "element_set_linear_id is not implemented for subelements yet.\n");
    standalone_scheme::element_set_linear_id (element_to_element (elem), level, id);
  }

  /** Compute the linear id of a given element in a hypothetical uniform
   * refinement of a given level.
   * \param [in] elem     The element whose id we compute.
   * \param [in] level    The level of the uniform refinement to consider.
   * \return              The linear id of the element.
   */
  static t8_linearidx_t
  element_get_linear_id (const t8_element_t *elem, const t8_element_level level) noexcept
  {
    SC_CHECK_ABORT (!element_is_subelement (elem), "element_get_linear_id is not implemented for subelements yet.\n");
    return standalone_scheme::element_get_linear_id (element_to_element (elem), level);
  }

  /** Construct the successor in a uniform refinement of a given element.
   * \param [in] elem1    The element whose successor should be constructed.
   * \param [in,out] elem2  The element whose entries will be set.
   */
  static void
  element_construct_successor (const t8_element_t *elem1, t8_element_t *elem2) noexcept
  {
    SC_CHECK_ABORT (!element_is_subelement (elem1),
                    "element_construct_successor is not implemented for subelements yet.\n");
    return standalone_scheme::element_construct_successor (element_to_element (elem1), element_to_element (elem2));
  }

  /** Count how many leaf descendants of a given uniform level an element would produce.
   * \param [in] elem  The element to be checked.
   * \param [in] level A refinement level.
   * \return Suppose \a elem is uniformly refined up to level \a level. The return value
   * is the resulting number of elements (of the given level).
   * If \a level < t8_element_level(t), the return value should be 0.
   *
   * Example: If \a elem is a line element that refines into 2 line elements on each level,
   *  then the return value is max(0, 2^{\a level - level(\a t)}).
   *  Thus, if \a elem's level is 0, and \a level = 3, the return value is 2^3 = 8.
   */
  static t8_gloidx_t
  element_count_leaves (const t8_element_t *elem, const t8_element_level level) noexcept
  {
    SC_CHECK_ABORT (!element_is_subelement (elem), "element_count_leaves is not implemented for subelements yet.\n");
    return standalone_scheme::element_count_leaves (element_to_element (elem), level);
  }

  /** Count how many leaf descendants of a given uniform level the root element will produce.
   * \param [in] level A refinement level.
   * \return The value of \ref t8_element_count_leaves if the input element
   *      is the root (level 0) element.
   *
   * This is a convenience function, and can be implemented via
   * \ref t8_element_count_leaves.
   */
  static t8_gloidx_t
  count_leaves_from_root (const t8_element_level level) noexcept
  {
    return standalone_scheme::count_leaves_from_root (level);
  }

  /** Compare two elements.
   * \param [in] elem1  The first element.
   * \param [in] elem2  The second element.
   * \return       negative if elem1 < elem2, zero if elem1 equals elem2
   *               and positive if elem1 > elem2.
   *  If elem2 is a copy of elem1 then the elements are equal.
   */
  static int
  element_compare (const t8_element_t *elem1, const t8_element_t *elem2) noexcept
  {
    SC_CHECK_ABORT (!element_is_subelement (elem1) && !element_is_subelement (elem2),
                    "element_compare is not implemented for subelements yet.\n");
    return standalone_scheme::element_compare (element_to_element (elem1), element_to_element (elem2));
  }

  // ################################################____VISUALIZATION____################################################

  /** Compute the coordinates of a given element vertex inside a reference tree
   *  that is embedded into [0,1]^d (d = dimension).
   *   \param [in] elem      The element to be considered.
   *   \param [in] vertex The id of the vertex whose coordinates shall be computed.
   *   \param [out] coords An array of at least as many doubles as the element's dimension
   *                      whose entries will be filled with the coordinates of \a vertex.
   */
  static void
  element_get_vertex_reference_coords (const t8_element_t *elem, const int vertex, double coords[]) noexcept
  {
    SC_CHECK_ABORT (!element_is_subelement (elem),
                    "element_get_vertex_reference_coords is not implemented for subelements yet.\n");
    return standalone_scheme::element_get_vertex_reference_coords (element_to_element (elem), vertex, coords);
  }

  /** Convert a point in the reference space of an element to a point in the
   *  reference space of the tree.
   *
   * \param [in] elem         The element.
   * \param [in] ref_coords   The coordinates of the point in the reference space of the element.
   * \param [in] num_coords   The number of coordinates to evaluate.
   * \param [out] out_coords  The coordinates of the point in the reference space of the tree.
   */
  static void
  element_get_reference_coords (const t8_element_t *elem, const double *ref_coords, const size_t num_coords,
                                double *out_coords) noexcept
  {
    SC_CHECK_ABORT (!element_is_subelement (elem),
                    "element_get_reference_coords is not implemented for subelements yet.\n");
    return standalone_scheme::element_get_reference_coords (element_to_element (elem), ref_coords, num_coords,
                                                            out_coords);
  }

  // ################################################____MEMORY____################################################

  /** Allocate memory for an array of elements of a given class and initialize them.
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
  /* TODO: would it be better to directly allocate an array of elements,
   *       not element pointers? */
  void
  element_new (const int length, t8_element_t **elems) const noexcept
  {
    /* allocate memory */
    T8_ASSERT (this->scheme_context != NULL);
    T8_ASSERT (0 <= length);
    T8_ASSERT (elems != NULL);

    for (int i = 0; i < length; ++i) {
      elems[i] = (t8_element_t *) sc_mempool_alloc ((sc_mempool_t *) this->scheme_context);
    }

/* in debug mode, set sensible default values. */
#if T8_ENABLE_DEBUG
    {
      for (int i = 0; i < length; i++) {
        element_init (1, elems[i]);
      }
    }
#endif
  }

  /** Initialize an array of allocated elements.
   * \param [in] length   The number of elements to be initialized.
   * \param [in,out] elems On input an array of \b length many allocated
   *                       elements.
   * \note In debugging mode, an element that was passed to \ref element_init
   * must pass \ref element_is_valid.
   * \note If an element was created by \ref element_new then \ref element_init
   * may not be called for it. Thus, \ref element_new should initialize an element
   * in the same way as a call to \ref element_init would.
   * \see element_new
   * \see element_is_valid
   */
  static void
  element_init ([[maybe_unused]] const int length, [[maybe_unused]] t8_element_t *elems) noexcept
  {
#if T8_ENABLE_DEBUG
    t8_subelement_element *subelement = (t8_subelement_element *) elems;
    for (int ielem = 0; ielem < length; ielem++) {
      reset_subelement_values (subelement + ielem);
      standalone_scheme::element_init (1, subelement_to_element (subelement + ielem));
      T8_ASSERT (element_is_valid ((t8_element_t *) (subelement + ielem)));
    }
#endif
  }

  /** Deinitialize an array of allocated elements.
   * \param [in] length   The number of elements to be deinitialized.
   * \param [in,out] elems On input an array of \a length many allocated
   *                       and initialized elements, on output an array of
   *                       \a length many allocated, but not initialized elements.
   * \note Call this function if you called element_init on the element pointers.
   * \see element_init
   */
  static constexpr void
  element_deinit ([[maybe_unused]] const int length, [[maybe_unused]] t8_element_t *elems) noexcept
  {
  }

  /** Deallocate an array of elements.
   * \param [in] length   The number of elements in the array.
   * \param [in,out] elems On input an array of \b length many allocated
   *                      element pointers.
   *                      On output all these pointers will be freed.
   *                      \b elems itself will not be freed by this function.
   */
  void
  element_destroy (const int length, t8_element_t **elems) const noexcept
  {
    T8_ASSERT (this->scheme_context != NULL);
    T8_ASSERT (0 <= length);
    T8_ASSERT (elems != NULL);
    for (int i = 0; i < length; ++i) {
      sc_mempool_free ((sc_mempool_t *) scheme_context, elems[i]);
    }
  }

  // ################################################____DEBUG____################################################

#if T8_ENABLE_DEBUG
  /** Query whether a given element can be considered as 'valid' and it is
   *  safe to perform any of the above algorithms on it.
   *  For example this could mean that all coordinates are in valid ranges
   *  and other membervariables do have meaningful values.
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
  static int
  element_is_valid (const t8_element_t *elem) noexcept
  {
    T8_ASSERT (elem != NULL);

    const t8_subelement_element *subelement = (const t8_subelement_element *) elem;
    int element_valid = standalone_scheme::element_is_valid (subelement_to_element (subelement));
    if (!element_is_subelement (elem)) {
      return element_valid;
    }
    bool subelement_valid = (subelement->subelement_type >= T8_SUB_QUAD_MIN_SUBELEMENT_TYPE
                             && subelement->subelement_type <= T8_SUB_QUAD_MAX_SUBELEMENT_TYPE)
                            && (subelement->subelement_id >= T8_SUB_QUAD_MIN_SUBELEMENT_ID
                                && subelement->subelement_id <= T8_SUB_QUAD_MAX_SUBELEMENT_ID);

    return subelement_valid && element_valid;
  }

  /**
  * Print a given element. For a example for a triangle print the coordinates
  * and the level of the triangle. This function is only available in the
  * debugging configuration.
  *
  * \param [in]        elem  The element to print
  */
  static void
  element_debug_print (const t8_element_t *elem) noexcept
  {
    const t8_subelement_element *subelement = (const t8_subelement_element *) elem;
    t8_debugf ("Subelement type: %i\n", subelement->subelement_type);
    t8_debugf ("Subelementid: %i\n", subelement->subelement_id);
    standalone_scheme::element_debug_print (subelement_to_element (subelement));
  }

#endif
  /**
 * Fill a string with readable information about the element
 * \param[in] elem The element to translate into human-readable information
 * \param[in, out] debug_string The string to fill.
 * \param[in] string_size Buffer size of c-string
 */
  static void
  element_to_string ([[maybe_unused]] const t8_element_t *elem, [[maybe_unused]] char *debug_string,
                     [[maybe_unused]] const int string_size) noexcept
  {
    SC_ABORT ("Not implemented.");
  }

  // ################################################____MPI____################################################

  /** Pack multiple elements into contiguous memory, so they can be sent via MPI.
     * \param [in] elements Array of elements that are to be packed
     * \param [in] count Number of elements to pack
     * \param [in,out] send_buffer Buffer in which to pack the elements
     * \param [in] buffer_size size of the buffer (in order to check that we don't access out of range)
     * \param [in, out] position the position of the first byte that is not already packed
     * \param [in] comm MPI Communicator
    */
  void
  element_MPI_Pack (t8_element_t **const elements, const unsigned int count, void *send_buffer, const int buffer_size,
                    int *position, sc_MPI_Comm comm) const noexcept

  {
    t8_subelement_element **els = (t8_subelement_element **) elements;
    standalone_scheme tmp {};
    for (unsigned int ielem = 0; ielem < count; ielem++) {
      t8_element_t *element = subelement_to_element (els[ielem]);
      tmp.element_MPI_Pack (&element, 1, send_buffer, buffer_size, position, comm);
      int mpiret = sc_MPI_Pack (&els[ielem]->subelement_type, 1, sc_MPI_INT, send_buffer, buffer_size, position, comm);
      SC_CHECK_MPI (mpiret);
      mpiret = sc_MPI_Pack (&els[ielem]->subelement_id, 1, sc_MPI_INT, send_buffer, buffer_size, position, comm);
      SC_CHECK_MPI (mpiret);
    }
  }

  /** Determine an upper bound for the size of the packed message of \a count elements
     * \param [in] count Number of elements to pack
     * \param [in] comm MPI Communicator
     * \param [out] pack_size upper bound on the message size
    */
  void
  element_MPI_Pack_size (const unsigned int count, sc_MPI_Comm comm, int *pack_size) const noexcept
  {
    // Get single size from standalone scheme.
    standalone_scheme tmp {};
    tmp.element_MPI_Pack_size (1, comm, pack_size);
    int singlesize = *pack_size;

    /* Type and id are both of type int. */
    int datasize = 0;
    int mpiret = sc_MPI_Pack_size (1, sc_MPI_INT, comm, &datasize);
    SC_CHECK_MPI (mpiret);
    singlesize += 2 * datasize;

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
  void
  element_MPI_Unpack (void *recvbuf, const int buffer_size, int *position, t8_element_t **elements,
                      const unsigned int count, sc_MPI_Comm comm) const noexcept
  {
    t8_subelement_element **els = (t8_subelement_element **) elements;
    standalone_scheme tmp {};
    for (unsigned int ielem = 0; ielem < count; ielem++) {
      t8_element_t *single = subelement_to_element (els[ielem]);
      tmp.element_MPI_Unpack (recvbuf, buffer_size, position, &single, 1, comm);
      int mpiret = sc_MPI_Unpack (recvbuf, buffer_size, position, &els[ielem]->subelement_type, 1, sc_MPI_INT, comm);
      SC_CHECK_MPI (mpiret);
      mpiret = sc_MPI_Unpack (recvbuf, buffer_size, position, &els[ielem]->subelement_id, 1, sc_MPI_INT, comm);
      SC_CHECK_MPI (mpiret);
    }
  }

  // --- Functions special for subelements ---

  static bool
  element_is_subelement (const t8_element_t *elem) noexcept
  {
    const t8_subelement_element *subelement = (const t8_subelement_element *) elem;
    return (subelement->subelement_type != 0);
  }
  static int
  element_get_number_of_subelements (int subelement_type)
  {

    int num_hanging_faces = 0;
    /* Count the number of ones of the binary subelement type. This number equals the number of hanging faces. */
    for (int i = 0; i < T8_ELEMENT_NUM_FACES[T8_ECLASS_QUAD]; ++i) {
      num_hanging_faces += (subelement_type & (1 << i)) >> i;
    }
    return T8_ELEMENT_NUM_FACES[T8_ECLASS_QUAD] + num_hanging_faces;
  }

  static void
  element_to_subelement (const t8_element_t *elem, int type, t8_element_t *c[])
  {
    const t8_subelement_element *element = (const t8_subelement_element *) elem;
    t8_subelement_element **subelements = (t8_subelement_element **) c;

    // const p4est_quadrant_t *q = &pquad_w_sub_elem->p4q;

    int num_subelements = element_get_number_of_subelements (type);

    T8_ASSERT (type >= T8_SUB_QUAD_MIN_SUBELEMENT_TYPE && type <= T8_SUB_QUAD_MAX_SUBELEMENT_TYPE);

    T8_ASSERT (!element_is_subelement (elem));
    T8_ASSERT (element_is_valid (elem));
#if T8_ENABLE_DEBUG
    {
      for (int j = 0; j < num_subelements; j++) {
        T8_ASSERT (element_is_valid (c[j]));
      }
    }
#endif

    /* Setting the parameter values for different subelements. 
   * The different subelement types (up to rotation) are:
   *                               
   *      x - - - - - - x         x - - - - - x        x - - - - - x        x - - - - - x        x - - x - - x        x - - x - - x
   *      |             |         | \   2   / |        | \       / |        | \       / |        | \   |   / |        | \   |   / |
   *      |             |         | 1 \   /   |        |   \   /   |        |   \   /   |        |   \ | /   |        |   \ | /   |
   *      |             |   -->   x - - X   3 |   or   x - - x     |   or   x - - x - - x   or   x - - x - - x   or   x - - x - - x
   *      |             |         | 0 /   \   |        |   / | \   |        |   /   \   |        |   /   \   |        |   / | \   |
   *      | elem        |         | /   4   \ |        | /   |   \ |        | /       \ |        | /       \ |        | /   |   \ |
   *      + - - - - - - x         x - - - - - x        x - - x - - x        x - - - - - x        x - - - - - x        x - - x - - x
   *           
   * Sub_ids are counted clockwise, starting with the (lower) left subelement with id 0.                    
   * Note, that we do not change the underlying quadrant. */

    for (int sub_id_counter = 0; sub_id_counter < num_subelements; sub_id_counter++) {

      standalone_scheme::element_copy (subelement_to_element (element),
                                       subelement_to_element (subelements[sub_id_counter]));
      subelements[sub_id_counter]->subelement_type = type;
      subelements[sub_id_counter]->subelement_id = sub_id_counter;
      T8_ASSERT (element_is_valid (c[sub_id_counter]));
    }
  }

 private:
  // PRIVATE HELPER
  static const t8_element_t *
  subelement_to_element (const t8_subelement_element *subelement) noexcept
  {
    return (const t8_element_t *) &subelement->element;
  }

  static t8_element_t *
  subelement_to_element (t8_subelement_element *subelement) noexcept
  {
    return (t8_element_t *) &subelement->element;
  }

  static const t8_element_t *
  element_to_element (const t8_element_t *element) noexcept
  {
    const t8_subelement_element *subelement = (const t8_subelement_element *) element;
    return subelement_to_element (subelement);
  }

  static t8_element_t *
  element_to_element (t8_element_t *element) noexcept
  {
    t8_subelement_element *subelement = (t8_subelement_element *) element;
    return subelement_to_element (subelement);
  }

  /** create the root element
   * \param [in,out] elem The element that is filled with the root
   */
  static void
  reset_subelement_values (t8_subelement_element *subelement) noexcept
  {
    subelement->subelement_type = 0;
    subelement->subelement_id = 0;
  }
};
