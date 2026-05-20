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
#include <t8_schemes/t8_standalone/t8_standalone.hxx>
#include <t8_schemes/t8_standalone/t8_standalone_implementation.hxx>
#include <t8_schemes/t8_subelement/t8_subelement_type.hxx>
#include <t8_schemes/t8_scheme_helpers.hxx>
#include <utility>
#include <algorithm>

/** A templated implementation of the scheme interface based on cutting planes. */
struct t8_subelementquad_scheme: public t8_scheme_helpers<T8_ECLASS_QUAD, t8_subelementquad_scheme>
{
 public:
  /** Constructor
  */
  t8_subelementquad_scheme () noexcept
    : m_element_size (sizeof (t8_subelement_element)), m_scheme_context (sc_mempool_new (m_element_size)) {};

 protected:
  // I am not sure why i even need this both variables.
  size_t m_element_size;  /**< The size in bytes of an element of class \a eclass */
  void *m_scheme_context; /**< Anonymous implementation context. */

 public:
  /** Destructor for all default schemes */
  ~t8_subelementquad_scheme ()
  {
    T8_ASSERT (m_scheme_context != NULL);
    SC_ASSERT (((sc_mempool_t *) m_scheme_context)->elem_count == 0);
    sc_mempool_destroy ((sc_mempool_t *) m_scheme_context);
  }

  /** Move constructor */
  t8_subelementquad_scheme (t8_subelementquad_scheme &&other) noexcept
    : m_element_size (other.m_element_size), m_scheme_context (std::exchange (other.m_scheme_context, nullptr))
  {
  }

  /** Move assignment operator */
  t8_subelementquad_scheme &
  operator= (t8_subelementquad_scheme &&other) noexcept
  {
    if (this != &other) {
      // Free existing resources of moved-to object
      if (m_scheme_context) {
        sc_mempool_destroy ((sc_mempool_t *) m_scheme_context);
      }

      // Transfer ownership of resources
      m_element_size = other.m_element_size;
      m_scheme_context = other.m_scheme_context;

      // Leave the source object in a valid state
      other.m_scheme_context = nullptr;
    }
    return *this;
  }

  /** Copy constructor */
  t8_subelementquad_scheme (const t8_subelementquad_scheme &other)
    : m_element_size (other.m_element_size), m_scheme_context (sc_mempool_new (other.m_element_size)) {};

  /** Copy assignment operator */
  t8_subelementquad_scheme &
  operator= (const t8_subelementquad_scheme &other)
  {
    if (this != &other) {
      // Free existing resources of assigned-to object
      if (m_scheme_context) {
        sc_mempool_destroy ((sc_mempool_t *) m_scheme_context);
      }

      // Copy the values from the source object
      m_element_size = other.m_element_size;
      m_scheme_context = sc_mempool_new (other.m_element_size);
    }
    return *this;
  }

  // ################################################____GENERAL INFO____################################################

  /** Return the size of any element of a given class.
   * \return                      The size of an element.
   */
  constexpr size_t
  get_element_size (void) const noexcept
  {
    return (sizeof (t8_subelement_element));
  }

  /** Returns true, if there is one element in the tree, that does not refine into 2^dim children.
   * Returns false otherwise.
   * \return                    non-zero if there is one element in the tree that does not refine into 2^dim children.
   */
  int
  refines_irregular (void) const noexcept
  {
    return true;  // Potentially there are subelements.
  }

  /** Return the maximum allowed level for any element of a given class.
   * \return                      The maximum allowed level for elements of class \b ts.
   */
  int
  get_maxlevel (void) const noexcept
  {
    return t8_standalone_scheme<T8_ECLASS_QUAD>::get_maxlevel ();
  }

  // ################################################____SHAPE INFORMATION____################################################

  /** Compute the number of corners of a given element.
   * \param [in] elem The element.
   * \return          The number of corners of \a elem.
   */
  int
  element_get_num_corners ([[maybe_unused]] const t8_element_t *elem) const noexcept
  {
    const t8_subelement_element *subelement = (const t8_subelement_element *) elem;
    T8_ASSERT (element_is_valid (elem));
    if (subelement->subelement_type == 0) {
      return t8_standalone_scheme<T8_ECLASS_QUAD>::element_get_num_corners (
        (const t8_element_t *) &subelement->element);
    }
    SC_ABORT ("This function is not implemented yet.\n");
  }

  /** Compute the number of faces of a given element.
   * \param [in] elem The element.
   * \return          The number of faces of \a elem.
   */
  int
  element_get_num_faces ([[maybe_unused]] const t8_element_t *elem) const noexcept
  {
    const t8_subelement_element *subelement = (const t8_subelement_element *) elem;
    T8_ASSERT (element_is_valid (elem));
    if (subelement->subelement_type == 0) {
      return t8_standalone_scheme<T8_ECLASS_QUAD>::element_get_num_faces ((const t8_element_t *) &subelement->element);
    }
    return T8_SUBELEMENT_FACES;
  }

  /** Compute the maximum number of faces of a given element and all of its
   *  descendants.
   * \param [in] elem The element.
   * \return          The maximum number of faces of \a elem and its descendants.
   */
  int
  element_get_max_num_faces ([[maybe_unused]] const t8_element_t *elem) const noexcept
  {
    const t8_subelement_element *subelement = (const t8_subelement_element *) elem;

    T8_ASSERT (element_is_valid (elem));

    if (subelement->subelement_type != 0) {
      return T8_SUBELEMENT_FACES;
    }
    return t8_standalone_scheme<T8_ECLASS_QUAD>::element_get_max_num_faces (
      (const t8_element_t *) &subelement->element);
  }

  /** Return the shape of an allocated element according its type.
   * For example, a child of an element can be an element of a different shape
   * and has to be handled differently - according to its shape.
   * \param [in] elem     The element to be considered
   * \return              The shape of the element as an eclass
   */
  t8_element_shape_t
  element_get_shape ([[maybe_unused]] const t8_element_t *elem) const noexcept
  {
    const t8_subelement_element *subelement = (const t8_subelement_element *) elem;

    T8_ASSERT (element_is_valid (elem));

    if (subelement->subelement_type == 0) {
      return t8_standalone_scheme<T8_ECLASS_QUAD>::element_get_shape ((const t8_element_t *) &subelement->element);
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
  int
  element_get_face_corner ([[maybe_unused]] const t8_element_t *element, const int face,
                           const int corner) const noexcept
  {
    const t8_subelement_element *subelement = (const t8_subelement_element *) element;

    T8_ASSERT (element_is_valid (element));

    if (subelement->subelement_type == 0) {
      return t8_standalone_scheme<T8_ECLASS_QUAD>::element_get_face_corner ((const t8_element_t *) &subelement->element,
                                                                            face, corner);
    }
    int t8_face_corners_subelement[3][2] = { { 0, 1 }, { 1, 2 }, { 2, 0 } };
    /*
     *
     *         x - - - - - x 1
     *         | \    f0 / |
     *         |   \ 0 /   |
     *         x - - x  el | f1
     *         |   /   \   |
     *         | /    f2 \ |
     *         x - - x - - x 2
     *               
     * The vertices of a subelement are enumerated clockwise, starting with the center vertex of the transition cell. 
     */

    T8_ASSERT (0 <= face && face < T8_SUBELEMENT_FACES);
    T8_ASSERT (0 <= corner && corner < 3);

    return t8_face_corners_subelement[face][corner];
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
  int
  element_get_corner_face (const t8_element_t *element, const int corner, const int face) const noexcept
  {
    const t8_subelement_element *subelement = (const t8_subelement_element *) element;
    T8_ASSERT (element_is_valid (element));
    if (subelement->subelement_type == 0) {
      return t8_standalone_scheme<T8_ECLASS_QUAD>::element_get_corner_face ((const t8_element_t *) &subelement->element,
                                                                            corner, face);
    }
    SC_ABORT ("This function is not implemented yet.\n");
  }

  /** Compute the shape of the face of an element.
   * \param [in] elem     The element.
   * \param [in] face     A face of \a elem.
   * \return              The element shape of the face.
   * I.e. T8_ECLASS_LINE for quads, T8_ECLASS_TRIANGLE for tets
   *      and depending on the face number either T8_ECLASS_QUAD or
   *      T8_ECLASS_TRIANGLE for prisms.
   */
  t8_element_shape_t
  element_get_face_shape ([[maybe_unused]] const t8_element_t *elem, [[maybe_unused]] const int face) const noexcept
  {
    T8_ASSERT (element_is_valid (elem));
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
  void
  element_copy (const t8_element_t *source, t8_element_t *dest) const noexcept
  {
    SC_ABORT ("This function is not implemented yet.\n");
  }

  /** Check if two elements are equal.
  * \param [in] elem1  The first element.
  * \param [in] elem2  The second element.
  * \return            true if the elements are equal, false if they are not equal
  */
  int
  element_is_equal (const t8_element_t *elem1, const t8_element_t *elem2) const noexcept
  {
    T8_ASSERT (element_is_valid (elem1));
    T8_ASSERT (element_is_valid (elem2));

    const t8_subelement_element *el1 = (const t8_subelement_element *) elem1;
    const t8_subelement_element *el2 = (const t8_subelement_element *) elem2;
    if (!t8_standalone_scheme<T8_ECLASS_QUAD>::element_is_equal ((const t8_element_t *) &el1->element,
                                                                 (const t8_element_t *) &el2->element))
      return 0;
    if (el1->subelement_type != el2->subelement_type) {
      return 0;
    }
    if (el1->subelement_id != el2->subelement_id) {
      return 0;
    }
    return 1;
  }

  // ################################################____ACCESSOR____################################################

  /** Return the level of a particular element.
   * \param [in] elem    The element whose level should be returned.
   * \return             The level of \b elem.
   */
  int
  element_get_level ([[maybe_unused]] const t8_element_t *elem) const noexcept
  {
    T8_ASSERT (element_is_valid (elem));
    return ((t8_subelement_element *) elem)->element.level;
  }

  // ################################################____REFINEMENT____################################################

  /** create the root element
   * \param [in,out] elem The element that is filled with the root
   */
  void
  set_to_root (t8_element_t *elem) const noexcept
  {
    t8_subelement_element *el = (t8_subelement_element *) elem;
    t8_standalone_scheme<T8_ECLASS_QUAD>::set_to_root ((t8_element_t *) &el->element);
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
  static void
  element_get_parent ([[maybe_unused]] const t8_element_t *elem, [[maybe_unused]] t8_element_t *parent) noexcept
  {
    SC_ABORT ("This function is not implemented yet.\n");
  }

  /** Compute the number of siblings of an element. That is the number of
   * elements with the same parent (if available).
   * \param [in] elem The element.
   * \return          The number of siblings of \a element.
   * Note that this number is >= 1, since we count the element itself as a sibling.
   * Note that the number of siblings is 1 for the root element.
   */
  int
  element_get_num_siblings (const t8_element_t *elem) const noexcept
  {
    const t8_subelement_element *subelement = (const t8_subelement_element *) elem;
    T8_ASSERT (element_is_valid (elem));
    if (subelement->subelement_type == 0) {
      return t8_standalone_scheme<T8_ECLASS_QUAD>::element_get_num_siblings (
        (const t8_element_t *) &subelement->element);
    }
    SC_ABORT ("This function is not implemented yet.\n");
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
  static void
  element_get_child ([[maybe_unused]] const t8_element_t *elem, [[maybe_unused]] const int childid,
                     [[maybe_unused]] t8_element_t *child) noexcept
  {
    SC_ABORT ("This function is not implemented yet.\n");
  }

  /** Return the number of children of an element when it is refined.
   * \param [in] elem   The element whose number of children is returned.
   * \return            The number of children of \a elem if it is to be refined.
   */
  static int
  element_get_num_children ([[maybe_unused]] const t8_element_t *elem) noexcept
  {
    SC_ABORT ("This function is not implemented yet.\n");
  }

  /** Return the max number of children of an eclass.
   * \return            The max number of children of \a element.
   */
  static int
  get_max_num_children () noexcept
  {
    SC_ABORT ("This function is not implemented yet.\n");
  }

  /**
   * Indicates if an element is refinable. Possible reasons for being not refinable could be
   * that the element has reached its max level.
   * \param [in] elem   The element to check.
   * \return            True if the element is refinable.
   */
  static bool
  element_is_refinable ([[maybe_unused]] const t8_element_t *elem) noexcept
  {
    SC_ABORT ("This function is not implemented yet.\n");
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
  element_get_children ([[maybe_unused]] const t8_element_t *elem, [[maybe_unused]] const int length,
                        [[maybe_unused]] t8_element_t *c[]) noexcept
  {
    SC_ABORT ("This function is not implemented yet.\n");
  }

  /** Compute the child id of an element.
   * \param [in] elem     This must be a valid element.
   * \return              The child id of elem.
   */
  static int
  element_get_child_id ([[maybe_unused]] const t8_element_t *elem) noexcept
  {
    SC_ABORT ("This function is not implemented yet.\n");
  }

  /** Compute the ancestor id of an element, that is the child id
   * at a given level.
   * \param [in] elem     This must be a valid element.
   * \param [in] level    A refinement level. Must satisfy \a level < elem.level
   * \return              The child_id of \a elem in regard to its \a level ancestor.
   */
  static int
  element_get_ancestor_id ([[maybe_unused]] const t8_element_t *elem,
                           [[maybe_unused]] const t8_element_level level) noexcept
  {
    SC_ABORT ("This function is not implemented yet.\n");
  }

  /** Query whether a given set of elements is a family or not.
   * \param [in] fam      An array of as many elements as an element of class
   *                      \b ts has siblings.
   * \return              Zero if \b fam is not a family, nonzero if it is.
   * \note level 0 elements do not form a family.
   */
  static int
  elements_are_family ([[maybe_unused]] t8_element_t *const *fam) noexcept
  {
    SC_ABORT ("This function is not implemented yet.\n");
  }

  // Note to devs: element_is_ancestor currently cannot be static
  //               since it uses the non-static function element_new
  /** Query whether element A is an ancestor of the element B.
   * An element A is ancestor of an element B if A == B or if B can 
   * be obtained from A via successive refinement.
   * \param [in] element_A An element of class \a eclass in scheme \a scheme.
   * \param [in] element_B An element of class \a eclass in scheme \a scheme.
   * \return     True if and only if \a element_A is an ancestor of \a element_B.
  */
  bool
  element_is_ancestor ([[maybe_unused]] const t8_element_t *element_A,
                       [[maybe_unused]] const t8_element_t *element_B) const noexcept
  {
    SC_ABORT ("This function is not implemented yet.\n");
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
  element_get_first_descendant ([[maybe_unused]] const t8_element_t *elem, [[maybe_unused]] t8_element_t *desc,
                                [[maybe_unused]] const t8_element_level level) noexcept
  {
    const t8_subelement_element *subelement = (const t8_subelement_element *) elem;
    t8_subelement_element *descsubelement = (t8_subelement_element *) desc;
    T8_ASSERT (element_is_valid (elem));
    if (subelement->subelement_type == 0) {
      t8_standalone_scheme<T8_ECLASS_QUAD>::element_get_first_descendant (
        (const t8_element_t *) &subelement->element, (t8_element_t *) &descsubelement->element, level);
      return;
    }
    SC_ABORT ("SUBELEMENTS: This function is not implemented yet.\n");

    // TODO: reset subelem values
  }

  /** Compute the last descendant of a given element.
   * \param [in] elem     The element whose descendant is computed.
   * \param [out] desc    The last element in a uniform refinement of \a elem
   *                      of the given level.
   * \param [in] level    The level, at which the descendant is computed.
   */
  static void
  element_get_last_descendant ([[maybe_unused]] const t8_element_t *elem, [[maybe_unused]] t8_element_t *desc,
                               [[maybe_unused]] const t8_element_level level) noexcept
  {
    const t8_subelement_element *subelement = (const t8_subelement_element *) elem;
    t8_subelement_element *descsubelement = (t8_subelement_element *) desc;
    T8_ASSERT (element_is_valid (elem));
    if (subelement->subelement_type == 0) {
      t8_standalone_scheme<T8_ECLASS_QUAD>::element_get_last_descendant (
        (const t8_element_t *) &subelement->element, (t8_element_t *) &descsubelement->element, level);
      return;
    }
    SC_ABORT ("SUBELEMENTS: This function is not implemented yet.\n");
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
    SC_ABORT ("This function is not implemented yet.\n");
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
  element_get_children_at_face ([[maybe_unused]] const t8_element_t *elem, [[maybe_unused]] const int face,
                                [[maybe_unused]] t8_element_t *children[], [[maybe_unused]] const int num_children,
                                [[maybe_unused]] int *child_indices) noexcept
  {
    SC_ABORT ("This function is not implemented yet.\n");
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
  element_face_get_child_face ([[maybe_unused]] const t8_element_t *elem, [[maybe_unused]] const int face,
                               [[maybe_unused]] const int face_child) noexcept
  {
    SC_ABORT ("This function is not implemented yet.\n");
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
  element_face_get_parent_face ([[maybe_unused]] const t8_element_t *elem, [[maybe_unused]] const int face) noexcept
  {
    SC_ABORT ("This function is not implemented yet.\n");
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
  element_get_first_descendant_face ([[maybe_unused]] const t8_element_t *elem, [[maybe_unused]] const int face,
                                     [[maybe_unused]] t8_element_t *first_desc,
                                     [[maybe_unused]] const t8_element_level level) noexcept
  {
    SC_ABORT ("This function is not implemented yet.\n");
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
  element_get_last_descendant_face ([[maybe_unused]] const t8_element_t *elem, [[maybe_unused]] const int face,
                                    [[maybe_unused]] t8_element_t *last_desc,
                                    [[maybe_unused]] const t8_element_level level) noexcept
  {
    SC_ABORT ("This function is not implemented yet.\n");
  }

  // ################################################____FACE NEIGHBOR____################################################

  /** Compute whether a given element shares a given face with its root tree.
   * \param [in] elem     The input element.
   * \param [in] face     A face of \a elem.
   * \return              True if \a face is a subface of the element's root element.
   * \note You can compute the corresponding face number of the tree via \ref element_get_tree_face.
   */
  static int
  element_is_root_boundary ([[maybe_unused]] const t8_element_t *elem, [[maybe_unused]] const int face) noexcept
  {
    SC_ABORT ("This function is not implemented yet.\n");
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
  element_get_tree_face ([[maybe_unused]] const t8_element_t *elem, [[maybe_unused]] const int face) noexcept
  {
    SC_ABORT ("This function is not implemented yet.\n");
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
  element_get_face_neighbor_inside ([[maybe_unused]] const t8_element_t *elem, [[maybe_unused]] t8_element_t *neigh,
                                    [[maybe_unused]] const int face, [[maybe_unused]] int *neigh_face) noexcept
  {
    SC_ABORT ("This function is not implemented yet.\n");
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
  static inline void
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
    SC_ABORT ("This function is not implemented yet.\n");
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
    t8_subelement_element *subelement = (t8_subelement_element *) elem;

    subelement->subelement_type = 0;
    subelement->subelement_id = 0;
    t8_standalone_scheme<T8_ECLASS_QUAD>::element_set_linear_id ((t8_element_t *) &subelement->element, level, id);
  }

  /** Compute the linear id of a given element in a hypothetical uniform
   * refinement of a given level.
   * \param [in] elem     The element whose id we compute.
   * \param [in] level    The level of the uniform refinement to consider.
   * \return              The linear id of the element.
   */
  static t8_linearidx_t
  element_get_linear_id ([[maybe_unused]] const t8_element_t *elem,
                         [[maybe_unused]] const t8_element_level level) noexcept
  {
    SC_ABORT ("This function is not implemented yet.\n");
  }

  /** Construct the successor in a uniform refinement of a given element.
   * \param [in] elem1    The element whose successor should be constructed.
   * \param [in,out] elem2  The element whose entries will be set.
   */
  static void
  element_construct_successor ([[maybe_unused]] const t8_element_t *elem1,
                               [[maybe_unused]] t8_element_t *elem2) noexcept
  {
    const t8_subelement_element *subelement1 = (const t8_subelement_element *) elem1;
    T8_ASSERT (element_is_valid (elem1));
    if (subelement1->subelement_type == 0) {
      t8_subelement_element *subelement2 = (t8_subelement_element *) elem2;
      t8_standalone_scheme<T8_ECLASS_QUAD>::element_construct_successor ((const t8_element_t *) &subelement1->element,
                                                                         (t8_element_t *) &subelement2->element);
      return;
    }
    SC_ABORT ("SUBELEMENTS: This function is not implemented yet.\n");
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
  t8_gloidx_t
  element_count_leaves (const t8_element_t *elem, const t8_element_level level) const noexcept
  {
    const t8_subelement_element *subelement = (const t8_subelement_element *) elem;
    T8_ASSERT (element_is_valid (elem));
    if (subelement->subelement_type == 0) {
      return t8_standalone_scheme<T8_ECLASS_QUAD>::element_count_leaves ((const t8_element_t *) &subelement->element,
                                                                         level);
    }
    SC_ABORT ("SUBELEMENTS: This function is not implemented yet.\n");
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
    return t8_standalone_scheme<T8_ECLASS_QUAD>::count_leaves_from_root (level);
  }

  /** Compare two elements.
   * \param [in] elem1  The first element.
   * \param [in] elem2  The second element.
   * \return       negative if elem1 < elem2, zero if elem1 equals elem2
   *               and positive if elem1 > elem2.
   *  If elem2 is a copy of elem1 then the elements are equal.
   */
  static int
  element_compare ([[maybe_unused]] const t8_element_t *elem1, [[maybe_unused]] const t8_element_t *elem2) noexcept
  {
    SC_ABORT ("This function is not implemented yet.\n");
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
  element_get_vertex_reference_coords ([[maybe_unused]] const t8_element_t *elem, [[maybe_unused]] const int vertex,
                                       [[maybe_unused]] double coords[]) noexcept
  {
    SC_ABORT ("This function is not implemented yet.\n");
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
  element_get_reference_coords ([[maybe_unused]] const t8_element_t *elem, [[maybe_unused]] const double *ref_coords,
                                [[maybe_unused]] const size_t num_coords, [[maybe_unused]] double *out_coords) noexcept
  {
    SC_ABORT ("This function is not implemented yet.\n");
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
  static void
  element_new ([[maybe_unused]] const int length, [[maybe_unused]] t8_element_t **elems) noexcept
  {
    t8_subelement_element *subelements = (t8_subelement_element *) elems;
    for (int i = 0; i < length; i++) {
      subelements[i].subelement_type = 0;
      subelements[i].subelement_id = 0;
      t8_standalone_scheme<T8_ECLASS_QUAD>::element_init (1, (t8_element_t *) &subelements[i].element);
    }
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
  static inline void
  element_init (const int length, t8_element_t *elems) noexcept
  {
    t8_subelement_element *subelements = (t8_subelement_element *) elems;
    for (int i = 0; i < length; i++) {
      subelements[i].subelement_type = 0;
      subelements[i].subelement_id = 0;
      t8_standalone_scheme<T8_ECLASS_QUAD>::element_init (1, (t8_element_t *) &subelements[i].element);
    }
  }

  /** Deinitialize an array of allocated elements.
   * \param [in] length   The number of elements to be deinitialized.
   * \param [in,out] elems On input an array of \a length many allocated
   *                       and initialized elements, on output an array of
   *                       \a length many allocated, but not initialized elements.
   * \note Call this function if you called element_init on the element pointers.
   * \see element_init
   */
  static void
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
  element_destroy ([[maybe_unused]] const int length, [[maybe_unused]] t8_element_t **elems) const noexcept
  {
    SC_ABORT ("This function is not implemented yet.\n");
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
  element_is_valid ([[maybe_unused]] const t8_element_t *elem) noexcept
  {
    const t8_subelement_element *subelement = (const t8_subelement_element *) elem;
    int element_valid
      = t8_standalone_scheme<T8_ECLASS_QUAD>::element_is_valid ((const t8_element_t *) &subelement->element);
    if (subelement->subelement_type == 0) {
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
  element_debug_print ([[maybe_unused]] const t8_element_t *elem) noexcept
  {

    SC_ABORT ("This function is not implemented yet.\n");
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
    SC_ABORT ("This function is not implemented yet.\n");
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
  element_MPI_Pack ([[maybe_unused]] t8_element_t **const elements, [[maybe_unused]] const unsigned int count,
                    [[maybe_unused]] void *send_buffer, [[maybe_unused]] const int buffer_size,
                    [[maybe_unused]] int *position, [[maybe_unused]] sc_MPI_Comm comm) const noexcept

  {
    SC_ABORT ("This function is not implemented yet.\n");
  }

  /** Determine an upper bound for the size of the packed message of \a count elements
     * \param [in] count Number of elements to pack
     * \param [in] comm MPI Communicator
     * \param [out] pack_size upper bound on the message size
    */
  void
  element_MPI_Pack_size ([[maybe_unused]] const unsigned int count, [[maybe_unused]] sc_MPI_Comm comm,
                         [[maybe_unused]] int *pack_size) const noexcept
  {
    SC_ABORT ("This function is not implemented yet.\n");
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
  element_MPI_Unpack ([[maybe_unused]] void *recvbuf, [[maybe_unused]] const int buffer_size,
                      [[maybe_unused]] int *position, [[maybe_unused]] t8_element_t **elements,
                      [[maybe_unused]] const unsigned int count, [[maybe_unused]] sc_MPI_Comm comm) const noexcept
  {
    SC_ABORT ("This function is not implemented yet.\n");
  }
};
