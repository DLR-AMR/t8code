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

/** \file t8_multilevel_implementation.hxx
 * The multilevel scheme converts any other scheme into a multilevel scheme.
 * This file contains the implementation of the multilevel scheme.
 */

#ifndef T8_MULTILEVEL_IMPLEMENTATION_HXX
#define T8_MULTILEVEL_IMPLEMENTATION_HXX

#include <t8_element.h>
#include <t8_schemes/t8_standalone/t8_standalone_elements.hxx>
#include <utility>

#if T8_ENABLE_DEBUG
#include <limits>
#endif

/* Forward declaration of the scheme so we can use it as an argument in the eclass schemes function. */
class t8_scheme;

template <typename TUnderlyingElementType>
struct t8_multilevel_element
{
  TUnderlyingElementType linear_element;
  bool is_child_of_itself; /** If true, the actual hierarchical level of the element is level(linear_element) + 1. */
};

/**
 * Converts an eclass scheme \a TUnderlyingEclassScheme into a multilevel scheme.
 * The complete tree is saved depth first inside this SFC. Also needs the type of the
 * underlying elements \a TUnderlyingElementType.
 * \tparam TUnderlyingEclassScheme  The scheme which is converted into a multilevel scheme.
 * \tparam TUnderlyingElementType   The type of the underlying elements.
 */
template <class TUnderlyingEclassScheme, typename TUnderlyingElementType>
class t8_multilevel_scheme: private TUnderlyingEclassScheme {

 private:
  using multilevel_element = t8_multilevel_element<TUnderlyingElementType>;
  using linear_element = TUnderlyingElementType;

 public:
  /**
   * Initialize the multilevel scheme with the adapted element size of the underlying scheme.
   * Forwards all arguments to the constructor of the underlying scheme.
   * \tparam _args        The types of the arguments.
   * \param [in] args     The arguments for the underlying scheme.
   */
  template <typename... _args>
  t8_multilevel_scheme (_args &&...args) noexcept
    : TUnderlyingEclassScheme (std::forward<_args> (args)...), multilevel_element_size (sizeof (multilevel_element)),
      multilevel_scheme_pool (sc_mempool_new (multilevel_element_size))
  {
  }

 protected:
  size_t multilevel_element_size; /**< The size in bytes of an element of
                                       t8_multilevel_element<TUnderlyingElementType> */
  void *multilevel_scheme_pool;   /**< Memory pool for multilevel elements. */

 public:
  ~t8_multilevel_scheme ()
  {
    T8_ASSERT (multilevel_scheme_pool != NULL);
    SC_ASSERT (((sc_mempool_t *) multilevel_scheme_pool)->elem_count == 0);
    sc_mempool_destroy ((sc_mempool_t *) multilevel_scheme_pool);
  }

  /** Move constructor */
  t8_multilevel_scheme (t8_multilevel_scheme &&other) noexcept
    : TUnderlyingEclassScheme (std::move (other)), multilevel_element_size (other.multilevel_element_size),
      multilevel_scheme_pool (std::exchange (other.multilevel_scheme_pool, nullptr))
  {
  }

  /** Move assignment operator */
  t8_multilevel_scheme &
  operator= (t8_multilevel_scheme &&other) noexcept
  {
    if (this != &other) {
      // Free existing resources of moved-to object
      if (multilevel_scheme_pool) {
        sc_mempool_destroy ((sc_mempool_t *) multilevel_scheme_pool);
      }

      // Transfer ownership of resources
      multilevel_element_size = other.multilevel_element_size;
      multilevel_scheme_pool = other.multilevel_scheme_pool;

      // Leave the source object in a valid state
      other.multilevel_scheme_pool = nullptr;
    }
    TUnderlyingEclassScheme::operator= (std::move (other));
    return *this;
  }

  /** Copy constructor */
  t8_multilevel_scheme (const t8_multilevel_scheme &other)
    : TUnderlyingEclassScheme (other), multilevel_element_size (other.multilevel_element_size),
      multilevel_scheme_pool (sc_mempool_new (other.multilevel_element_size))
  {
  }

  /** Copy assignment operator */
  t8_multilevel_scheme &
  operator= (const t8_multilevel_scheme &other)
  {
    if (this != &other) {
      // Free existing resources of assigned-to object
      if (multilevel_scheme_pool) {
        sc_mempool_destroy ((sc_mempool_t *) multilevel_scheme_pool);
      }

      // Copy the values from the source object
      multilevel_element_size = other.multilevel_element_size;
      multilevel_scheme_pool = sc_mempool_new (other.multilevel_element_size);
    }
    TUnderlyingEclassScheme::operator= (other);
    return *this;
  }

 public:
  // ################################################____GENERAL INFO____################################################

  /** Return the tree class of this scheme.
   * \return The tree class of this scheme.
   */
  inline t8_eclass_t
  get_eclass (void) const
  {
    return TUnderlyingEclassScheme::get_eclass ();
  }

  /** Return the size of any element of a given class.
   * \return                      The size of an element.
   */
  inline size_t
  get_element_size (void) const
  {
    return multilevel_element_size;
  }

  /** Returns true, if there is one element in the tree, that does not refine into 2^dim children.
   * Returns false otherwise.
   * \return                    non-zero if there is one element in the tree that does not refine into 2^dim children.
   */
  inline int
  refines_irregular (void) const
  {
    return TUnderlyingEclassScheme::refines_irregular ();
  }

  /** Return the maximum allowed level for any element of a given class.
   * \return                      The maximum allowed level for elements of class \b ts.
   */
  inline int
  get_maxlevel (void) const
  {
    /* The maxlevel is limited by the size of the linear id datatype.
    So check if the level is valid. We use a double to calculate the max id, since it should
    always be lange enough. */
#if T8_ENABLE_DEBUG
    T8_ASSERTF (!TUnderlyingEclassScheme::refines_irregular (),
                "Multilevel scheme conversion currently does not work with irregular schemes.\n");
    const t8_element_level maxlevel = TUnderlyingEclassScheme::get_maxlevel ();
    double count = 0; /* Use double because it should hold bigger numbers than t8_linearidx_t */
    const int dim = t8_eclass_to_dimension[get_eclass ()];
    for (size_t i_level = 0; i_level <= maxlevel; ++i_level) {
      count += pow (2.0, dim * i_level * 1.0);
    }
    T8_ASSERTF (std::numeric_limits<t8_linearidx_t>::max () >= count,
                "t8_linearidx_t cannot hold enough elements for multilevel conversion.\n");
#endif
    return TUnderlyingEclassScheme::get_maxlevel ();
  }

  // ################################################____SHAPE INFORMATION____################################################

  /** Compute the number of corners of a given element.
   * \param [in] elem The element.
   * \return          The number of corners of \a elem.
   */
  inline int
  element_get_num_corners (const t8_element_t *elem) const
  {
    T8_ASSERT (element_is_valid (elem));
    const multilevel_element *elem_m = (const multilevel_element *) elem;
    return TUnderlyingEclassScheme::element_get_num_corners ((const t8_element_t *) &elem_m->linear_element);
  }

  /** Compute the number of faces of a given element.
   * \param [in] elem The element.
   * \return          The number of faces of \a elem.
   */
  inline int
  element_get_num_faces (const t8_element_t *elem) const
  {
    T8_ASSERT (element_is_valid (elem));
    const multilevel_element *elem_m = (const multilevel_element *) elem;
    return TUnderlyingEclassScheme::element_get_num_faces ((const t8_element_t *) &elem_m->linear_element);
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
    const multilevel_element *elem_m = (const multilevel_element *) elem;
    return TUnderlyingEclassScheme::element_get_max_num_faces ((const t8_element_t *) &elem_m->linear_element);
  }

  /** Return the shape of an allocated element according its type.
   * For example, a child of an element can be an element of a different shape
   * and has to be handled differently - according to its shape.
   * \param [in] elem     The element to be considered
   * \return              The shape of the element as an eclass
   */
  inline t8_element_shape_t
  element_get_shape (const t8_element_t *elem) const
  {
    T8_ASSERT (element_is_valid (elem));
    const multilevel_element *elem_m = (const multilevel_element *) elem;
    return TUnderlyingEclassScheme::element_get_shape ((const t8_element_t *) &elem_m->linear_element);
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
  inline int
  element_get_face_corner (const t8_element_t *elem, const int face, const int corner) const
  {
    T8_ASSERT (element_is_valid (elem));
    const multilevel_element *elem_m = (const multilevel_element *) elem;
    return TUnderlyingEclassScheme::element_get_face_corner ((const t8_element_t *) &elem_m->linear_element, face,
                                                             corner);
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
  inline int
  element_get_corner_face (const t8_element_t *elem, const int corner, const int face) const
  {
    T8_ASSERT (element_is_valid (elem));
    const multilevel_element *elem_m = (const multilevel_element *) elem;
    return TUnderlyingEclassScheme::element_get_corner_face ((const t8_element_t *) &elem_m->linear_element, corner,
                                                             face);
  }

  /** Compute the shape of the face of an element.
   * \param [in] elem     The element.
   * \param [in] face     A face of \a elem.
   * \return              The element shape of the face.
   * I.e. T8_ECLASS_LINE for quads, T8_ECLASS_TRIANGLE for tets
   *      and depending on the face number either T8_ECLASS_QUAD or
   *      T8_ECLASS_TRIANGLE for prisms.
   */
  inline t8_element_shape_t
  element_get_face_shape (const t8_element_t *elem, const int face) const
  {
    T8_ASSERT (element_is_valid (elem));
    const multilevel_element *elem_m = (const multilevel_element *) elem;
    return TUnderlyingEclassScheme::element_get_face_shape ((const t8_element_t *) &elem_m->linear_element, face);
  }

  // ################################################____GENERAL HELPER____################################################

  /** Copy all entries of \b source to \b dest. \b dest must be an existing
   *  element. No memory is allocated by this function.
   * \param [in] source The element whose entries will be copied to \b dest.
   * \param [in,out] dest This element's entries will be overwrite with the
   *                    entries of \b source.
   * \note \a source and \a dest may point to the same element.
   */
  inline void
  element_copy (const t8_element_t *source, t8_element_t *dest) const
  {
    T8_ASSERT (element_is_valid (source));
    T8_ASSERT (element_is_valid (dest));
    const multilevel_element *source_m = (const multilevel_element *) source;
    multilevel_element *dest_m = (multilevel_element *) dest;
    dest_m->is_child_of_itself = source_m->is_child_of_itself;
    TUnderlyingEclassScheme::element_copy ((const t8_element_t *) &source_m->linear_element,
                                           (t8_element_t *) &dest_m->linear_element);
  }

  /** Check if two elements are equal.
  * \param [in] elem1  The first element.
  * \param [in] elem2  The second element.
  * \return            true if the elements are equal, false if they are not equal
  */
  inline int
  element_is_equal (const t8_element_t *elem1, const t8_element_t *elem2) const
  {
    T8_ASSERT (element_is_valid (elem1));
    T8_ASSERT (element_is_valid (elem2));
    if (element_get_level (elem1) != element_get_level (elem2))
      return 0;
    const multilevel_element *elem_m1 = (const multilevel_element *) elem1;
    const multilevel_element *elem_m2 = (const multilevel_element *) elem2;
    return TUnderlyingEclassScheme::element_is_equal ((const t8_element_t *) &elem_m1->linear_element,
                                                      (const t8_element_t *) &elem_m2->linear_element);
  }

  // ################################################____ACCESSOR____################################################

  /** Return the level of a particular element.
   * \param [in] elem    The element whose level should be returned.
   * \return             The level of \b elem.
   */
  inline int
  element_get_level (const t8_element_t *elem) const
  {
    T8_ASSERT (element_is_valid (elem));
    const multilevel_element *elem_m = (const multilevel_element *) elem;
    return TUnderlyingEclassScheme::element_get_level ((const t8_element_t *) &elem_m->linear_element)
           + elem_m->is_child_of_itself;
  }

  // ################################################____REFINEMENT____################################################

  /** Create the root element
   * \param [in,out] elem The element that is filled with the root
   */
  inline void
  set_to_root (t8_element_t *elem) const
  {
    multilevel_element *elem_m = (multilevel_element *) elem;
    elem_m->is_child_of_itself = 0;
    return TUnderlyingEclassScheme::set_to_root ((t8_element_t *) &elem_m->linear_element);
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
    const multilevel_element *elem_m = (multilevel_element *) elem;
    multilevel_element *parent_m = (multilevel_element *) parent;
    if (elem_m->is_child_of_itself) {
      /* Elem is child of itself, so return itself. */
      TUnderlyingEclassScheme::element_copy ((const t8_element_t *) &elem_m->linear_element,
                                             (t8_element_t *) &parent_m->linear_element);
    }
    else {
      TUnderlyingEclassScheme::element_get_parent ((const t8_element_t *) &elem_m->linear_element,
                                                   (t8_element_t *) &parent_m->linear_element);
    }
    /* Parent is never child of itself. */
    parent_m->is_child_of_itself = 0;
  }

  /** Compute the number of siblings of an element. That is the number of
   * elements with the same parent (if available).
   * \param [in] elem The element.
   * \return          The number of siblings of \a element.
   * Note that this number is >= 1, since we count the element itself as a sibling.
   * Note that the number of siblings is 1 for the root element.
   */
  inline int
  element_get_num_siblings (const t8_element_t *elem) const
  {
    T8_ASSERT (element_is_valid (elem));
    if (element_get_level (elem) == 0)
      return 1;
    /* In this scheme every element except the root has one more sibling. */
    if (TUnderlyingEclassScheme::refines_irregular ()) {
      SC_ABORT ("Not implemented yet.\n");
    }
    linear_element root;
    TUnderlyingEclassScheme::element_init (1, (t8_element_t *) &root);
    TUnderlyingEclassScheme::set_to_root ((t8_element_t *) &root);
    const t8_child_id num_children = TUnderlyingEclassScheme::element_get_num_children ((t8_element_t *) &root);
    TUnderlyingEclassScheme::element_deinit (1, (t8_element_t *) &root);
    return num_children + 1;
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
  inline void
  element_get_sibling (const t8_element_t *elem, const int sibid, t8_element_t *sibling) const
  {
    T8_ASSERT (element_is_valid (elem));
    const multilevel_element *elem_m = (multilevel_element *) elem;
    multilevel_element *sibling_m = (multilevel_element *) sibling;
    if (sibid == 0) {
      /* The first sibling is the parent as child of itself. */
      sibling_m->is_child_of_itself = 1;
      return TUnderlyingEclassScheme::element_get_parent ((const t8_element_t *) &elem_m->linear_element,
                                                          (t8_element_t *) &sibling_m->linear_element);
    }
    else {
      /* All other siblings are shiftet up one id. */
      sibling_m->is_child_of_itself = 0;
      return TUnderlyingEclassScheme::element_get_sibling ((const t8_element_t *) &elem_m->linear_element, sibid - 1,
                                                           (t8_element_t *) &sibling_m->linear_element);
    }
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
  inline void
  element_get_child (const t8_element_t *elem, const int childid, t8_element_t *child) const
  {
    T8_ASSERT (element_is_valid (elem));
    T8_ASSERT (element_is_valid (child));
    const multilevel_element *elem_m = (const multilevel_element *) elem;
    T8_ASSERT (elem_m->is_child_of_itself == 0);
    multilevel_element *child_m = (multilevel_element *) child;
    if (childid == 0) {
      /* The first child is the element itself. */
      element_copy (elem, child);
      child_m->is_child_of_itself = 1;
    }
    else {
      /* The other children are the normal children shifted by one. */
      TUnderlyingEclassScheme::element_get_child ((const t8_element_t *) &elem_m->linear_element, childid - 1,
                                                  (t8_element_t *) &child_m->linear_element);
      child_m->is_child_of_itself = 0;
    }
  }

  /** Return the number of children of an element when it is refined.
   * \param [in] elem   The element whose number of children is returned.
   * \return            The number of children of \a elem if it is to be refined.
   */
  inline int
  element_get_num_children (const t8_element_t *elem) const
  {
    T8_ASSERT (element_is_valid (elem));
    const multilevel_element *elem_m = (const multilevel_element *) elem;
    if (elem_m->is_child_of_itself)
      /* If an element is child of itself it cannot be refined anymore. */
      return 1;
    /* Increase the number of children by one so that an element becomes child of itself. */
    return 1 + TUnderlyingEclassScheme::element_get_num_children ((const t8_element_t *) &elem_m->linear_element);
  }

  /** Return the max number of children of an eclass.
   * \return            The max number of children of \a element.
   */
  inline int
  get_max_num_children () const
  {
    return 1 + TUnderlyingEclassScheme::get_max_num_children ();
  }

  /**
   * Indicates if an element is refinable. Possible reasons for being not refinable could be
   * that the element has reached its max level or that it is a subelement.
   * \param [in] elem   The element to check.
   * \return            True if the element is refinable.
   */
  inline bool
  element_is_refinable (const t8_element_t *elem) const
  {
    T8_ASSERT (element_is_valid (elem));
    const multilevel_element *elem_m = (const multilevel_element *) elem;
    /* An element which is child of itself cannot be refined. */
    if (elem_m->is_child_of_itself)
      return false;
    if (element_get_level (elem) >= get_maxlevel ())
      return false;
    return TUnderlyingEclassScheme::element_is_refinable ((const t8_element_t *) &elem_m->linear_element);
  }

  /** Construct all children of a given element.
   * \param [in] elem           This must be a valid element, bigger than maxlevel.
   * \param [in] length         The length of the output array \a c must match
   *                            the number of children.
   * \param [in,out] children   The storage for these \a length elements must exist
   *                            and match the element class in the children's ordering.
   *                            On output, all children are valid.
   * It is valid to call this function with elem = c[0].
   * \see t8_element_num_children
   * \see t8_element_child_eclass
   */
  inline void
  element_get_children (const t8_element_t *elem, const int length, t8_element_t *children[]) const
  {
    T8_ASSERT (element_is_valid (elem));
    T8_ASSERT (length == element_get_num_children (elem));
    multilevel_element *elem_m = (multilevel_element *) elem;
    T8_ASSERT (elem_m->is_child_of_itself == 0);
    multilevel_element **children_m = (multilevel_element **) children;

    /* The first child is the element itself. */
    T8_ASSERT (element_is_valid (*children));
    element_copy (elem, *children);
    children_m[0]->is_child_of_itself = 1;

    /* The rest are the normal children. */
    T8_ASSERT (TUnderlyingEclassScheme::element_get_num_children ((const t8_element_t *) &elem_m->linear_element)
               == length - 1);
    for (t8_child_id child_id = 1; child_id < length; ++child_id) {
      TUnderlyingEclassScheme::element_get_child ((const t8_element_t *) &elem_m->linear_element, child_id - 1,
                                                  (t8_element_t *) &children_m[child_id]->linear_element);
      children_m[child_id]->is_child_of_itself = 0;
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
    multilevel_element *elem_m = (multilevel_element *) elem;
    /* If the element is child if itself it has id 0. */
    if (elem_m->is_child_of_itself) {
      return 0;
    }
    /* All other children are shifted by one to make space for the first child. */
    return 1 + TUnderlyingEclassScheme::element_get_child_id ((const t8_element_t *) &elem_m->linear_element);
  }

  /** Compute the ancestor id of an element, that is the child id
   * at a given level.
   * \param [in] elem     This must be a valid element.
   * \param [in] level    A refinement level. Must satisfy \a level < elem.level
   * \return              The child_id of \a elem in regard to its \a level ancestor.
   */
  inline int
  element_get_ancestor_id (const t8_element_t *elem, const t8_element_level level) const
  {
    T8_ASSERT (element_is_valid (elem));
    multilevel_element *elem_m = (multilevel_element *) elem;
    const int elem_level = element_get_level (elem);
    T8_ASSERT (level <= elem_level);
    /* If the element is child of itself the id is always 0. */
    if (elem_m->is_child_of_itself) {
      return 0;
    }
    /* All other children are shifted by one to make space for the first child. */
    return 1 + TUnderlyingEclassScheme::element_get_ancestor_id ((const t8_element_t *) &elem_m->linear_element, level);
  }

  /** Query whether a given set of elements is a family or not.
   * \param [in] fam      An array of as many elements as an element of class
   *                      \b ts has siblings.
   * \return              Zero if \b fam is not a family, nonzero if it is.
   * \note level 0 elements do not form a family.
   */
  inline int
  elements_are_family (t8_element_t *const *fam) const
  {
#if T8_ENABLE_DEBUG
    const int num_siblings = element_get_num_siblings (fam[0]);
    for (int isib = 0; isib < num_siblings; isib++) {
      T8_ASSERT (element_is_valid (fam[isib]));
    }
#endif
    /* The root is a family */
    if (element_get_level (fam[0]) == 0)
      return 1;

    multilevel_element *const *fam_m = (multilevel_element *const *) fam;

    /* The first element should be parent of the other elements. */
    linear_element parent;
    TUnderlyingEclassScheme::element_init (1, (t8_element_t *) &parent);
    TUnderlyingEclassScheme::element_get_parent ((const t8_element_t *) &fam_m[1]->linear_element,
                                                 (t8_element_t *) &parent);
    const bool is_parent = TUnderlyingEclassScheme::element_is_equal ((const t8_element_t *) &fam_m[0]->linear_element,
                                                                      (t8_element_t *) &parent);
    TUnderlyingEclassScheme::element_deinit (1, (t8_element_t *) &parent);
    if (!is_parent)
      return 0;

    /* The other elements should be siblings. */
    const t8_child_id num_underlying_siblings
      = TUnderlyingEclassScheme::element_get_num_siblings ((const t8_element_t *) &fam_m[1]->linear_element);
    t8_element_t **siblings = T8_ALLOC (t8_element_t *, num_underlying_siblings);
    TUnderlyingEclassScheme::element_new (num_underlying_siblings, siblings);
    for (t8_child_id i_sibling = 0; i_sibling < num_underlying_siblings; ++i_sibling) {
      TUnderlyingEclassScheme::element_copy ((const t8_element_t *) &fam_m[i_sibling + 1]->linear_element,
                                             siblings[i_sibling]);
    }
    const bool is_fam = TUnderlyingEclassScheme::elements_are_family (siblings);
    TUnderlyingEclassScheme::element_destroy (num_underlying_siblings, siblings);
    T8_FREE (siblings);
    return is_fam;
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
  inline void
  element_get_nca (const t8_element_t *elem1, const t8_element_t *elem2, t8_element_t *nca) const
  {
    T8_ASSERT (element_is_valid (elem1));
    T8_ASSERT (element_is_valid (elem2));
    T8_ASSERT (element_is_valid (nca));
    const multilevel_element *elem_m1 = (multilevel_element *) elem1;
    const multilevel_element *elem_m2 = (multilevel_element *) elem2;
    multilevel_element *nca_m = (multilevel_element *) nca;
    /* The nca does not change due to the tree being multilevelized */
    TUnderlyingEclassScheme::element_get_nca ((const t8_element_t *) &elem_m1->linear_element,
                                              (const t8_element_t *) &elem_m2->linear_element,
                                              (t8_element_t *) &nca_m->linear_element);
    /* The ancestor cannot be child of itself since then it cannot get any ancestors. */
    nca_m->is_child_of_itself = 0;
  }

  /** Compute the first descendant of a given element.
   * \param [in] elem     The element whose descendant is computed.
   * \param [out] desc    The first element in a uniform refinement of \a elem
   *                      of the given level.
   * \param [in] level    The level, at which the descendant is computed.
   */
  inline void
  element_get_first_descendant (const t8_element_t *elem, t8_element_t *desc, const t8_element_level level) const
  {
    multilevel_element *desc_m = (multilevel_element *) desc;
    T8_ASSERT (element_is_valid (elem));
    const t8_element_level elem_level = element_get_level (elem);
    T8_ASSERT (level >= elem_level);

    /* The first descendant is the first descendant one level above as child of itself. */
    element_copy (elem, desc);

    /* The element is only child of itself if the requested level is higher than the elem level. */
    if (level > elem_level)
      desc_m->is_child_of_itself = 1;
    else
      desc_m->is_child_of_itself = 0;
  }

  /** Compute the last descendant of a given element.
   * \param [in] elem     The element whose descendant is computed.
   * \param [out] desc    The last element in a uniform refinement of \a elem
   *                      of the given level.
   * \param [in] level    The level, at which the descendant is computed.
   */
  inline void
  element_get_last_descendant (const t8_element_t *elem, t8_element_t *desc, const t8_element_level level) const
  {
    T8_ASSERT (element_is_valid (elem));
    const multilevel_element *elem_m = (const multilevel_element *) elem;
    multilevel_element *desc_m = (multilevel_element *) desc;
    const t8_element_level elem_level = element_get_level (elem);
    T8_ASSERT (level >= elem_level);

    if (elem_m->is_child_of_itself) {
      element_copy (elem, desc);
      return;
    }
    /* The last descendant is given by the underlying scheme. */
    TUnderlyingEclassScheme::element_get_last_descendant ((const t8_element_t *) &elem_m->linear_element,
                                                          (t8_element_t *) &desc_m->linear_element, level);
    desc_m->is_child_of_itself = 0;
  }

  // ################################################____FACE REFINEMENT____################################################

  /** Return the number of children of an element's face when the element is refined.
   * \param [in] elem   The element whose face is considered.
   * \param [in] face   A face of \a elem.
   * \return            The number of children of \a face if \a elem is to be refined.
   */
  inline int
  element_get_num_face_children ([[maybe_unused]] const t8_element_t *elem, [[maybe_unused]] const int face) const
  {
    SC_ABORTF ("Not implemented.");
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
  inline void
  element_get_children_at_face ([[maybe_unused]] const t8_element_t *elem, [[maybe_unused]] const int face,
                                [[maybe_unused]] t8_element_t *children[], [[maybe_unused]] const int num_children,
                                [[maybe_unused]] int *child_indices) const
  {
    SC_ABORTF ("Not implemented.");
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
  inline int
  element_face_get_child_face ([[maybe_unused]] const t8_element_t *elem, [[maybe_unused]] const int face,
                               [[maybe_unused]] const int face_child) const
  {
    SC_ABORTF ("Not implemented.");
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
  inline int
  element_face_get_parent_face ([[maybe_unused]] const t8_element_t *elem, [[maybe_unused]] const int face) const
  {
    SC_ABORTF ("Not implemented.");
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
  inline void
  element_get_first_descendant_face ([[maybe_unused]] const t8_element_t *elem, [[maybe_unused]] const int face,
                                     [[maybe_unused]] t8_element_t *first_desc,
                                     [[maybe_unused]] const t8_element_level level) const
  {
    SC_ABORTF ("Not implemented.");
  }

  /** Construct the last descendant of an element at a given level that touches a given face.
   * \param [in] elem      The input element.
   * \param [in] face      A face of \a elem.
   * \param [in, out] last_desc An allocated element. This element's data will be
   *                       filled with the data of the last descendant of \a elem
   *                       that shares a face with \a face.
   * \param [in] level     The level, at which the last descendant is constructed
   */
  inline void
  element_get_last_descendant_face ([[maybe_unused]] const t8_element_t *elem, [[maybe_unused]] const int face,
                                    [[maybe_unused]] t8_element_t *last_desc,
                                    [[maybe_unused]] const t8_element_level level) const
  {
    SC_ABORTF ("Not implemented.");
  }

  // ################################################____FACE NEIGHBOR____################################################

  /** Compute whether a given element shares a given face with its root tree.
   * \param [in] elem     The input element.
   * \param [in] face     A face of \a elem.
   * \return              True if \a face is a subface of the element's root element.
   * \note You can compute the corresponding face number of the tree via \ref t8_element_tree_face.
   */
  inline int
  element_is_root_boundary (const t8_element_t *elem, const int face) const
  {
    T8_ASSERT (element_is_valid (elem));
    const multilevel_element *elem_m = (const multilevel_element *) elem;
    return TUnderlyingEclassScheme::element_is_root_boundary ((const t8_element_t *) &elem_m->linear_element, face);
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
  inline int
  element_get_tree_face (const t8_element_t *elem, const int face) const
  {
    T8_ASSERT (element_is_valid (elem));
    const multilevel_element *elem_m = (const multilevel_element *) elem;
    return TUnderlyingEclassScheme::element_get_tree_face ((t8_element_t *) &elem_m->linear_element, face);
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
  inline int
  element_get_face_neighbor_inside ([[maybe_unused]] const t8_element_t *elem, [[maybe_unused]] t8_element_t *neigh,
                                    [[maybe_unused]] const int face, [[maybe_unused]] int *neigh_face) const
  {
    SC_ABORTF ("Not implemented.");
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
  inline void
  element_transform_face ([[maybe_unused]] const t8_element_t *elem1, [[maybe_unused]] t8_element_t *elem2,
                          [[maybe_unused]] const int orientation, [[maybe_unused]] const int sign,
                          [[maybe_unused]] const int is_smaller_face) const
  {
    SC_ABORTF ("Not implemented.");
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
  inline int
  element_extrude_face ([[maybe_unused]] const t8_element_t *face, [[maybe_unused]] t8_element_t *elem,
                        [[maybe_unused]] const int root_face, [[maybe_unused]] const t8_scheme *scheme) const
  {
    SC_ABORTF ("Not implemented.");
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
  inline void
  element_get_boundary_face ([[maybe_unused]] const t8_element_t *elem, [[maybe_unused]] const int face,
                             [[maybe_unused]] t8_element_t *boundary, [[maybe_unused]] const t8_scheme *scheme) const
  {
    SC_ABORTF ("Not implemented.");
  }

  // ################################################____LINEAR ID____################################################

  /** Initialize the entries of an allocated element according to a
   *  given linear id in a uniform refinement.
   * \param [in,out] elem The element whose entries will be set.
   * \param [in] level    The level of the uniform refinement to consider.
   * \param [in] id       The linear id.
   *                      id must fulfil 0 <= id < 'number of leaves in the uniform refinement'
   */
  inline void
  element_set_linear_id (t8_element_t *elem, const t8_element_level uniform_level, t8_linearidx_t id) const
  {
    multilevel_element *elem_m = (multilevel_element *) elem;
    const uint8_t dim = t8_eclass_to_dimension[TUnderlyingEclassScheme::get_eclass ()];
    const t8_element_level maxlvl = get_maxlevel ();
#if T8_ENABLE_DEBUG
    const t8_linearidx_t id_max = get_num_elem_in_regular_subtree (dim, maxlvl);
    T8_ASSERT (id < id_max);
#endif
    int level = 0;                 /* current operating level */
    t8_linearidx_t id_linear = id; /* linear id */
    int id_in_subtree = id_linear; /* id in subtree */
    t8_linearidx_t subtree_id;     /* id of the subtree */
    for (; level < maxlvl; ++level) {
      /* if id in subtree is 0 this is the root */
      if (id_in_subtree == 0) {
        break;
      }
      /* Subtract the root of the current subtree */
      id_linear--;
      /* Subtract the subtrees before */
      if (level < maxlvl - 1) {
        /* compute the subtree id */
        subtree_id = (id_in_subtree - 1) / get_num_elem_in_regular_subtree (dim, maxlvl - level - 1);
        /* compute next id in subtree. For this we subtract all the subtrees before and the root of the current subtree */
        id_in_subtree -= subtree_id * get_num_elem_in_regular_subtree (dim, maxlvl - level - 1) + 1;
        /* subtract the multilevel elements of the subtrees before of the current subtree from the linear id */
        id_linear -= subtree_id * get_num_elem_in_regular_subtree (dim, maxlvl - level - 2);
      }
    }
    T8_ASSERT (level <= uniform_level);
    elem_m->is_child_of_itself = level < uniform_level;
    TUnderlyingEclassScheme::element_set_linear_id ((t8_element_t *) &elem_m->linear_element, level, id_linear);
  }

  /** Compute the linear id of a given element in a hypothetical uniform
   * refinement of a given level.
   * \param [in] elem     The element whose id we compute.
   * \param [in] level    The level of the uniform refinement to consider.
   * \return              The linear id of the element.
   */
  inline t8_linearidx_t
  element_get_linear_id (const t8_element_t *elem, const t8_element_level level) const
  {
    T8_ASSERT (element_is_valid (elem));
    const int maxlevel = get_maxlevel ();
    T8_ASSERT (0 <= level && level <= maxlevel);
    const multilevel_element *elem_m = (const multilevel_element *) elem;
    const int underlying_element_level
      = TUnderlyingEclassScheme::element_get_level ((t8_element_t *) &elem_m->linear_element);
    if (level < underlying_element_level) {
      /* The level of the element is higher than the requested level.
         coarsen element until requested level is reached. */
      multilevel_element parent;
      element_init (1, (t8_element_t *) &parent);
      element_get_parent (elem, (t8_element_t *) &parent);
      return element_get_linear_id ((const t8_element_t *) &parent, level);
    }

    const int dim = t8_eclass_to_dimension[get_eclass ()];
    const t8_linearidx_t id_linear
      = TUnderlyingEclassScheme::element_get_linear_id ((const t8_element_t *) &elem_m->linear_element, maxlevel);

    /* The multilevel conversion happens via the following formula:
     * #\f$\mathrm{id_{multilevel}} (\mathrm{id_{linear}, lvl}) = \mathrm{lvl} + \sum_{n = 0}^{\mathrm{lvl_{max}-1}} \lfloor \mathrm{id_{linear}} / 2^{n \cdot d} \rfloor \f$
     */
    t8_linearidx_t id_multilevel = underlying_element_level;  //-elem_m->is_child_of_itself;
    for (int i_level = 0; i_level < maxlevel; i_level++) {
      /* This is just id_multilevel += id_linear / sc_intpow (2, i_level * dim); */
      id_multilevel += id_linear >> (i_level * dim);
    }
    return id_multilevel;
  }

  /** Construct the successor in a uniform refinement of a given element.
   * \param [in] elem           The element whose successor should be constructed.
   * \param [in] uniform_level  The level of the uniform refinement.
   * \param [in,out] succ       The element whose entries will be set.
   */
  inline void
  element_construct_successor (const t8_element_t *elem, const t8_element_level uniform_level, t8_element_t *succ) const
  {
    T8_ASSERT (element_is_valid (elem));
    multilevel_element *succ_m = (multilevel_element *) succ;

    if (element_get_level (elem) == 0) {
      element_get_child (elem, 0, succ);
      T8_ASSERT (element_is_valid (succ));
      return;
    }

    element_copy (elem, succ);

    const t8_child_id child_id = element_get_child_id (elem);
    const int num_siblings = element_get_num_siblings (elem);
    T8_ASSERT (0 <= child_id && child_id < num_siblings);
    /* If the element is the last child of the parent, we need to go to the parent's successor (go to a coarser level)*/
    if (child_id == num_siblings - 1) {
      element_get_parent (succ, succ);
      element_construct_successor (succ, uniform_level, succ);
    }
    else {
      element_get_parent (succ, succ);
      element_get_child (succ, child_id + 1, succ);
      if (!succ_m->is_child_of_itself && element_get_level (succ) < uniform_level) {
        element_get_child (succ, 0, succ);
      }
      T8_ASSERT (succ_m->is_child_of_itself || element_get_level (succ) == uniform_level);
    }
    T8_ASSERT (element_is_valid (succ));
  }

  /** Count how many leaf descendants of a given uniform level an element would produce.
   * \param [in] elem  The element to be checked.
   * \param [in] level A refinement level.
   * \return Suppose \a t is uniformly refined up to level \a level. The return value
   * is the resulting number of elements (of the given level).
   * If \a level < t8_element_level(t), the return value should be 0.
   *
   * Example: If \a t is a line element that refines into 2 line elements on each level,
   *  then the return value is max(0, 2^{\a level - level(\a t)}).
   *  Thus, if \a t's level is 0, and \a level = 3, the return value is 2^3 = 8.
   */
  inline t8_gloidx_t
  element_count_leaves (const t8_element_t *elem, const t8_element_level level) const
  {
    T8_ASSERT (element_is_valid (elem));
    const multilevel_element *elem_m = (const multilevel_element *) elem;
    const t8_element_level elem_level = element_get_level (elem);
    T8_ASSERT (level >= elem_level);
    if (elem_m->is_child_of_itself) {
      return 0;
    }
    if (TUnderlyingEclassScheme::refines_irregular ()) {
      SC_ABORT ("Not implemented yet.\n");
    }
    else {
      t8_gloidx_t count_leaves
        = TUnderlyingEclassScheme::element_count_leaves ((const t8_element_t *) &elem_m->linear_element, level);
      const t8_child_id num_children
        = TUnderlyingEclassScheme::element_get_num_children ((const t8_element_t *) &elem_m->linear_element);
      /* For each level this scheme adds num_children^ilevel leaves. */
      for (int ilevel = 0; ilevel < level - elem_level; ++ilevel) {
        count_leaves += sc_intpow (num_children, ilevel);
      }
      return count_leaves;
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
  inline t8_gloidx_t
  count_leaves_from_root (const t8_element_level level) const
  {
    if (TUnderlyingEclassScheme::refines_irregular ()) {
      SC_ABORT ("Not implemented yet.\n");
    }
    else {
      t8_gloidx_t count_leaves = TUnderlyingEclassScheme::count_leaves_from_root (level);
      /* Get number of children of an element. */
      multilevel_element root;
      set_to_root ((t8_element_t *) &root);
      const t8_child_id num_children
        = TUnderlyingEclassScheme::element_get_num_children ((const t8_element_t *) &root.linear_element);
      /* For each level this scheme adds num_children^ilevel leaves. */
      for (size_t ilevel = 0; ilevel < level; ++ilevel) {
        count_leaves += sc_intpow (num_children, ilevel);
      }
      return count_leaves;
    }
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
    else if (id1 == id2) {
      /* The linear ids are the same, the element with the smaller level is considered smaller */
      T8_ASSERT (element_get_level (elem1) != element_get_level (elem2) || element_is_equal (elem1, elem2));
      return element_get_level (elem1) - element_get_level (elem2);
    }
    return 0;
  }

  // ################################################____VISUALIZATION____################################################

  /** Compute the coordinates of a given element vertex inside a reference tree
   *  that is embedded into [0,1]^d (d = dimension).
   *   \param [in] elem      The element to be considered.
   *   \param [in] vertex The id of the vertex whose coordinates shall be computed.
   *   \param [out] coords An array of at least as many doubles as the element's dimension
   *                      whose entries will be filled with the coordinates of \a vertex.
   */
  inline void
  element_get_vertex_reference_coords (const t8_element_t *elem, const int vertex, double coords[]) const
  {
    T8_ASSERT (element_is_valid (elem));
    const multilevel_element *elem_m = (const multilevel_element *) elem;
    TUnderlyingEclassScheme::element_get_vertex_reference_coords ((const t8_element_t *) &elem_m->linear_element,
                                                                  vertex, coords);
  }

  /** Convert a point in the reference space of an element to a point in the
   *  reference space of the tree.
   *
   * \param [in] elem         The element.
   * \param [in] coords_input The coordinates of the point in the reference space of the element.
   * \param [in] user_data    User data.
   * \param [out] out_coords  The coordinates of the point in the reference space of the tree.
   */
  inline void
  element_get_reference_coords (const t8_element_t *elem, const double *ref_coords, const size_t num_coords,
                                double *out_coords) const
  {
    T8_ASSERT (element_is_valid (elem));
    const multilevel_element *elem_m = (const multilevel_element *) elem;
    TUnderlyingEclassScheme::element_get_reference_coords ((const t8_element_t *) &elem_m->linear_element, ref_coords,
                                                           num_coords, out_coords);
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
  inline void
  element_new (int length, t8_element_t **elem) const
  {
    /* allocate memory */
    T8_ASSERT (this->multilevel_scheme_pool != NULL);
    T8_ASSERT (0 <= length);
    T8_ASSERT (elem != NULL);
    const multilevel_element **elem_m = (const multilevel_element **) elem;

    for (int i_elem = 0; i_elem < length; ++i_elem) {
      elem[i_elem] = (t8_element_t *) sc_mempool_alloc ((sc_mempool_t *) this->multilevel_scheme_pool);
      /* Init element with underlying scheme. */
      TUnderlyingEclassScheme::element_init (1, (t8_element_t *) &elem_m[i_elem]->linear_element);
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
  inline void
  element_init (const int length, t8_element_t *elem) const
  {
    const multilevel_element *elem_m = (const multilevel_element *) elem;
    for (int i_elem = 0; i_elem < length; ++i_elem) {
      TUnderlyingEclassScheme::element_init (1, (t8_element_t *) &elem_m[i_elem].linear_element);
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
  inline void
  element_deinit (const int length, t8_element_t *elem) const
  {
    const multilevel_element *elem_m = (const multilevel_element *) elem;
    for (int i_elem = 0; i_elem < length; ++i_elem) {
      TUnderlyingEclassScheme::element_deinit (1, (t8_element_t *) &elem_m[i_elem].linear_element);
    }
  }

  /** Deallocate an array of elements.
   * \param [in] length   The number of elements in the array.
   * \param [in,out] elems On input an array of \b length many allocated
   *                      element pointers.
   *                      On output all these pointers will be freed.
   *                      \b elem itself will not be freed by this function.
   */
  inline void
  element_destroy (const int length, t8_element_t **elem) const
  {
    T8_ASSERT (this->multilevel_scheme_pool != NULL);
    T8_ASSERT (0 <= length);
    T8_ASSERT (elem != NULL);
    const multilevel_element **elem_m = (const multilevel_element **) elem;
    for (int i_elem = 0; i_elem < length; ++i_elem) {
      TUnderlyingEclassScheme::element_deinit (1, (t8_element_t *) &elem_m[i_elem]->linear_element);
      sc_mempool_free ((sc_mempool_t *) multilevel_scheme_pool, elem[i_elem]);
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
  inline int
  element_is_valid (const t8_element_t *elem) const
  {
    const multilevel_element *elem_m = (const multilevel_element *) elem;
    return TUnderlyingEclassScheme::element_is_valid ((const t8_element_t *) &elem_m->linear_element);
  }

  /**
   * Print a given element. For a example for a triangle print the coordinates
   * and the level of the triangle. This function is only available in the
   * debugging configuration.
   *
   * \param [in]        elem  The element to print
   */
  inline void
  element_debug_print (const t8_element_t *elem) const
  {
    const multilevel_element *elem_m = (const multilevel_element *) elem;
    t8_debugf ("is_child_of_itself: %i\n", elem_m->is_child_of_itself);
    TUnderlyingEclassScheme::element_debug_print ((const t8_element_t *) &elem_m->linear_element);
  }

  /**
   * Print a given element. For a example for a triangle print the coordinates
   * and the level of the triangle. This function is only available in the
   * debugging configuration.
   *
   * \param [in]        elem  The element to print
   */
  inline void
  element_to_string (const t8_element_t *elem, char *debug_string, const int string_size) const
  {
    T8_ASSERT (element_is_valid (elem));
    const multilevel_element *elem_m = (const multilevel_element *) elem;
    T8_ASSERT (debug_string != NULL);
    snprintf (debug_string, string_size, "child_of_itself: %i, ", elem_m->is_child_of_itself);
    TUnderlyingEclassScheme::element_to_string ((const t8_element_t *) &elem_m->linear_element, debug_string,
                                                string_size);
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
  inline void
  element_MPI_Pack ([[maybe_unused]] t8_element_t **const elements, [[maybe_unused]] const unsigned int count,
                    [[maybe_unused]] void *send_buffer, [[maybe_unused]] const int buffer_size,
                    [[maybe_unused]] int *position, [[maybe_unused]] sc_MPI_Comm comm) const

  {
    SC_ABORTF ("Not implemented.");
  }

  /** Determine an upper bound for the size of the packed message of \a count elements
   * \param [in] count Number of elements to pack
   * \param [in] comm MPI Communicator
   * \param [out] pack_size upper bound on the message size
  */
  inline void
  element_MPI_Pack_size ([[maybe_unused]] const unsigned int count, [[maybe_unused]] sc_MPI_Comm comm,
                         [[maybe_unused]] int *pack_size) const
  {
    SC_ABORTF ("Not implemented.");
  }

  /** Unpack multiple elements from contiguous memory that was received via MPI.
   * \param [in] recvbuf Buffer from which to unpack the elements
   * \param [in] buffer_size size of the buffer (in order to check that we don't access out of range)
   * \param [in, out] position the position of the first byte that is not already packed
   * \param [in] elements Array of initialised elements that is to be filled from the message
   * \param [in] count Number of elements to unpack
   * \param [in] comm MPI Communicator
  */
  inline void
  element_MPI_Unpack ([[maybe_unused]] void *recvbuf, [[maybe_unused]] const int buffer_size,
                      [[maybe_unused]] int *position, [[maybe_unused]] t8_element_t **elements,
                      [[maybe_unused]] const unsigned int count, [[maybe_unused]] sc_MPI_Comm comm) const
  {
    SC_ABORTF ("Not implemented.");
  }

 private:
  // ################################################____HELPER____################################################

  /**
   * Compute the amount of elements in a level \a level multilevel subtree with dimendion \a dim.
   * \param [in] dim    The dimension of the subtree.
   * \param [in] level  The level of the subtree.
   * \return            The number of elements in that subtree.
  */
  inline t8_linearidx_t
  get_num_elem_in_regular_subtree (const int dim, const t8_element_level level) const
  {
    t8_linearidx_t count = 0;
    for (size_t i_level = 0; i_level <= level; ++i_level) {
      // The following is "count += sc_intpow (2, dim * i_level);" in bitshift
      count |= 1ULL << i_level * dim;
    }
    return count;
  }
};

#endif /* !T8_MULTILEVEL_IMPLEMENTATION_HXX */
