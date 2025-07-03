/*
  This file is part of t8code.
  t8code is a C library to manage a collection (a forest) of multiple
  connected adaptive space-trees of general element classes in parallel.

  Copyright (C) 2024 the developers

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

/** \file t8_scheme.hxx
 * This file defines the interface of eclass scheme implementations. Further, it defines the t8_scheme class, which
 * holds one or more element schemes. It also relays the function calls to the specific schemes.
 */

#ifndef T8_SCHEME_HXX
#define T8_SCHEME_HXX

#include <variant>
#include <vector>
#include <t8_refcount.h>
#include <t8_eclass.h>
#include <t8_schemes/t8_default/t8_default.hxx>
#include <t8_schemes/t8_default/t8_default_vertex/t8_default_vertex.hxx>
#include <t8_schemes/t8_default/t8_default_line/t8_default_line.hxx>
#include <t8_schemes/t8_default/t8_default_quad/t8_default_quad.hxx>
#include <t8_schemes/t8_default/t8_default_hex/t8_default_hex.hxx>
#include <t8_schemes/t8_default/t8_default_tri/t8_default_tri.hxx>
#include <t8_schemes/t8_default/t8_default_tet/t8_default_tet.hxx>
#include <t8_schemes/t8_default/t8_default_prism/t8_default_prism.hxx>
#include <t8_schemes/t8_default/t8_default_pyramid/t8_default_pyramid.hxx>
#include <t8_schemes/t8_standalone/t8_standalone.hxx>
#include <t8_schemes/t8_standalone/t8_standalone_implementation.hxx>
#include <string>
#if T8_ENABLE_DEBUG
// Only needed for t8_debug_print_type
#include <typeinfo>

/**
 * Get the type of the template parameter as a string.
 * \tparam TType          The template parameter to get the type of.
 * \return std::string&   The type of the template parameter as a string.
 */
template <typename TType>
inline std::string &
t8_debug_print_type ()
{
  static std::string type_name = typeid (TType).name ();
  return type_name;
}
#endif  // T8_ENABLE_DEBUG

/** This class holds one or more element schemes.
 * It also relays the function calls to the specific schemes. */
class t8_scheme {
  friend class t8_scheme_builder;

 public:
  t8_scheme ()
  {
    t8_refcount_init (&rc);
  };

  ~t8_scheme ()
  {
    if (sc_refcount_is_active (&rc)) {
      T8_ASSERT (t8_refcount_is_last (&rc));
      t8_refcount_unref (&rc);
    }
    t8_debugf ("Deleted the scheme.\n");
  };

  /* clang-format off */
  
  /** Variant to hold an eclass scheme. */
  using scheme_var = std::variant<
                                /* Default schemes */
                                t8_default_scheme_vertex,
                                t8_default_scheme_line,
                                t8_default_scheme_quad,
                                t8_default_scheme_tri,
                                t8_default_scheme_hex,
                                t8_default_scheme_tet,
                                t8_default_scheme_prism,
                                t8_default_scheme_pyramid,
                                t8_standalone_scheme<T8_ECLASS_VERTEX>,
                                t8_standalone_scheme<T8_ECLASS_LINE>,
                                t8_standalone_scheme<T8_ECLASS_QUAD>,
                                t8_standalone_scheme<T8_ECLASS_HEX>
                                >;
  /* clang-format on */

  using scheme_container = std::vector<scheme_var>; /**< Container type for holding eclass schemes. */

 private:
  scheme_container eclass_schemes; /**< The container holding the eclass schemes. */
  mutable t8_refcount_t
    rc; /**< The reference count of the scheme. Mutable so that the class can be const and the ref counter is still mutable. TODO: Replace by shared_ptr when forest becomes a class. */

 public:
  /**
   * Increase the reference count of the scheme.
   */
  inline void
  ref () const
  {
    t8_refcount_ref (&rc);
  }

  /**
   * Decrease the reference count of the scheme.
   * If the reference count reaches zero, the scheme is deleted.
   * \return The remaining reference count. If 0 the scheme was deleted.
   */
  inline int
  unref () const
  {
    const int remaining = rc.refcount - 1;
    if (t8_refcount_unref (&rc)) {
      t8_debugf ("Deleting the scheme.\n");
      delete this;
    }
    return remaining;
  }

  /** Get the number of eclass schemes inside the scheme.
   * \return The number of eclass schemes.
   */
  inline size_t
  get_num_eclass_schemes () const
  {
    return eclass_schemes.size ();
  }

  /** Check if the scheme is of a specific type.
   * \tparam TEclass_Scheme The type of the scheme to check for.
   * \param [in] tree_class    The eclass of the current tree.
   * \return                    True if the scheme is of type \a TEclassScheme, false otherwise.
   */
  template <class TEclassScheme>
  inline bool
  check_eclass_scheme_type (const t8_eclass_t tree_class) const
  {
    return std::holds_alternative<TEclassScheme> (eclass_schemes[tree_class]);
  }

  /** Get the eclass an eclass scheme is valid for. \Note: This function should return the input value as long as the
   * eclass schemes are soreted correctly. In the future, the trees will access the schemes by a key and then this
   * function will make more sense.
   * \param [in] tree_class     The eclass of the current tree.
   * \return                    The valid tree class for the eclass scheme.
   */
  inline t8_eclass_t
  get_eclass_scheme_eclass (const t8_eclass_t tree_class) const
  {
    return std::visit ([&] (auto &&scheme) { return scheme.get_eclass (); }, eclass_schemes[tree_class]);
  }

  /** Return the size of any element of a given class.
   * \param [in] tree_class    The eclass of the current tree.
   * \return                      The size of an element of class \a tree_class.
   * We provide a default implementation of this routine that should suffice
   * for most use cases.
   */
  inline size_t
  get_element_size (const t8_eclass_t tree_class) const
  {
    return std::visit ([&] (auto &&scheme) { return scheme.get_element_size (); }, eclass_schemes[tree_class]);
  };

  /** Returns true, if there is one element in the tree, that does not refine into 2^dim children.
   * Returns false otherwise.
   * \param [in] tree_class    The eclass of the current tree.
   * \return                   true if there is one element in the tree that does not refine into 2^dim children.
   */
  inline bool
  refines_irregular (const t8_eclass_t tree_class) const
  {
    return std::visit ([&] (auto &&scheme) { return scheme.refines_irregular (); }, eclass_schemes[tree_class]);
  };

  /** Return the maximum allowed level for any element of a given class.
   * \param [in] tree_class    The eclass of the current tree.
   * \return                      The maximum allowed level for elements of class \a tree_class.
   */
  inline int
  get_maxlevel (const t8_eclass_t tree_class) const
  {
    return std::visit ([&] (auto &&scheme) { return scheme.get_maxlevel (); }, eclass_schemes[tree_class]);
  };

  /** Return the level of a particular element.
   * \param [in] tree_class    The eclass of the current tree.
   * \param [in] elem    The element whose level should be returned.
   * \return             The level of \a elem.
   */
  inline int
  element_get_level (const t8_eclass_t tree_class, const t8_element_t *elem) const
  {
    return std::visit ([&] (auto &&scheme) { return scheme.element_get_level (elem); }, eclass_schemes[tree_class]);
  };

  /** Copy all entries of \a source to \a dest. \a dest must be an existing
   *  element. No memory is allocated by this function.
   * \param [in] tree_class    The eclass of the current tree.
   * \param [in] source The element whose entries will be copied to \a dest.
   * \param [in,out] dest This element's entries will be overwritten with the
   *                    entries of \a source.
   * \note \a source and \a dest may point to the same element.
   */
  inline void
  element_copy (const t8_eclass_t tree_class, const t8_element_t *source, t8_element_t *dest) const
  {
    return std::visit ([&] (auto &&scheme) { return scheme.element_copy (source, dest); }, eclass_schemes[tree_class]);
  };

  /** Compare two elements.
   * \param [in] tree_class    The eclass of the current tree.
   * \param [in] elem1  The first element.
   * \param [in] elem2  The second element.
   * \return       negative if elem1 < elem2, zero if elem1 equals elem2
   *               and positive if elem1 > elem2.
   *  If elem2 is a copy of elem1 then the elements are equal.
   */
  inline int
  element_compare (const t8_eclass_t tree_class, const t8_element_t *elem1, const t8_element_t *elem2) const
  {
    return std::visit ([&] (auto &&scheme) { return scheme.element_compare (elem1, elem2); },
                       eclass_schemes[tree_class]);
  };

  /** Check if two elements are equal.
   * \param [in] tree_class    The eclass of the current tree.
  * \param [in] elem1  The first element.
  * \param [in] elem2  The second element.
  * \return            true if the elements are equal, false if they are not equal
  */
  inline bool
  element_is_equal (const t8_eclass_t tree_class, const t8_element_t *elem1, const t8_element_t *elem2) const
  {
    return std::visit ([&] (auto &&scheme) { return scheme.element_is_equal (elem1, elem2); },
                       eclass_schemes[tree_class]);
  };

  /**
   * Indicates if an element is refinable. Possible reasons for being not refinable could be
   * that the element has reached its max level.
   * \param [in] elem   The element to check.
   * \return            True if the element is refinable.
   */
  inline bool
  element_is_refinable (const t8_eclass_t tree_class, const t8_element_t *elem) const
  {
    return std::visit ([&] (auto &&scheme) { return scheme.element_is_refinable (elem); }, eclass_schemes[tree_class]);
  };

  /** Compute the parent of a given element \a elem and store it in \a parent.
   *  \a parent needs to be an existing element. No memory is allocated by this function.
   *  \a elem and \a parent can point to the same element, then the entries of
   *  \a elem are overwritten by the ones of its parent.
   * \param [in] tree_class    The eclass of the current tree.
   * \param [in] elem   The element whose parent will be computed.
   * \param [in,out] parent This element's entries will be overwritten by those
   *                    of \a elem's parent.
   *                    The storage for this element must exist
   *                    and match the element class of the parent.
   *                    For a pyramid, for example, it may be either a
   *                    tetrahedron or a pyramid depending on \a elem's childid.
   */
  inline void
  element_get_parent (const t8_eclass_t tree_class, const t8_element_t *elem, t8_element_t *parent) const
  {
    return std::visit ([&] (auto &&scheme) { return scheme.element_get_parent (elem, parent); },
                       eclass_schemes[tree_class]);
  };

  /** Compute the number of siblings of an element. That is the number of 
   * Children of its parent.
   * \param [in] tree_class    The eclass of the current tree.
   * \param [in] elem The element.
   * \return          The number of siblings of \a element.
   * Note that this number is >= 1, since we count the element itself as a sibling.
   */
  inline int
  element_get_num_siblings (const t8_eclass_t tree_class, const t8_element_t *elem) const
  {
    return std::visit ([&] (auto &&scheme) { return scheme.element_get_num_siblings (elem); },
                       eclass_schemes[tree_class]);
  };

  /** Compute a specific sibling of a given element \a elem and store it in \a sibling.
   *  \a sibling needs to be an existing element. No memory is allocated by this function.
   *  \a elem and \a sibling can point to the same element, then the entries of
   *  \a elem are overwritten by the ones of its sibid-th sibling.
   * \param [in] tree_class    The eclass of the current tree.
   * \param [in] elem   The element whose sibling will be computed.
   * \param [in] sibid  The id of the sibling computed.
   * \param [in,out] sibling This element's entries will be overwritten by those
   *                    of \a elem's sibid-th sibling.
   *                    The storage for this element must exist
   *                    and match the element class of the sibling.
   */
  inline void
  element_get_sibling (const t8_eclass_t tree_class, const t8_element_t *elem, const int sibid,
                       t8_element_t *sibling) const
  {
    return std::visit ([&] (auto &&scheme) { return scheme.element_get_sibling (elem, sibid, sibling); },
                       eclass_schemes[tree_class]);
  };

  /** Compute the number of corners of a given element.
   * \param [in] tree_class    The eclass of the current tree.
   * \param [in] elem The element.
   * \return          The number of corners of \a elem.
   */
  inline int
  element_get_num_corners (const t8_eclass_t tree_class, const t8_element_t *elem) const
  {
    return std::visit ([&] (auto &&scheme) { return scheme.element_get_num_corners (elem); },
                       eclass_schemes[tree_class]);
  };

  /** Compute the number of faces of a given element.
   * \param [in] tree_class    The eclass of the current tree.
   * \param [in] elem The element.
   * \return          The number of faces of \a elem.
   */
  inline int
  element_get_num_faces (const t8_eclass_t tree_class, const t8_element_t *elem) const
  {
    return std::visit ([&] (auto &&scheme) { return scheme.element_get_num_faces (elem); }, eclass_schemes[tree_class]);
  };

  /** Compute the maximum number of faces of a given element and all of its
   *  descendants.
   * \param [in] tree_class    The eclass of the current tree.
   * \param [in] elem The element.
   * \return          The maximum number of faces of \a elem and its descendants.
   */
  inline int
  element_get_max_num_faces (const t8_eclass_t tree_class, const t8_element_t *elem) const
  {
    return std::visit ([&] (auto &&scheme) { return scheme.element_get_max_num_faces (elem); },
                       eclass_schemes[tree_class]);
  };

  /** Return the number of children of an element when it is refined.
   * \param [in] tree_class    The eclass of the current tree.
   * \param [in] elem   The element whose number of children is returned.
   * \return            The number of children of \a elem if it is to be refined.
   */
  inline int
  element_get_num_children (const t8_eclass_t tree_class, const t8_element_t *elem) const
  {
    return std::visit ([&] (auto &&scheme) { return scheme.element_get_num_children (elem); },
                       eclass_schemes[tree_class]);
  };

  /** Return the max number of children of an eclass.
   * \param [in] tree_class    The eclass of tree the elements are part of.
   * \return            The max number of children of \a element.
   */
  inline int
  get_max_num_children (const t8_eclass_t tree_class) const
  {
    return std::visit ([&] (auto &&scheme) { return scheme.get_max_num_children (); }, eclass_schemes[tree_class]);
  };

  /** Return the number of children of an element's face when the element is refined.
   * \param [in] tree_class    The eclass of the current tree.
   * \param [in] elem   The element whose face is considered.
   * \param [in] face   A face of \a elem.
   * \return            The number of children of \a face if \a elem is to be refined.
   */
  inline int
  element_get_num_face_children (const t8_eclass_t tree_class, const t8_element_t *elem, const int face) const
  {
    return std::visit ([&] (auto &&scheme) { return scheme.element_get_num_face_children (elem, face); },
                       eclass_schemes[tree_class]);
  };

  /** Return the corner number of an element's face corner.
   * Example quad: 2 x --- x 3
   *                 |     |
   *                 |     |   face 1
   *               0 x --- x 1
   *      Thus for face = 1 the output is: corner=0 : 1, corner=1: 3
   *
   * \param [in] tree_class    The eclass of the current tree.
   * \param [in] element  The element.
   * \param [in] face     A face index for \a element.
   * \param [in] corner   A corner index for the face 0 <= \a corner < num_face_corners.
   * \return              The corner number of the \a corner-th vertex of \a face.
   */
  inline int
  element_get_face_corner (const t8_eclass_t tree_class, const t8_element_t *element, const int face,
                           const int corner) const
  {
    return std::visit ([&] (auto &&scheme) { return scheme.element_get_face_corner (element, face, corner); },
                       eclass_schemes[tree_class]);
  };

  /** Return the face numbers of the faces sharing an element's corner.
   * Example quad: 2 x --- x 3
   *                 |     |
   *                 |     |   face 1
   *               0 x --- x 1
   *                  face 2
   * Thus for corner = 1 the output is: face=0 : 2, face=1: 1
   * \param [in] tree_class    The eclass of the current tree.
   * \param [in] element  The element.
   * \param [in] corner   A corner index for the face.
   * \param [in] face     A face index for \a corner.
   * \return              The face number of the \a face-th face at \a corner.
   */
  inline int
  element_get_corner_face (const t8_eclass_t tree_class, const t8_element_t *element, const int corner,
                           const int face) const
  {
    return std::visit ([&] (auto &&scheme) { return scheme.element_get_corner_face (element, corner, face); },
                       eclass_schemes[tree_class]);
  };

  /** Construct the child element of a given number.
   * \param [in] tree_class    The eclass of the current tree.
   * \param [in] elem     This must be a valid element, bigger than maxlevel.
   * \param [in] childid  The number of the child to construct.
   * \param [in,out] child        The storage for this element must exist.
   *                              On output, a valid element.
   * It is valid to call this function with elem = child.
   */
  inline void
  element_get_child (const t8_eclass_t tree_class, const t8_element_t *elem, const int childid,
                     t8_element_t *child) const
  {
    return std::visit ([&] (auto &&scheme) { return scheme.element_get_child (elem, childid, child); },
                       eclass_schemes[tree_class]);
  };

  /** Construct all children of a given element.
   * \param [in] tree_class    The eclass of the current tree.
   * \param [in] elem     This must be a valid element, bigger than maxlevel.
   * \param [in] length   The length of the output array \a c must match
   *                      the number of children.
   * \param [in,out] c    The storage for these \a length elements must exist.
   *                      On output, all children are valid.
   * It is valid to call this function with elem = c[0].
   * \see t8_element_num_children
   */
  inline void
  element_get_children (const t8_eclass_t tree_class, const t8_element_t *elem, const int length,
                        t8_element_t *c[]) const
  {
    return std::visit ([&] (auto &&scheme) { return scheme.element_get_children (elem, length, c); },
                       eclass_schemes[tree_class]);
  };

  /** Compute the child id of an element.
   * \param [in] tree_class    The eclass of the current tree.
   * \param [in] elem     This must be a valid element.
   * \return              The child id of elem.
   */
  inline int
  element_get_child_id (const t8_eclass_t tree_class, const t8_element_t *elem) const
  {
    return std::visit ([&] (auto &&scheme) { return scheme.element_get_child_id (elem); }, eclass_schemes[tree_class]);
  };

  /** Compute the ancestor id of an element, that is the child id
   * at a given level.
   * \param [in] tree_class    The eclass of the current tree.
   * \param [in] elem     This must be a valid element.
   * \param [in] level    A refinement level. Must satisfy \a level <= elem.level
   * \return              The child_id of \a elem in regard to its \a level ancestor.
   * \note The ancestor id at elem.level is the same as the child id.
   */
  inline int
  element_get_ancestor_id (const t8_eclass_t tree_class, const t8_element_t *elem, const int level) const
  {
    return std::visit ([&] (auto &&scheme) { return scheme.element_get_ancestor_id (elem, level); },
                       eclass_schemes[tree_class]);
  };

  /** Query whether element A is an ancestor of the element B.
   * An element A is ancestor of an element B if A == B or if B can 
   * be obtained from A via successive refinement.
   * \param [in] scheme A scheme.
   * \param [in] eclass An eclass.
   * \param [in] element_A An element of class \a eclass in scheme \a scheme.
   * \param [in] element_B An element of class \a eclass in scheme \a scheme.
   * \return     True if and only if \a element_A is an ancestor of \a element_B.
  */
  inline bool
  element_is_ancestor (const t8_eclass_t tree_class, const t8_element_t *element_A, const t8_element_t *element_B) const
  {
    return std::visit ([&] (auto &&scheme) { return scheme.element_is_ancestor (element_A, element_B); },
                       eclass_schemes[tree_class]);
  }

  /** Query whether a given set of elements is a family or not.
   * \param [in] tree_class    The eclass of the current tree.
   * \param [in] fam      An array of as many elements as an element of class
   *                      \a tree_class has siblings.
   * \return              Zero if \a fam is not a family, nonzero if it is.
   * \note level 0 elements do not form a family.
   */
  inline bool
  elements_are_family (const t8_eclass_t tree_class, t8_element_t *const *fam) const
  {
    return std::visit ([&] (auto &&scheme) { return scheme.elements_are_family (fam); }, eclass_schemes[tree_class]);
  };

  /** Compute the nearest common ancestor of two elements. That is,
   * the element with highest level that still has both given elements as
   * descendants.
   * \param [in] tree_class    The eclass of the current tree.
   * \param [in] elem1    The first of the two input elements.
   * \param [in] elem2    The second of the two input elements.
   * \param [in,out] nca  The storage for this element must exist
   *                      and match the element class of the child.
   *                      On output the unique nearest common ancestor of
   *                      \a elem1 and \a elem2.
   */
  inline void
  element_get_nca (const t8_eclass_t tree_class, const t8_element_t *elem1, const t8_element_t *elem2,
                   t8_element_t *const nca) const
  {
    std::visit ([&] (auto &&scheme) { return scheme.element_get_nca (elem1, elem2, nca); }, eclass_schemes[tree_class]);
    T8_ASSERT (element_is_ancestor (tree_class, nca, elem1));
    T8_ASSERT (element_is_ancestor (tree_class, nca, elem2));
  };

  /** Compute the shape of the face of an element.
   * \param [in] tree_class    The eclass of the current tree.
   * \param [in] elem     The element.
   * \param [in] face     A face of \a elem.
   * \return              The element shape of the face.
   * I.e. T8_ECLASS_LINE for quads, T8_ECLASS_TRIANGLE for tets
   *      and depending on the face number either T8_ECLASS_QUAD or
   *      T8_ECLASS_TRIANGLE for prisms.
   */
  inline t8_element_shape_t
  element_get_face_shape (const t8_eclass_t tree_class, const t8_element_t *elem, const int face) const
  {
    return std::visit ([&] (auto &&scheme) { return scheme.element_get_face_shape (elem, face); },
                       eclass_schemes[tree_class]);
  };

  /** Given an element and a face of the element, compute all children of
   * the element that touch the face.
   * \param [in] tree_class    The eclass of the current tree.
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
  element_get_children_at_face (const t8_eclass_t tree_class, const t8_element_t *elem, const int face,
                                t8_element_t *children[], int num_children, int *child_indices) const
  {
    return std::visit (
      [&] (auto &&scheme) {
        return scheme.element_get_children_at_face (elem, face, children, num_children, child_indices);
      },
      eclass_schemes[tree_class]);
  };

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

   * \param [in] tree_class    The eclass of the current tree.
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
  element_face_get_child_face (const t8_eclass_t tree_class, const t8_element_t *elem, const int face,
                               const int face_child) const
  {
    return std::visit ([&] (auto &&scheme) { return scheme.element_face_get_child_face (elem, face, face_child); },
                       eclass_schemes[tree_class]);
  };

  /** Given a face of an element return the face number
   * of the parent of the element that matches the element's face. Or return -1 if
   * no face of the parent matches the face.
   * \param [in] tree_class    The eclass of the current tree.
   * \param [in]  elem    The element.
   * \param [in]  face    Then number of the face.
   * \return              If \a face of \a elem is also a face of \a elem's parent,
   *                      the face number of this face. Otherwise -1.
   * \note For the root element this function always returns \a face.
   */
  inline int
  element_face_get_parent_face (const t8_eclass_t tree_class, const t8_element_t *elem, const int face) const
  {
    return std::visit ([&] (auto &&scheme) { return scheme.element_face_get_parent_face (elem, face); },
                       eclass_schemes[tree_class]);
  };

  /** Given an element and a face of this element. If the face lies on the
   * tree boundary, return the face number of the tree face.
   * If not the return value is arbitrary.
   * You can call \ref t8_element_is_root_boundary to query whether the face is
   * at the tree boundary.
   * \param [in] tree_class    The eclass of the current tree.
   * \param [in] elem     The element.
   * \param [in] face     The index of a face of \a elem.
   * \return The index of the tree face that \a face is a subface of, if
   *         \a face is on a tree boundary.
   *         Any arbitrary integer if \a is not at a tree boundary.
   * \warning The return value may look like a valid face of the tree even if 
   *   the element does not lie on the root boundary.
   */
  inline int
  element_get_tree_face (const t8_eclass_t tree_class, const t8_element_t *elem, const int face) const
  {
    return std::visit ([&] (auto &&scheme) { return scheme.element_get_tree_face (elem, face); },
                       eclass_schemes[tree_class]);
  };

  /** Suppose we have two trees that share a common face f.
   * Given an element e that is a subface of f in one of the trees
   * and given the orientation of the tree connection, construct the face
   * element of the respective tree neighbor that logically coincides with e
   * but lies in the coordinate system of the neighbor tree.
   * \param [in] tree_class    The eclass of the current tree.
   * \param [in] elem1     The face element.
   * \param [in,out] elem2 On return the face element  \a elem1 with respect
   *                       to the coordinate system of the other tree.
   * \param [in] orientation The orientation of the tree-tree connection.
   *                       \see t8_cmesh_set_join
   * \param [in] sign      Depending on the topological orientation of the two tree faces,
   *                       either 0 (both faces have opposite orientation)
   *                       or 1 (both faces have the same top. orientattion).
   *                       \ref t8_eclass_face_orientation
   * \param [in] is_smaller_face Flag to declare whether \a elem1 belongs to
   *                       the smaller face. A face f of tree T is smaller than
   *                       f' of T' if either the eclass of T is smaller or if
   *                       the classes are equal and f<f'. The orientation is
   *                       defined in relation to the smaller face.
   * \note \a elem1 and \a elem2 may point to the same element.
   */
  inline void
  element_transform_face (const t8_eclass_t tree_class, const t8_element_t *elem1, t8_element_t *elem2,
                          const int orientation, const int sign, const int is_smaller_face) const
  {
    return std::visit (
      [&] (auto &&scheme) { return scheme.element_transform_face (elem1, elem2, orientation, sign, is_smaller_face); },
      eclass_schemes[tree_class]);
  };

  /** Given a boundary face inside a root tree's face construct
   *  the element inside the root tree that has the given face as a
   *  face.
   * \param [in] tree_class    The eclass of the current tree.
   * \param [in] face     A face element.
   * \param [in,out] elem An allocated element. The entries will be filled with
   *                      the data of the element that has \a face as a face and
   *                      lies within the root tree.
   * \param [in] root_face The index of the face of the root tree in which \a face
   *                      lies.
   * \return              The face number of the face of \a elem that coincides
   *                      with \a face.
   */
  inline int
  element_extrude_face (const t8_eclass_t tree_class, const t8_element_t *face, t8_element_t *elem,
                        const int root_face) const
  {
    return std::visit ([&] (auto &&scheme) { return scheme.element_extrude_face (face, elem, root_face, this); },
                       eclass_schemes[tree_class]);
  };

  /** Construct the boundary element at a specific face.
   * \param [in] tree_class    The eclass of the current tree.
   * \param [in] elem     The input element.
   * \param [in] face     The index of the face of which to construct the
   *                      boundary element.
   * \param [in,out] boundary An allocated element of dimension of \a element
   *                      minus 1. The entries will be filled with the entries
   *                      of the face of \a element.
   * If \a elem is of class T8_ECLASS_VERTEX, then \a boundary must be NULL
   * and will not be modified.
   */
  inline void
  element_get_boundary_face (const t8_eclass_t tree_class, const t8_element_t *elem, const int face,
                             t8_element_t *boundary) const
  {
    return std::visit ([&] (auto &&scheme) { return scheme.element_get_boundary_face (elem, face, boundary, this); },
                       eclass_schemes[tree_class]);
  };

  /** Construct the first descendant of an element at a given level that touches a given face.
   * \param [in] tree_class    The eclass of the current tree.
   * \param [in] elem      The input element.
   * \param [in] face      A face of \a elem.
   * \param [in, out] first_desc An allocated element. This element's data will be
   *                       filled with the data of the first descendant of \a elem
   *                       that shares a face with \a face.
   * \param [in] level     The level, at which the first descendant is constructed
   */
  inline void
  element_get_first_descendant_face (const t8_eclass_t tree_class, const t8_element_t *elem, const int face,
                                     t8_element_t *first_desc, const int level) const
  {
    return std::visit (
      [&] (auto &&scheme) { return scheme.element_get_first_descendant_face (elem, face, first_desc, level); },
      eclass_schemes[tree_class]);
  };

  /** Construct the last descendant of an element at a given level that touches a given face.
   * \param [in] tree_class    The eclass of the current tree.
   * \param [in] elem      The input element.
   * \param [in] face      A face of \a elem.
   * \param [in, out] last_desc An allocated element. This element's data will be
   *                       filled with the data of the last descendant of \a elem
   *                       that shares a face with \a face.
   * \param [in] level     The level, at which the last descendant is constructed
   */
  inline void
  element_get_last_descendant_face (const t8_eclass_t tree_class, const t8_element_t *elem, const int face,
                                    t8_element_t *last_desc, const int level) const
  {
    return std::visit (
      [&] (auto &&scheme) { return scheme.element_get_last_descendant_face (elem, face, last_desc, level); },
      eclass_schemes[tree_class]);
  };

  /** Compute whether a given element shares a given face with its root tree.
   * \param [in] tree_class    The eclass of the current tree.
   * \param [in] elem     The input element.
   * \param [in] face     A face of \a elem.
   * \return              True if \a face is a subface of the element's root element.
   * \note You can compute the corresponding face number of the tree via \ref element_get_tree_face.
   */
  inline bool
  element_is_root_boundary (const t8_eclass_t tree_class, const t8_element_t *elem, const int face) const
  {
    return std::visit ([&] (auto &&scheme) { return scheme.element_is_root_boundary (elem, face); },
                       eclass_schemes[tree_class]);
  };

  /** Construct the face neighbor of a given element if this face neighbor
   * is inside the root tree. Return 0 otherwise.
   * \param [in] tree_class    The eclass of the current tree.
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
  element_get_face_neighbor_inside (const t8_eclass_t tree_class, const t8_element_t *elem, t8_element_t *neigh,
                                    const int face, int *neigh_face) const
  {
    return std::visit (
      [&] (auto &&scheme) { return scheme.element_get_face_neighbor_inside (elem, neigh, face, neigh_face); },
      eclass_schemes[tree_class]);
  };

  /** Return the shape of an allocated element according its type.
   * For example, a child of an element can be an element of a different shape
   * and has to be handled differently - according to its shape.
   * \param [in] tree_class    The eclass of the current tree.
   * \param [in] elem     The element to be considered
   * \return              The shape of the element as an eclass
   */
  inline t8_element_shape_t
  element_get_shape (const t8_eclass_t tree_class, const t8_element_t *elem) const
  {
    return std::visit ([&] (auto &&scheme) { return scheme.element_get_shape (elem); }, eclass_schemes[tree_class]);
  };

  /** Initialize the entries of an allocated element according to a
   * given linear id in a uniform refinement.
   * \param [in] tree_class    The eclass of the current tree.
   * \param [in,out] elem The element whose entries will be set.
   * \param [in] level    The level of the uniform refinement to consider.
   * \param [in] id       The linear id.
   *                      id must fulfil 0 <= id < 'number of leaves in the uniform refinement'
   */
  inline void
  element_set_linear_id (const t8_eclass_t tree_class, t8_element_t *elem, const int level,
                         const t8_linearidx_t id) const
  {
    return std::visit ([&] (auto &&scheme) { return scheme.element_set_linear_id (elem, level, id); },
                       eclass_schemes[tree_class]);
  };

  /** Compute the linear id of a given element in a hypothetical uniform
   * refinement of a given level.
   * \param [in] tree_class    The eclass of the current tree.
   * \param [in] elem     The element whose id we compute.
   * \param [in] level    The level of the uniform refinement to consider.
   * \return              The linear id of the element.
   */
  inline t8_linearidx_t
  element_get_linear_id (const t8_eclass_t tree_class, const t8_element_t *elem, const int level) const
  {
    return std::visit ([&] (auto &&scheme) { return scheme.element_get_linear_id (elem, level); },
                       eclass_schemes[tree_class]);
  };

  /** Compute the first descendant of a given element.
   * \param [in] tree_class    The eclass of the current tree.
   * \param [in] elem     The element whose descendant is computed.
   * \param [out] desc    The first element in a uniform refinement of \a elem
   *                      of the given level.
   * \param [in] level    The level, at which the descendant is computed.
   */
  inline void
  element_get_first_descendant (const t8_eclass_t tree_class, const t8_element_t *elem, t8_element_t *desc,
                                const int level) const
  {
    return std::visit ([&] (auto &&scheme) { return scheme.element_get_first_descendant (elem, desc, level); },
                       eclass_schemes[tree_class]);
  };

  /** Compute the last descendant of a given element.
   * \param [in] tree_class    The eclass of the current tree.
   * \param [in] elem     The element whose descendant is computed.
   * \param [out] desc    The last element in a uniform refinement of \a elem
   *                      of the given level.
   * \param [in] level    The level, at which the descendant is computed.
   */
  inline void
  element_get_last_descendant (const t8_eclass_t tree_class, const t8_element_t *elem, t8_element_t *desc,
                               const int level) const
  {
    return std::visit ([&] (auto &&scheme) { return scheme.element_get_last_descendant (elem, desc, level); },
                       eclass_schemes[tree_class]);
  };

  /** Construct the successor in a uniform refinement of a given element.
   * \param [in] tree_class    The eclass of the current tree.
   * \param [in] elem1    The element whose successor should be constructed.
   * \param [in,out] elem2  The element whose entries will be set.
   */
  inline void
  element_construct_successor (const t8_eclass_t tree_class, const t8_element_t *t, t8_element_t *s) const
  {
    return std::visit ([&] (auto &&scheme) { return scheme.element_construct_successor (t, s); },
                       eclass_schemes[tree_class]);
  };

  /** Compute the coordinates of a given element vertex inside a reference tree
   * that is embedded into [0,1]^d (d = dimension).
   * \param [in] tree_class    The eclass of the current tree.
   * \param [in] t      The element to be considered.
   * \param [in] vertex The id of the vertex whose coordinates shall be computed.
   * \param [out] coords An array of at least as many doubles as the element's dimension
   *                    whose entries will be filled with the coordinates of \a vertex.
   * \warning           coords should be zero-initialized, as only the first d coords will be set, but when used elsewhere
   *                    all coords might be used. 
   */
  inline void
  element_get_vertex_reference_coords (const t8_eclass_t tree_class, const t8_element_t *t, const int vertex,
                                       double coords[]) const
  {
    return std::visit ([&] (auto &&scheme) { return scheme.element_get_vertex_reference_coords (t, vertex, coords); },
                       eclass_schemes[tree_class]);
  };

  /** Convert points in the reference space of an element to points in the
   *  reference space of the tree.
   * \param [in] tree_class    The eclass of the current tree.
   * \param [in] elem         The element.
   * \param [in] coords_input The coordinates \f$ [0,1]^\mathrm{dim} \f$ of the point
   *                          in the reference space of the element.
   * \param [in] num_coords   Number of \f$ dim\f$-sized coordinates to evaluate.
   * \param [out] out_coords  The coordinates of the points in the
   *                          reference space of the tree.
   */
  inline void
  element_get_reference_coords (const t8_eclass_t tree_class, const t8_element_t *elem, const double *ref_coords,
                                const size_t num_coords, double *out_coords) const
  {
    return std::visit (
      [&] (auto &&scheme) { return scheme.element_get_reference_coords (elem, ref_coords, num_coords, out_coords); },
      eclass_schemes[tree_class]);
  };

  /** Count how many leaf descendants of a given uniform level an element would produce.
   * \param [in] tree_class    The eclass of the current tree.
   * \param [in] t     The element to be checked.
   * \param [in] level A refinement level.
   * \return Suppose \a t is uniformly refined up to level \a level. The return value
   * is the resulting number of elements (of the given level).
   * If \a level < element_get_level(t), the return value should be 0.
   *
   * Example: If \a t is a line element that refines into 2 line elements on each level,
   *  then the return value is max(0, 2^{\a level - level(\a t)}).
   *  Thus, if \a t's level is 0, and \a level = 3, the return value is 2^3 = 8.
   */
  inline t8_gloidx_t
  element_count_leaves (const t8_eclass_t tree_class, const t8_element_t *t, const int level) const
  {
    return std::visit ([&] (auto &&scheme) { return scheme.element_count_leaves (t, level); },
                       eclass_schemes[tree_class]);
  };

  /** Count how many leaf descendants of a given uniform level the root element will produce.
   * \param [in] tree_class    The eclass of the current tree.
   * \param [in] level A refinement level.
   * \return The value of \ref t8_element_count_leaves if the input element
   *      is the root (level 0) element.
   *
   * This is a convenience function, and can be implemented via
   * \ref t8_element_count_leaves.
   */
  inline t8_gloidx_t
  count_leaves_from_root (const t8_eclass_t tree_class, const int level) const
  {
    return std::visit ([&] (auto &&scheme) { return scheme.count_leaves_from_root (level); },
                       eclass_schemes[tree_class]);
  };

#if T8_ENABLE_DEBUG
  /** Query whether a given element can be considered as 'valid' and it is
   *  safe to perform any of the above algorithms on it.
   *  For example this could mean that all coordinates are in valid ranges
   *  and other membervariables do have meaningful values.
   * \param [in] tree_class    The eclass of the current tree.
   * \param [in] elem  The element to be checked.
   * \return          True if \a elem is safe to use. False otherwise.
   * \note            An element that is constructed with \ref t8_element_new
   *                  must pass this test.
   * \note            An element for which \ref t8_element_init was called must pass
   *                  this test.
   * \note            This function is used for debugging to catch certain errors.
   *                  These can for example occur when an element points to a region
   *                  of memory which should not be interpreted as an element.
   * \note            We recommend to use the assertion T8_ASSERT (t8_element_is_valid (elem))
   *                  in the implementation of each of the functions in this file.
   */
  inline int
  element_is_valid (const t8_eclass_t tree_class, const t8_element_t *elem) const
  {
    return std::visit ([&] (auto &&scheme) { return scheme.element_is_valid (elem); }, eclass_schemes[tree_class]);
  };

  /**
 * Print a given element. For a example for a triangle print the coordinates
 * and the level of the triangle. This function is only available in the
 * debugging configuration.
 * \param [in] tree_class    The eclass of the current tree. 
 * \param [in] elem  The element to print
 */
  inline void
  element_debug_print (const t8_eclass_t tree_class, const t8_element_t *elem) const
  {
    return std::visit ([&] (auto &&scheme) { return scheme.element_debug_print (elem); }, eclass_schemes[tree_class]);
  };

  /**
 * Fill a string with readable information about the element
 * \param [in] tree_class    The eclass of the current tree.
 * \param[in] elem The element to translate into human-readable information
 * \param[in, out] debug_string The string to fill. 
 */
  inline void
  element_to_string (const t8_eclass_t tree_class, const t8_element_t *elem, char *debug_string,
                     const int string_size) const
  {
    return std::visit ([&] (auto &&scheme) { return scheme.element_to_string (elem, debug_string, string_size); },
                       eclass_schemes[tree_class]);
  };
#endif

  /** Allocate memory for \a length many elements of a given class and initialize them,
   * and put pointers to the elements in the provided array.
   * \param [in] tree_class    The eclass of the current tree.
   * \param [in] length   The number of elements to be allocated.
   * \param [in,out] elems On input an array of \a length many element pointers.
   *                      On output all these pointers will point to an allocated
   *                      and initialized element.
   * \note There are two ways to create multiple elements of the same type. Create an
   * array of element pointers and fill it with t8_element_new, or allocate memory
   * for \a length times \a element_size many bytes, and fill them with t8_element_init.
   * To access a specific element, offset calculation needs to be done manually, as
   * t8_element_t is incomplete.
   * \note In debugging mode, an element that was created with \ref t8_element_new
   * must pass \ref t8_element_is_valid (for example the root element).
   * \note If an element was created by \ref t8_element_new then \ref t8_element_init
   * may not be called for it. Thus, \ref t8_element_new should initialize an element
   * in the same way as a call to \ref t8_element_init would.
   * \note Every call to \ref t8_element_new must be matched by a call to \ref t8_element_destroy
   * \see t8_element_destroy
   * \see t8_element_init
   * \see t8_element_is_valid
   */
  inline void
  element_new (const t8_eclass_t tree_class, const int length, t8_element_t **elem) const
  {
    return std::visit ([&] (auto &&scheme) { return scheme.element_new (length, elem); }, eclass_schemes[tree_class]);
  };

  /** Initialize an array of allocated elements.
   * \param [in] tree_class    The eclass of the current tree.
   * \param [in] length   The number of elements to be initialized.
   * \param [in,out] elems On input an array of \a length many allocated
   *                       elements.
   * \note In debugging mode, an element that was passed to \ref t8_element_init
   * must pass \ref t8_element_is_valid.
   * \note If an element was created by \ref t8_element_new then \ref t8_element_init
   * may not be called for it. Thus, \ref t8_element_init should initialize an element
   * in the same way as a call to \ref t8_element_new would.
   * \note Every call to \ref t8_element_init must be matched by a call to \ref t8_element_deinit
   * \see t8_element_deinit
   * \see t8_element_new
   * \see t8_element_is_valid
   */
  inline void
  element_init (const t8_eclass_t tree_class, const int length, t8_element_t *elem) const
  {
    return std::visit ([&] (auto &&scheme) { return scheme.element_init (length, elem); }, eclass_schemes[tree_class]);
  };

  /** Deinitialize an array of allocated elements.
   * \param [in] tree_class    The eclass of the current tree.
   * \param [in] length   The number of elements to be deinitialized.
   * \param [in,out] elems On input an array of \a length many allocated
   *                       and initialized elements, on output an array of
   *                       \a length many allocated, but not initialized elements.
   * \note Call this function if you called t8_element_init on the element pointers.
   * \see t8_element_init
   */
  inline void
  element_deinit (const t8_eclass_t tree_class, const int length, t8_element_t *elem) const
  {
    return std::visit ([&] (auto &&scheme) { return scheme.element_deinit (length, elem); },
                       eclass_schemes[tree_class]);
  };

  /** Deallocate an array of elements.
   * \param [in] tree_class    The eclass of the current tree.
   * \param [in] length   The number of elements in the array.
   * \param [in,out] elem On input an array of \a length many allocated
   *                      element pointers.
   *                      On output all these pointers will be freed.
   *                      \a elem itself will not be freed by this function.
   * \see t8_element_new
   */
  inline void
  element_destroy (const t8_eclass_t tree_class, const int length, t8_element_t **elem) const
  {
    return std::visit ([&] (auto &&scheme) { return scheme.element_destroy (length, elem); },
                       eclass_schemes[tree_class]);
  };

  /** create the root element
   * \param [in] tree_class    The eclass of the current tree.
   * \param [in,out] elem The element that is filled with the root
   */
  inline void
  set_to_root (const t8_eclass_t tree_class, t8_element_t *elem) const
  {
    return std::visit ([&] (auto &&scheme) { return scheme.set_to_root (elem); }, eclass_schemes[tree_class]);
  };

  /** Pack multiple elements into contiguous memory, so they can be sent via MPI.
   * \param [in] tree_class    The eclass of the current tree.
   * \param [in] elements Array of elements that are to be packed
   * \param [in] count Number of elements to pack
   * \param [in,out] send_buffer Buffer in which to pack the elements
   * \param [in] buffer_size size of the buffer (in order to check that we don't access out of range)
   * \param [in, out] position the position of the first byte that is not already packed
   * \param [in] comm MPI Communicator
  */
  inline void
  element_MPI_Pack (const t8_eclass_t tree_class, t8_element_t **const elements, const unsigned int count,
                    void *send_buffer, int buffer_size, int *position, sc_MPI_Comm comm) const
  {
    return std::visit (
      [&] (auto &&scheme) {
        return scheme.element_MPI_Pack (elements, count, send_buffer, buffer_size, position, comm);
      },
      eclass_schemes[tree_class]);
  };

  /** Determine an upper bound for the size of the packed message of \a count elements
   * \param [in] tree_class    The eclass of the current tree.
   * \param [in] count Number of elements to pack
   * \param [in] comm MPI Communicator
   * \param [out] pack_size upper bound on the message size
  */
  inline void
  element_MPI_Pack_size (const t8_eclass_t tree_class, const unsigned int count, sc_MPI_Comm comm, int *pack_size) const
  {
    return std::visit ([&] (auto &&scheme) { return scheme.element_MPI_Pack_size (count, comm, pack_size); },
                       eclass_schemes[tree_class]);
  };

  /** Unpack multiple elements from contiguous memory that was received via MPI.
   * \param [in] tree_class    The eclass of the current tree.
   * \param [in] recvbuf Buffer from which to unpack the elements
   * \param [in] buffer_size size of the buffer (in order to check that we don't access out of range)
   * \param [in, out] position the position of the first byte that is not already packed
   * \param [in] elements Array of initialised elements that is to be filled from the message
   * \param [in] count Number of elements to unpack
   * \param [in] comm MPI Communicator
  */
  inline void
  element_MPI_Unpack (t8_eclass_t tree_class, void *recvbuf, const int buffer_size, int *position,
                      t8_element_t **elements, const unsigned int count, sc_MPI_Comm comm) const
  {
    return std::visit (
      [&] (auto &&scheme) { return scheme.element_MPI_Unpack (recvbuf, buffer_size, position, elements, count, comm); },
      eclass_schemes[tree_class]);
  };
};

#endif /* !T8_SCHEME_HXX */
