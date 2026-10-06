/*
  This file is part of t8code.
  t8code is a C library to manage a collection (a forest) of multiple
  connected adaptive space-trees of general element classes in parallel.

  Copyright (C) 2026 the developers

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

/** \file t8_extruded_scheme.hxx
 * Extruded (quasi-3D) schemes.
 * A 3D tree is treated as the Cartesian product of a 2D base tree and a line in z-direction.
 * All refinement is dictated by the 2D base scheme: a 3D element refines into as many children as its 2D base element,
 * and every element spans the whole height of its tree. The SFC is the SFC of the base scheme.
 *
 * Face numbering: The lateral faces 0, ..., F-1 are the extrusions of the F faces of the base element.
 * Face F is the bottom and face F + 1 the top face. Vertex v is the base vertex v % V at the bottom (v < V) or
 * the top (v >= V), where V is the number of base vertices. This matches the default hex (and prism) conventions.
 *
 * The face elements of lateral faces are quads in the scheme of the face eclass. Since the element spans the whole tree
 * height, the lateral face element stores the in-plane coordinate in x and 0 in y. The z-part is ignored when
 * extruding the face again. This is only consistent, if all trees glued at lateral faces have parallel extrusion
 * directions, see \ref t8_cmesh_is_extrusion_compatible.
 */

#pragma once

#include <t8.h>
#include <sc_containers.h>
#include <t8_element/t8_element.h>
#include <t8_eclass/t8_eclass.h>
#include <t8_schemes/t8_scheme_helpers.hxx>
#include <t8_schemes/t8_default/t8_default_line/t8_dline.h>
#include <t8_schemes/t8_default/t8_default_quad/t8_default_quad.hxx>
#include <cstdio>
#include <cstring>
#include <utility>

/** Scheme for extruded elements.
 * \tparam TEclass     The 3D element class of the extruded trees.
 * \tparam TBaseScheme The scheme of the 2D base elements.
 * \tparam TBaseElem   The element type of \a TBaseScheme.
 */
template <t8_eclass_t TEclass, class TBaseScheme, class TBaseElem>
struct t8_extruded_scheme: public t8_scheme_helpers<TEclass, t8_extruded_scheme<TEclass, TBaseScheme, TBaseElem>>
{
 public:
  /** The eclass of the base elements. */
  static constexpr t8_eclass_t base_eclass = TBaseScheme::get_eclass ();
  /** The number of faces of a base element, which is also the number of lateral faces. */
  static constexpr int num_lateral_faces = t8_eclass_num_faces[base_eclass];
  /** The number of vertices of a base element. */
  static constexpr int num_base_vertices = t8_eclass_num_vertices[base_eclass];
  /** The face number of the bottom face. The top face is \a bottom_face + 1. */
  static constexpr int bottom_face = num_lateral_faces;

  // TODO extruded: do we use static_assert?
  static_assert (t8_eclass_to_dimension[base_eclass] == 2, "The base scheme of an extruded scheme must be 2D.");
  static_assert (t8_eclass_to_dimension[TEclass] == 3, "An extruded scheme must be 3D.");
  static_assert (t8_eclass_num_faces[TEclass] == num_lateral_faces + 2,
                 "The eclass must have the base faces plus a bottom and a top face.");
  static_assert (t8_eclass_num_vertices[TEclass] == 2 * num_base_vertices,
                 "The eclass must have twice as many vertices as the base eclass.");

  // TODO extruded: do we use private? or protected?
 private:
  TBaseScheme base_scheme; /**< The scheme of the 2D base elements. */
  void *scheme_context;    /**< The sc_mempool_t of the extruded elements. */

 public:
  // #################################____CONSTRUCTORS & DESTRUCTOR____#################################################

  /** Constructor. */
  t8_extruded_scheme () noexcept: base_scheme (), scheme_context (sc_mempool_new (sizeof (TBaseElem))) {};

  /** Destructor. */
  ~t8_extruded_scheme ()
  {
    if (scheme_context != nullptr) {
      SC_ASSERT (((sc_mempool_t *) scheme_context)->elem_count == 0);
      sc_mempool_destroy ((sc_mempool_t *) scheme_context);
    }
  }

  /** Move constructor */
  t8_extruded_scheme (t8_extruded_scheme &&other) noexcept
    : base_scheme (std::move (other.base_scheme)), scheme_context (std::exchange (other.scheme_context, nullptr))
  {
  }

  /** Move assignment operator */
  t8_extruded_scheme &
  operator= (t8_extruded_scheme &&other) noexcept
  {
    if (this != &other) {
      // Free existing resources of moved-to object
      if (scheme_context) {
        sc_mempool_destroy ((sc_mempool_t *) scheme_context);
      }
      // TODO extruded: move correct?
      base_scheme = std::move (other.base_scheme);
      scheme_context = std::exchange (other.scheme_context, nullptr);
    }
    return *this;
  }

  /** Copy constructor */
  t8_extruded_scheme (const t8_extruded_scheme &other)
    : base_scheme (other.base_scheme), scheme_context (sc_mempool_new (sizeof (TBaseElem))) {};

  /** Copy assignment operator */
  t8_extruded_scheme &
  operator= (const t8_extruded_scheme &other)
  {
    if (this != &other) {
      // Free existing resources of assigned-to object
      if (scheme_context) {
        sc_mempool_destroy ((sc_mempool_t *) scheme_context);
      }
      base_scheme = other.base_scheme;
      scheme_context = sc_mempool_new (sizeof (TBaseElem));
    }
    return *this;
  }

 private:
  // ################################################____HELPERS____####################################################

  /** Return true if \a face is a lateral face. */
  static constexpr bool
  is_lateral_face (const int face)
  {
    return face < num_lateral_faces;
  }

 public:
  // ################################################____GENERAL INFO____###############################################

  /** Return the size of an extruded element.
   * \return The size of an element.
   */
  static constexpr size_t
  get_element_size (void) noexcept
  {
    return sizeof (TBaseElem);
  }

  /** Returns true, if there is one element in the tree, that does not refine into 2^dim children.
   * \return Dictated by the base scheme because no refinement occurs in extruded direction.
   */
  int
  refines_irregular (void) const
  {
    return base_scheme.refines_irregular ();
  }

  /** Return the maximum allowed level for any element of this class.
   * \return The maximum level of the base scheme.
   */
  inline int
  get_maxlevel (void) const
  {
    return base_scheme.get_maxlevel ();
  }

  /** Return the level of a particular element.
   * \param [in] elem    The element whose level should be returned.
   * \return             The level of \a elem according to the base scheme.
   */
  inline int
  element_get_level (const t8_element_t *elem) const
  {
    T8_ASSERT (element_is_valid (elem));
    return base_scheme.element_get_level (elem);
  }

  /** Return the shape of an element.
   * \param [in] elem The element.
   * \return The shape of the element.
   */
  static constexpr t8_element_shape_t
  element_get_shape ([[maybe_unused]] const t8_element_t *elem)
  {
    return TEclass;
  }

  /** Compute the number of corners of a given element.
   * \param [in] elem The element.
   * \return The number of corners of the element.
   */
  static constexpr int
  element_get_num_corners ([[maybe_unused]] const t8_element_t *elem)
  {
    return t8_eclass_num_vertices[TEclass];
  }

  /** Compute the number of faces of a given element.
   * \param [in] elem The element.
   * \return The number of faces of the element.
   */
  static constexpr int
  element_get_num_faces ([[maybe_unused]] const t8_element_t *elem)
  {
    return t8_eclass_num_faces[TEclass];
  }

  /** Compute the maximum number of faces of a given element and all of its descendants.
   * \param [in] elem The element.
   * \return The maximum number of faces of the element and its descendants.
   */
  static constexpr int
  element_get_max_num_faces ([[maybe_unused]] const t8_element_t *elem)
  {
    return t8_eclass_num_faces[TEclass];
  }

  /** Compute the shape of the face of an element.
   * \param [in] elem     The element.
   * \param [in] face     A face of \a elem.
   * \return              The element shape of the face.
   */
  static constexpr t8_element_shape_t
  element_get_face_shape ([[maybe_unused]] const t8_element_t *elem, const int face)
  {
    return (t8_element_shape_t) t8_eclass_face_types[TEclass][face];
  }

  // ################################################____ORDER____######################################################

  /** Copy all entries of \a source to \a dest.
   * \param [in] source The element whose entries will be copied to \a dest.
   * \param [in,out] dest This element's entries will be overwritten with the entries of \a source.
   * \note \a source and \a dest may point to the same element.
   */
  inline void
  element_copy (const t8_element_t *source, t8_element_t *dest) const
  {
    T8_ASSERT (element_is_valid (source));
    T8_ASSERT (element_is_valid (dest));
    if (source == dest) {
      return;
    }
    base_scheme.element_copy (source, dest);
  }

  /** Compare two elements in the order of the SFC of the base scheme.
   * \param [in] elem1  The first element.
   * \param [in] elem2  The second element.
   * \return            negative if elem1 < elem2, zero if elem1 equals elem2 and positive if elem1 > elem2.
   */
  inline int
  element_compare (const t8_element_t *elem1, const t8_element_t *elem2) const
  {
    T8_ASSERT (element_is_valid (elem1));
    T8_ASSERT (element_is_valid (elem2));
    return base_scheme.element_compare (elem1, elem2);
  }

  /** Check if two elements are equal.
   * \param [in] elem1  The first element.
   * \param [in] elem2  The second element.
   * \return            1 if the elements are equal, 0 if they are not equal
   */
  inline int
  element_is_equal (const t8_element_t *elem1, const t8_element_t *elem2) const
  {
    T8_ASSERT (element_is_valid (elem1));
    T8_ASSERT (element_is_valid (elem2));
    return base_scheme.element_is_equal (elem1, elem2);
  }

  // ################################################____FAMILY____#####################################################

  /** Indicates if an element is refinable.
   * \param [in] elem   The element to check.
   * \return            True if the element in the base scheme is refinable.
   */
  inline bool
  element_is_refinable (const t8_element_t *elem) const
  {
    T8_ASSERT (element_is_valid (elem));
    return base_scheme.element_is_refinable (elem);
  }

  /** Compute the parent of a given element \a elem and store it in \a parent.
   * \param [in] elem   The element whose parent will be computed.
   * \param [in,out] parent This element's entries will be overwritten by those of \a elem's parent.
   * \note \a elem and \a parent can point to the same element.
   */
  inline void
  element_get_parent (const t8_element_t *elem, t8_element_t *parent) const
  {
    T8_ASSERT (element_is_valid (elem));
    base_scheme.element_get_parent (elem, parent);
    T8_ASSERT (element_is_valid (parent));
  }

  /** Compute the number of siblings of an element. That is the number of children of its parent.
   * \param [in] elem The element.
   * \return          The number of siblings of \a elem according to the base scheme.
   */
  inline int
  element_get_num_siblings (const t8_element_t *elem) const
  {
    T8_ASSERT (element_is_valid (elem));
    return base_scheme.element_get_num_siblings (elem);
  }

  /** Compute a specific sibling of a given element \a elem and store it in \a sibling.
   * \param [in] elem   The element whose sibling will be computed.
   * \param [in] sibid  The id of the sibling computed.
   * \param [in,out] sibling This element's entries will be overwritten by those of \a elem's sibid-th sibling.
   * \note \a elem and \a sibling can point to the same element.
   */
  inline void
  element_get_sibling (const t8_element_t *elem, const int sibid, t8_element_t *sibling) const
  {
    T8_ASSERT (element_is_valid (elem));
    base_scheme.element_get_sibling (elem, sibid, sibling);
    T8_ASSERT (element_is_valid (sibling));
  }

  /** Return the number of children of an element when it is refined.
   * \param [in] elem   The element whose number of children is returned.
   * \return            The number of children of \a elem if it is to be refined according to the base scheme.
   */
  inline int
  element_get_num_children (const t8_element_t *elem) const
  {
    T8_ASSERT (element_is_valid (elem));
    return base_scheme.element_get_num_children (elem);
  }

  /** Return the max number of children of an element of this class.
   * \return            The max number of children according to the base scheme.
   */
  inline int
  get_max_num_children () const
  {
    return base_scheme.get_max_num_children ();
  }

  /** Construct the child element of a given number. It is valid to call this function with elem = child.
   * \param [in] elem     This must be a valid element, bigger than maxlevel.
   * \param [in] childid  The number of the child to construct.
   * \param [in,out] child  The storage for this element must exist. On output, a valid element.
   */
  inline void
  element_get_child (const t8_element_t *elem, const int childid, t8_element_t *child) const
  {
    T8_ASSERT (element_is_valid (elem));
    base_scheme.element_get_child (elem, childid, child);
    T8_ASSERT (element_is_valid (child));
  }

  /** Construct all children of a given element. It is valid to call this function with elem = children[0].
   * \param [in] elem          This must be a valid element, bigger than maxlevel.
   * \param [in] length        The length of the output array \a c must match the number of children.
   * \param [in,out] children  The storage for these \a length elements must exist. On output, all children are valid.
   */
  inline void
  element_get_children (const t8_element_t *elem, const int length, t8_element_t *children[]) const
  {
    T8_ASSERT (element_is_valid (elem));
    T8_ASSERT (length == element_get_num_children (elem));
    base_scheme.element_get_children (elem, length, children);
    for (int ichild = 0; ichild < length; ++ichild) {
      T8_ASSERT (element_is_valid (children[ichild]));
    }
  }

  /** Compute the child id of an element.
   * \param [in] elem     This must be a valid element.
   * \return              The child id of elem according to the base scheme.
   */
  inline int
  element_get_child_id (const t8_element_t *elem) const
  {
    T8_ASSERT (element_is_valid (elem));
    return base_scheme.element_get_child_id (elem);
  }

  /** Compute the ancestor id of an element, that is the child id at a given level.
   * \param [in] elem     This must be a valid element.
   * \param [in] level    A refinement level. Must satisfy \a level <= elem.level
   * \return              The child_id of \a elem in regard to its \a level ancestor according to the base scheme.
   */
  inline int
  element_get_ancestor_id (const t8_element_t *elem, const int level) const
  {
    T8_ASSERT (element_is_valid (elem));
    return base_scheme.element_get_ancestor_id (elem, level);
  }

  /** Query whether element A is an ancestor of the element B.
   * \param [in] element_A An element.
   * \param [in] element_B An element.
   * \return     True if and only if \a element_A is an ancestor of \a element_B.
   */
  inline bool
  element_is_ancestor (const t8_element_t *element_A, const t8_element_t *element_B) const
  {
    T8_ASSERT (element_is_valid (element_A));
    T8_ASSERT (element_is_valid (element_B));
    return base_scheme.element_is_ancestor (element_A, element_B);
  }

  /** Query whether a given set of elements is a family or not.
   * \param [in] fam      An array of as many elements as an element of this class has siblings.
   * \return              Zero if \a fam is not a family, nonzero if it is.
   */
  inline int
  elements_are_family (const t8_element_t *const *fam) const
  {
    const int num_siblings = element_get_num_siblings (fam[0]);
    return base_scheme.elements_are_family (fam);
  }

  /** Compute the nearest common ancestor of two elements.
   * \param [in] elem1    The first of the two input elements.
   * \param [in] elem2    The second of the two input elements.
   * \param [in,out] nca  The storage for this element must exist. On output the nearest common ancestor.
   */
  inline void
  element_get_nca (const t8_element_t *elem1, const t8_element_t *elem2, t8_element_t *nca) const
  {
    T8_ASSERT (element_is_valid (elem1));
    T8_ASSERT (element_is_valid (elem2));
    base_scheme.element_get_nca (elem1, elem2, nca);
    T8_ASSERT (element_is_valid (nca));
  }

  // ################################################____FACES____######################################################

  /** Return the number of children of an element's face when the element is refined.
   * \param [in] elem   The element whose face is considered.
   * \param [in] face   A face of \a elem.
   * \return            The number of children of \a face if \a elem is to be refined.
   */
  inline int
  element_get_num_face_children (const t8_element_t *elem, const int face) const
  {
    T8_ASSERT (element_is_valid (elem));
    T8_ASSERT (0 <= face && face < element_get_num_faces (elem));
    if (is_lateral_face (face)) {
      return base_scheme.element_get_num_face_children (elem, face);
    }
    /* All children touch the bottom and the top face. */
    return element_get_num_children (elem);
  }

  /** Return the corner number of an element's face corner.
   * \param [in] elem     The element.
   * \param [in] face     A face index for \a elem.
   * \param [in] corner   A corner index for the face 0 <= \a corner < num_face_corners.
   * \return              The corner number of the \a corner-th vertex of \a face.
   */
  inline int
  element_get_face_corner (const t8_element_t *elem, const int face, const int corner) const
  {
    T8_ASSERT (0 <= face && face < element_get_num_faces (elem));
    if (is_lateral_face (face)) {
      /* A face of the base scheme is always a line, an extruded line is always a quad */
      T8_ASSERT (0 <= corner && corner < 4);
      /* The lateral face is the base face at the bottom (corner 0, 1) and at the top (corner 2, 3). */
      return base_scheme.element_get_face_corner (elem, face, corner % 2) + num_base_vertices * (corner / 2);
    }
    T8_ASSERT (0 <= corner && corner < num_base_vertices);
    return corner + num_base_vertices * (face - bottom_face);
  }

  /** Return the face numbers of the faces sharing an element's corner.
   * \param [in] elem     The element.
   * \param [in] corner   A corner index for the element.
   * \param [in] face     A face index for \a corner, 0 <= \a face < 3.
   * \return              The face number of the \a face-th face at \a corner.
   */
  inline int
  element_get_corner_face (const t8_element_t *elem, const int corner, const int face) const
  {
    T8_ASSERT (0 <= corner && corner < element_get_num_corners (elem));
    T8_ASSERT (0 <= face && face < 3);
    if (face < 2) {
      return base_scheme.element_get_corner_face (elem, corner % num_base_vertices, face);
    }
    return bottom_face + corner / num_base_vertices;
  }

  /** Given an element and a face of the element, compute all children of the element that touch the face.
   * It is valid to call this function with elem = children[0].
   * \param [in] elem     The element.
   * \param [in] face     A face of \a elem.
   * \param [in,out] children Allocated elements, in which the children of \a elem that share a face with \a face
   *                      are stored. They will be stored in order of their linear id.
   * \param [in] num_children The number of elements in \a children. Must match the number of children that touch
   *                      \a face. \ref element_get_num_face_children
   * \param [in,out] child_indices If not NULL, an array of num_children integers must be given,
   *                      on output its i-th entry is the child_id of the i-th face_child.
   */
  inline void
  element_get_children_at_face (const t8_element_t *elem, const int face, t8_element_t *children[],
                                const int num_children, int *child_indices) const
  {
    T8_ASSERT (element_is_valid (elem));
    T8_ASSERT (num_children == element_get_num_face_children (elem, face));
    if (!is_lateral_face (face)) {
      /* All children touch the bottom and the top face. */
      element_get_children (elem, num_children, children);
      if (child_indices != NULL) {
        for (int ichild = 0; ichild < num_children; ++ichild) {
          child_indices[ichild] = ichild;
        }
      }
      return;
    }
    base_scheme.element_get_children_at_face (elem, face, children, num_children, child_indices);
    for (int ichild = 0; ichild < num_children; ++ichild) {
      T8_ASSERT (element_is_valid (children[ichild]));
    }
  }

  /** Given a face of an element and a child number of a child of that face, return the face number
   * of the child of the element that matches the child face.
   * \param [in]  elem    The element.
   * \param [in]  face    The number of the face.
   * \param [in]  face_child A number 0 <= \a face_child < num_face_children.
   * \return              The face number of the face of a child of \a elem that coincides with \a face_child.
   */
  inline int
  element_face_get_child_face (const t8_element_t *elem, const int face, const int face_child) const
  {
    T8_ASSERT (element_is_valid (elem));
    if (is_lateral_face (face)) {
      return base_scheme.element_face_get_child_face (elem, face, face_child);
    }
    return face;
  }

  /** Given a face of an element return the face number of the parent of the element that matches the element's face.
   * Or return -1 if no face of the parent matches the face.
   * \param [in]  elem    The element.
   * \param [in]  face    The number of the face.
   * \return              If \a face of \a elem is also a face of \a elem's parent, the face number of this face.
   *                      Otherwise -1.
   */
  inline int
  element_face_get_parent_face (const t8_element_t *elem, const int face) const
  {
    T8_ASSERT (element_is_valid (elem));
    if (is_lateral_face (face)) {
      return base_scheme.element_face_get_parent_face (elem, face);
    }
    return face;
  }

  /** Given a face of an element and a level coarser than (or equal to) the element's level, return the face number
   * of the ancestor of the element that matches the element's face. Or return -1 if no face of the ancestor matches
   * the face.
   * \param [in]  elem    The element.
   * \param [in]  ancestor_level A refinement level smaller than (or equal to) \a elem's level.
   * \param [in]  face    The number of a face of \a elem.
   * \return              The face number of the ancestor's face, or -1.
   */
  inline int
  element_face_get_ancestor_face (const t8_element_t *elem, const int ancestor_level, const int face) const
  {
    T8_ASSERT (element_is_valid (elem));
    if (is_lateral_face (face)) {
      return base_scheme.element_face_get_ancestor_face (elem, ancestor_level, face);
    }
    return face;
  }

  /** Given an element and a face of this element. If the face lies on the tree boundary, return the face number of
   * the tree face. If not the return value is arbitrary.
   * \param [in] elem     The element.
   * \param [in] face     The index of a face of \a elem.
   * \return The index of the tree face that \a face is a subface of, if \a face is on a tree boundary.
   */
  inline int
  element_get_tree_face (const t8_element_t *elem, const int face) const
  {
    T8_ASSERT (element_is_valid (elem));
    if (is_lateral_face (face)) {
      return base_scheme.element_get_tree_face (elem, face);
    }
    return face;
  }

  /** Transform a face element between two trees.
   * \note Not implemented
   * \param [in] elem1     The face element.
   * \param [in,out] elem2 On return the face element \a elem1 with respect to the coordinate system of the other tree.
   * \param [in] orientation The orientation of the tree-tree connection.
   * \param [in] sign      The topological orientation of the two tree faces.
   * \param [in] is_smaller_face Flag to declare whether \a elem1 belongs to the smaller face.
   */
  inline void
  element_transform_face ([[maybe_unused]] const t8_element_t *elem1, [[maybe_unused]] t8_element_t *elem2,
                          [[maybe_unused]] const int orientation, [[maybe_unused]] const int sign,
                          [[maybe_unused]] const int is_smaller_face) const
  {
    // TODO extruded: "an extruded element is never a face."
    SC_ABORT ("An extruded element is never a face element.\n");
  }

  /** Given a boundary face inside a root tree's face construct the element inside the root tree that has the given
   * face as a face.
   * \note For lateral faces only the in-plane coordinate (x) of the face element is used, since the element spans the whole
   * tree height.
   * \param [in] face     A face element.
   * \param [in,out] elem An allocated element. The entries will be filled with the data of the element that has
   *                      \a face as a face and lies within the root tree.
   * \param [in] root_face The index of the face of the root tree in which \a face lies.
   * \param [in] scheme   The scheme collection with a scheme for the eclass of the face.
   * \return              The face number of the face of \a elem that coincides with \a face.
   */
  inline int
  element_extrude_face (const t8_element_t *face, t8_element_t *elem, const int root_face,
                        const t8_scheme *scheme) const
  {
    T8_ASSERT (0 <= root_face && root_face < element_get_num_faces (elem));
    if (!is_lateral_face (root_face)) {
      /* The bottom and top faces are elements of the base scheme. */
      base_scheme.element_copy (face, elem);
      T8_ASSERT (element_is_valid (elem));
      return root_face;
    }
    // TODO extruded:
    /* Build the boundary line of the base element from the in-plane coordinate of the face and extrude it. */
    const p4est_quadrant_t *face_quad = (const p4est_quadrant_t *) face;
    t8_dline_t base_face;
    base_face.level = face_quad->level;
    base_face.x = ((int64_t) face_quad->x * T8_DLINE_ROOT_LEN) / P4EST_ROOT_LEN;
    const int base_face_number
      = base_scheme.element_extrude_face ((const t8_element_t *) &base_face, elem, root_face, scheme);
    T8_ASSERT (element_is_valid (elem));
    return base_face_number;
  }

  /** Construct the first descendant of an element at a given level that touches a given face.
   * \param [in] elem      The input element.
   * \param [in] face      A face of \a elem.
   * \param [in, out] first_desc An allocated element. This element's data will be filled with the data of the first
   *                       descendant of \a elem that shares a face with \a face.
   * \param [in] level     The level, at which the first descendant is constructed
   */
  inline void
  element_get_first_descendant_face (const t8_element_t *elem, const int face, t8_element_t *first_desc,
                                     const int level) const
  {
    T8_ASSERT (element_is_valid (elem));
    if (is_lateral_face (face)) {
      base_scheme.element_get_first_descendant_face (elem, face, first_desc, level);
      T8_ASSERT (element_is_valid (first_desc));
    }
    else {
      element_get_first_descendant (elem, first_desc, level);
    }
  }

  /** Construct the last descendant of an element at a given level that touches a given face.
   * \param [in] elem      The input element.
   * \param [in] face      A face of \a elem.
   * \param [in, out] last_desc An allocated element. This element's data will be filled with the data of the last
   *                       descendant of \a elem that shares a face with \a face.
   * \param [in] level     The level, at which the last descendant is constructed
   */
  inline void
  element_get_last_descendant_face (const t8_element_t *elem, const int face, t8_element_t *last_desc,
                                    const int level) const
  {
    T8_ASSERT (element_is_valid (elem));
    if (is_lateral_face (face)) {
      base_scheme.element_get_last_descendant_face (elem, face, last_desc, level);
      T8_ASSERT (element_is_valid (last_desc));
    }
    else {
      element_get_last_descendant (elem, last_desc, level);
    }
  }

  /** Construct the boundary element at a specific face.
   * \note The bottom and top face elements are copies of the base element. The lateral face element is a quad with the
   * in-plane coordinate in x and 0 in y.
   * \param [in] elem     The input element.
   * \param [in] face     The index of the face of which to construct the boundary element.
   * \param [in,out] boundary An allocated element of dimension of \a elem minus 1. The entries will be filled with the
   *                      entries of the face of \a elem.
   * \param [in] scheme   The scheme containing an eclass scheme for the boundary face.
   */
  inline void
  element_get_boundary_face (const t8_element_t *elem, const int face, t8_element_t *boundary,
                             const t8_scheme *scheme) const
  {
    T8_ASSERT (element_is_valid (elem));
    T8_ASSERT (0 <= face && face < element_get_num_faces (elem));
    if (!is_lateral_face (face)) {
      base_scheme.element_copy (elem, boundary);
      return;
    }
    // TODO extruded: ??
    /* Compute the boundary line of the base element and use its coordinate as the in-plane coordinate. */
    t8_dline_t base_face;
    base_face.x = 0;
    base_face.level = 0;
    base_scheme.element_get_boundary_face (elem, face, (t8_element_t *) &base_face, scheme);
    p4est_quadrant_t *face_quad = (p4est_quadrant_t *) boundary;
    face_quad->level = base_face.level;
    face_quad->x = ((int64_t) base_face.x * P4EST_ROOT_LEN) / T8_DLINE_ROOT_LEN;
    face_quad->y = 0;
  }

  /** Compute whether a given element shares a given face with its root tree.
   * \param [in] elem     The input element.
   * \param [in] face     A face of \a elem.
   * \return              True if \a face is a subface of the element's root element.
   */
  inline int
  element_is_root_boundary (const t8_element_t *elem, const int face) const
  {
    T8_ASSERT (element_is_valid (elem));
    if (is_lateral_face (face)) {
      return base_scheme.element_is_root_boundary (elem, face);
    }
    /* Every element spans the whole tree height. */
    return true;
  }

  /** Construct the face neighbor of a given element if this face neighbor is inside the root tree.
   * \param [in] elem The element to be considered.
   * \param [in,out] neigh If the face neighbor of \a elem along \a face is inside the root tree, this element's data
   *                  is filled with the data of the face neighbor. Otherwise the data is the neighbor outside of the
   *                  root tree.
   * \param [in] face The number of the face along which the neighbor should be constructed.
   * \param [out] neigh_face The number of \a face as viewed from \a neigh.
   * \return          True if \a neigh is inside the root tree.
   */
  inline int
  element_get_face_neighbor_inside (const t8_element_t *elem, t8_element_t *neigh, const int face,
                                    int *neigh_face) const
  {
    T8_ASSERT (element_is_valid (elem));
    T8_ASSERT (element_is_valid (neigh));
    T8_ASSERT (neigh_face != NULL);
    if (is_lateral_face (face)) {
      const int is_inside = base_scheme.element_get_face_neighbor_inside (elem, neigh, face, neigh_face);
      T8_ASSERT (element_is_valid (neigh));
      return is_inside;
    }
    // TODO extruded: ??
    base_scheme.element_copy (elem, neigh);
    *neigh_face = face == bottom_face ? bottom_face + 1 : bottom_face;
    /* Since elements span the whole tree height, the neighbor is always outside. */
    return false;
  }

  // ################################################____SFC____########################################################

  /** Initialize the entries of an allocated element according to a linear id in a uniform refinement
      given by the base scheme.
   * \param [in,out] elem The element whose entries will be set.
   * \param [in] level    The level of the uniform refinement to consider.
   * \param [in] id       The linear id. 0 <= id < 'number of leaves in the uniform refinement'
   */
  inline void
  element_set_linear_id (t8_element_t *elem, const int level, const t8_linearidx_t id) const
  {
    base_scheme.element_set_linear_id (elem, level, id);
    T8_ASSERT (element_is_valid (elem));
  }

  /** Compute the linear id of a given element in a hypothetical uniform refinement of a given level.
   * \param [in] elem     The element whose id we compute.
   * \param [in] level    The level of the uniform refinement to consider.
   * \return              The linear id of the element.
   */
  inline t8_linearidx_t
  element_get_linear_id (const t8_element_t *elem, const int level) const
  {
    T8_ASSERT (element_is_valid (elem));
    return base_scheme.element_get_linear_id (elem, level);
  }

  /** Compute the first descendant of a given element.
   * \param [in] elem     The element whose descendant is computed.
   * \param [out] desc    The first element in a uniform refinement of \a elem of the given level.
   * \param [in] level    The level, at which the descendant is computed.
   */
  inline void
  element_get_first_descendant (const t8_element_t *elem, t8_element_t *desc, const int level) const
  {
    T8_ASSERT (element_is_valid (elem));
    base_scheme.element_get_first_descendant (elem, desc, level);
    T8_ASSERT (element_is_valid (desc));
  }

  /** Compute the last descendant of a given element.
   * \param [in] elem     The element whose descendant is computed.
   * \param [out] desc    The last element in a uniform refinement of \a elem of the given level.
   * \param [in] level    The level, at which the descendant is computed.
   */
  inline void
  element_get_last_descendant (const t8_element_t *elem, t8_element_t *desc, const int level) const
  {
    T8_ASSERT (element_is_valid (elem));
    base_scheme.element_get_last_descendant (elem, desc, level);
    T8_ASSERT (element_is_valid (desc));
  }

  /** Construct the successor in a uniform refinement of a given element.
   * \param [in] elem      The element whose successor should be constructed.
   * \param [in,out] succ  The successor element whose entries will be set.
   */
  inline void
  element_construct_successor (const t8_element_t *elem, t8_element_t *succ) const
  {
    T8_ASSERT (element_is_valid (elem));
    base_scheme.element_construct_successor (elem, succ);
    T8_ASSERT (element_is_valid (succ));
  }

  /** Count how many leaf descendants of a given uniform level an element would produce.
   * \param [in] elem      The element to be checked.
   * \param [in] level     A refinement level.
   * \return The number of elements of level \a level in a uniform refinement of \a elem.
   */
  inline t8_gloidx_t
  element_count_leaves (const t8_element_t *elem, const int level) const
  {
    return base_scheme.element_count_leaves (elem, level);
  }

  /** Count how many leaf descendants of a given uniform level the root element will produce.
   * \param [in] level A refinement level.
   * \return The number of elements of level \a level in a uniform refinement of the root.
   */
  inline t8_gloidx_t
  count_leaves_from_root (const int level) const
  {
    return base_scheme.count_leaves_from_root (level);
  }

  // ################################################____COORDINATES____################################################

  /** Compute the coordinates of a given element vertex inside a reference tree that is embedded into [0,1]^3.
   * \param [in] elem    The element to be considered.
   * \param [in] vertex  The id of the vertex whose coordinates shall be computed.
   * \param [out] coords An array of at least 3 doubles, filled with the coordinates of \a vertex.
   */
  inline void
  element_get_vertex_reference_coords (const t8_element_t *elem, const int vertex, double coords[]) const
  {
    T8_ASSERT (element_is_valid (elem));
    T8_ASSERT (0 <= vertex && vertex < element_get_num_corners (elem));
    base_scheme.element_get_vertex_reference_coords (elem, vertex % num_base_vertices, coords);
    coords[2] = vertex / num_base_vertices;
  }

  /** Convert points in the reference space of an element to points in the reference space of the tree.
   * \param [in] elem         The element.
   * \param [in] ref_coords   The coordinates in [0,1]^3 of the points in the reference space of the element.
   * \param [in] num_coords   Number of 3-sized coordinates to evaluate.
   * \param [out] out_coords  The coordinates of the points in the reference space of the tree.
   */
  inline void
  element_get_reference_coords (const t8_element_t *elem, const double *ref_coords, const size_t num_coords,
                                double *out_coords) const
  {
    T8_ASSERT (element_is_valid (elem));
    for (size_t icoord = 0; icoord < num_coords; ++icoord) {
      const double *ref = ref_coords + 3 * icoord;
      double *out = out_coords + 3 * icoord;
      /* The base scheme reads 3D input and writes 2D output coordinates. */
      double base_out[2];
      base_scheme.element_get_reference_coords (elem, ref, 1, base_out);
      out[0] = base_out[0];
      out[1] = base_out[1];
      out[2] = ref[2];
    }
  }

  // ################################################____DEBUG____######################################################

#if T8_ENABLE_DEBUG
  /** Query whether a given element can be considered as 'valid'.
   * \param [in]      elem  The element to be checked.
   * \return          True if \a elem is safe to use. False otherwise.
   */
  inline int
  element_is_valid (const t8_element_t *elem) const
  {
    return base_scheme.element_is_valid (elem);
  }

  /** Print a given element.
   * \param [in]        elem  The element to print
   */
  inline void
  element_debug_print (const t8_element_t *elem) const
  {
    char debug_string[BUFSIZ];
    element_to_string (elem, debug_string, BUFSIZ);
    t8_debugf ("%s\n", debug_string);
  }
#else
  /** Query whether a given element can be considered as 'valid'. Only checked in debugging mode.
   * \param [in]      elem  The element to be checked.
   * \return          Always true.
   */
  static constexpr int
  element_is_valid ([[maybe_unused]] const t8_element_t *elem)
  {
    return true;
  }
#endif

  /** Fill a string with the coordinates and levels of an element.
   * \param [in]        elem  The element to print
   * \param [in]        debug_string  String to fill
   * \param [in]        string_size  String size of \a debug_string.
   */
  inline void
  element_to_string (const t8_element_t *elem, char *debug_string, const int string_size) const
  {
    T8_ASSERT (debug_string != NULL);
    base_scheme.element_to_string (elem, debug_string, string_size);
    const size_t base_length = strlen (debug_string);
    if (base_length < (size_t) string_size) {
      snprintf (debug_string + base_length, string_size - base_length, ", z-level: 0");
    }
  }

  // ################################################____MEMORY____#####################################################

  /** Allocate memory for an array of elements and initialize them to the root element.
   * \param [in] length   The number of elements to be allocated.
   * \param [in,out] elem On input an array of \a length many unallocated element pointers.
   *                      On output all these pointers will point to an allocated and initialized element.
   */
  inline void
  element_new (const int length, t8_element_t **elem) const
  {
    T8_ASSERT (0 <= length);
    for (int ielem = 0; ielem < length; ++ielem) {
      elem[ielem] = (t8_element_t *) sc_mempool_alloc ((sc_mempool_t *) scheme_context);
      base_scheme.element_init (1, elem[ielem]);
      set_to_root (elem[ielem]);
    }
  }

  /** Initialize an array of allocated elements.
   * \param [in] length   The number of elements to be initialized.
   * \param [in,out] elem On input an array of \a length many allocated elements (contiguous memory).
   */
  inline void
  element_init (const int length, t8_element_t *elem) const
  {
    TBaseElem *elements = (TBaseElem *) elem;
    for (int ielem = 0; ielem < length; ++ielem) {
      t8_element_t *element = (t8_element_t *) (elements + ielem);
      base_scheme.element_init (1, element);
    }
  }

  /** Deinitialize an array of allocated elements.
   * \param [in] length   The number of elements to be deinitialized.
   * \param [in,out] elem On input an array of \a length many allocated and initialized elements (contiguous memory).
   */
  inline void
  element_deinit (const int length, t8_element_t *elem) const
  {
    TBaseElem *elements = (TBaseElem *) elem;
    for (int ielem = 0; ielem < length; ++ielem) {
      base_scheme.element_deinit (1, (t8_element_t *) (elements + ielem));
    }
  }

  /** Deallocate an array of elements.
   * \param [in] length   The number of elements in the array.
   * \param [in,out] elem On input an array of \a length many allocated element pointers.
   */
  inline void
  element_destroy (const int length, t8_element_t **elem) const
  {
    T8_ASSERT (0 <= length);
    for (int ielem = 0; ielem < length; ++ielem) {
      sc_mempool_free ((sc_mempool_t *) scheme_context, elem[ielem]);
    }
  }

  /** Fill an element with the root element.
   * \param [in,out] elem   The element to be filled with root.
   */
  inline void
  set_to_root (t8_element_t *elem) const
  {
    base_scheme.set_to_root (elem);
  }

  // ################################################____MPI____########################################################

  /** Pack multiple elements into contiguous memory, so they can be sent via MPI.
   * \note First all base elements are packed, then all lines.
   * \param [in] elements        Array of elements that are to be packed
   * \param [in] count           Number of elements to pack
   * \param [in,out] send_buffer Buffer in which to pack the elements
   * \param [in] buffer_size     size of the buffer (in order to check that we don't access out of range)
   * \param [in, out] position   the position of the first byte that is not already packed
   * \param [in] comm            MPI Communicator
   */
  inline void
  element_MPI_Pack (t8_element_t **const elements, const unsigned int count, void *send_buffer, const int buffer_size,
                    int *position, sc_MPI_Comm comm) const
  {
    base_scheme.element_MPI_Pack (elements, count, send_buffer, buffer_size, position, comm);
  }

  /** Determine an upper bound for the size of the packed message of \a count elements
   * \param [in] count      Number of elements to pack
   * \param [in] comm       MPI Communicator
   * \param [out] pack_size upper bound on the message size
   */
  inline void
  element_MPI_Pack_size (const unsigned int count, sc_MPI_Comm comm, int *pack_size) const
  {
    base_scheme.element_MPI_Pack_size (count, comm, pack_size);
    int datasize = 0;
    int line_size = 0;
    int mpiret = sc_MPI_Pack_size (1, sc_MPI_INT, comm, &datasize);
    SC_CHECK_MPI (mpiret);
    line_size += datasize;
    mpiret = sc_MPI_Pack_size (1, sc_MPI_INT8_T, comm, &datasize);
    SC_CHECK_MPI (mpiret);
    line_size += datasize;
    *pack_size += count * line_size;
  }

  /** Unpack multiple elements from contiguous memory that was received via MPI.
   * \param [in] recvbuf        Buffer from which to unpack the elements
   * \param [in] buffer_size    size of the buffer (in order to check that we don't access out of range)
   * \param [in, out] position  the position of the first byte that is not already packed
   * \param [in] elements       Array of initialised elements that is to be filled from the message
   * \param [in] count          Number of elements to unpack
   * \param [in] comm           MPI Communicator
   */
  inline void
  element_MPI_Unpack (void *recvbuf, const int buffer_size, int *position, t8_element_t **elements,
                      const unsigned int count, sc_MPI_Comm comm) const
  {
    base_scheme.element_MPI_Unpack (recvbuf, buffer_size, position, elements, count, comm);
  }
};

/** The extruded hex scheme: default quads extruded in z-direction. */
using t8_extruded_scheme_hex = t8_extruded_scheme<T8_ECLASS_HEX, t8_default_scheme_quad, p4est_quadrant_t>;
