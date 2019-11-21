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

/** \file t8_element_cxx.hxx
 * This file defines basic operations on an element in a refinement tree.
 *
 * All operations work for all element classes by providing a virtual function table.
 * For each element class, one implementation of the type and virtual table is required.
 */

#ifndef T8_ELEMENT_CXX_HXX
#define T8_ELEMENT_CXX_HXX

#include <sc_refcount.h>
#include <t8_eclass.h>
#include <t8_element.h>

T8_EXTERN_C_BEGIN ();

/* TODO: Implement a set of rules that have to hold between different eclass,
 *       i.e. lines must have a greater or equal maxlevel than quads and triangles.
 *       Check whether this rules are fulfilled in the construction of a scheme.
 */
/* TODO: Implement a test that boundary and extrude leads to the original element. */

/** This struct holds virtual functions for a particular element class. */
struct t8_eclass_scheme
{
  /** This scheme defines the operations for a particular element class. */
protected:
  size_t element_size;                          /**< The size in bytes of an element of class \a eclass */
  void               *ts_context;               /**< Anonymous implementation context. */

public:
                      t8_eclass_t eclass;
                              /**< The element class */

  /** The destructor. It does nothing but has to be defined since
   * we may want to delete an eclass_scheme that is actually inherited
   * (for example t8_default_scheme_quad) and providing and implementation
   * for the destructor ensures that the
   * destructor of the child class will be executed. */
                      virtual ~ t8_eclass_scheme ()
  {
  }

  /** The virtual table for a particular implementation of an element class. */

  /** Return the size of any element of a given class.
   * \return                      The size of an element of class \b ts.
   * We provide a default implementation of this routine that should suffice
   * for most use cases.
   */
  virtual size_t      t8_element_size (void);

  /** Return the maximum allowed level for any element of a given class.
   * \return                      The maximum allowed level for elements of class \b ts.
   */
  virtual int         t8_element_maxlevel (void) = 0;

  /** Return the type of each child in the ordering of the implementation.
   * \param [in] childid  Must be between 0 and the number of children (exclusive).
   *                      The number of children is defined in \a t8_element_num_children.
   * \return              The type for the given child.
   */
  virtual t8_eclass_t t8_element_child_eclass (int childid) = 0;

  /** Return the level of a particular element.
   * \param [in] elem    The element whose level should be returned.
   * \return             The level of \b elem.
   */
  virtual int         t8_element_level (const t8_element_t * elem) = 0;

  /** Copy all entries of \b source to \b dest. \b dest must be an existing
   *  element. No memory is allocated by this function.
   * \param [in] source The element whose entries will be copied to \b dest.
   * \param [in,out] dest This element's entries will be overwritted with the
   *                    entries of \b source.
   * \note \a source and \a dest may point to the same element.
   */
  virtual void        t8_element_copy (const t8_element_t * source,
                                       t8_element_t * dest) = 0;

  /** Compare two elements.
   * \param [in] elem1  The first element.
   * \param [in] elem2  The second element.
   * \return       negativ if elem1 < elem2, zero if elem1 equals elem2
   *               and positiv if elem1 > elem2.
   *  If elem2 is a copy of elem1 then the elements are equal.
   */
  virtual int         t8_element_compare (const t8_element_t * elem1,
                                          const t8_element_t * elem2) = 0;

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
  virtual void        t8_element_parent (const t8_element_t * elem,
                                         t8_element_t * parent) = 0;

  /** Compute a specific sibling of a given element \b elem and store it in \b sibling.
   *  \b sibling needs to be an existing element. No memory is allocated by this function.
   *  \b elem and \b sibling can point to the same element, then the entries of
   *  \b elem are overwritten by the ones of its i-th sibling.
   * \param [in] elem   The element whose parent will be computed.
   * \param [in] sibid  The id of the sibling computed.
   * \param [in,out] sibling This element's entries will be overwritten by those
   *                    of \b elem's sibid-th sibling.
   *                    The storage for this element must exist
   *                    and match the element class of the sibling.
   *                    For a pyramid, for example, it may be either a
   *                    tetrahedron or a pyramid depending on \b sibid.
   */
  virtual void        t8_element_sibling (const t8_element_t * elem,
                                          int sibid,
                                          t8_element_t * sibling) = 0;

  /** Compute the number of corners of a given element.
   * \param [in] elem The element.
   * \return          The number of corners of \a elem.
   */
  virtual int         t8_element_num_corners (const t8_element_t * elem) = 0;

  /** Compute the number of faces of a given element.
   * \param [in] elem The element.
   * \return          The number of faces of \a elem.
   */
  virtual int         t8_element_num_faces (const t8_element_t * elem) = 0;

  /** Compute the maximum number of faces of a given element and all of its
   *  descendants.
   * \param [in] elem The element.
   * \return          The maximum number of faces of \a elem and its descendants.
   */
  virtual int         t8_element_max_num_faces (const t8_element_t * elem) =
    0;

  /** Return the number of children of an element when it is refined.
   * \param [in] elem   The element whose number of children is returned.
   * \return            The number of children of \a elem if it is to be refined.
   */
  virtual int         t8_element_num_children (const t8_element_t * elem) = 0;

  /** Return the number of children of an element's face when the element is refined.
   * \param [in] elem   The element whose face is considered.
   * \param [in] face   A face of \a elem.
   * \return            The number of children of \a face if \a elem is to be refined.
   */
  virtual int         t8_element_num_face_children (const t8_element_t *
                                                    elem, int face) = 0;

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
  /* TODO: Prism order, Pyramid order. */
  virtual int         t8_element_get_face_corner (const t8_element_t *
                                                  element, int face,
                                                  int corner) = 0;

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
  virtual int         t8_element_get_corner_face (const t8_element_t *
                                                  element, int corner,
                                                  int face) = 0;

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
  virtual void        t8_element_child (const t8_element_t * elem,
                                        int childid, t8_element_t * child) =
    0;

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
  virtual void        t8_element_children (const t8_element_t * elem,
                                           int length, t8_element_t * c[]) =
    0;

  /** Compute the child id of an element.
   * \param [in] elem     This must be a valid element.
   * \return              The child id of elem.
   */
  virtual int         t8_element_child_id (const t8_element_t * elem) = 0;

  /** Compute the ancestor id of an element, that is the child id
   * at a given level.
   * \param [in] elem     This must be a valid element.
   * \param [in] level    A refinement level. Must satisfy \a level < elem.level
   * \return              The child_id of \a elem in regard to its \a level ancestor.
   */
  virtual int         t8_element_ancestor_id (const t8_element_t * elem,
                                              int level) = 0;

  /** Query whether a given set of elements is a family or not.
   * \param [in] fam      An array of as many elements as an element of class
   *                      \b ts has children.
   * \return              Zero if \b fam is not a family, nonzero if it is.
   */
  virtual int         t8_element_is_family (t8_element_t ** fam) = 0;

  /* TODO: This could be problematic for pyramids, since elem1 and elem2
   *       could be of different classes. Would need two eclass_schemes as input */
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
  virtual void        t8_element_nca (const t8_element_t * elem1,
                                      const t8_element_t * elem2,
                                      t8_element_t * nca) = 0;

  /** Compute the elmement class of the face of an element.
   * \param [in] elem     The element.
   * \param [in] face     A face of \a elem.
   * \return              The element class of the face.
   * I.e. T8_ECLASS_LINE for quads, T8_ECLASS_TRIANGLE for tets
   *      and depending on the face number either T8_ECLASS_QUAD or
   *      T8_ECLASS_TRIANGLE for prisms.
   */
  virtual t8_eclass_t t8_element_face_class (const t8_element_t * elem,
                                             int face) = 0;

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
  virtual void        t8_element_children_at_face (const t8_element_t * elem,
                                                   int face,
                                                   t8_element_t * children[],
                                                   int num_children,
                                                   int *child_indices) = 0;

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
   *                      that conincides with \a face_child.
   */
  virtual int         t8_element_face_child_face (const t8_element_t * elem,
                                                  int face, int face_child) =
    0;

    /** Given a face of an element return the face number
     * of the parent of the element that matches the element's face. Or return -1 if
     * no face of the parent matches the face.

     * \param [in]  elem    The element.
     * \param [in]  face    Then number of the face.
     * \return              If \a face of \a elem is also a face of \a elem's parent,
     *                      the face number of this face. Otherwise -1.
     * \note For the root element this function always returns \a face.
     */
  virtual int         t8_element_face_parent_face (const t8_element_t * elem,
                                                   int face) = 0;

  /** Given an element and a face of this element. If the face lies on the
   *  tree boundary, return the face number of the tree face.
   *  If not the return value is arbitrary.
   * \param [in] elem     The element.
   * \param [in] face     The index of a face of \a elem.
   * \return The index of the tree face that \a face is a subface of, if
   *         \a face is on a tree boundary.
   *         Any arbitrary integer if \a is not at a tree boundary.
   */
  virtual int         t8_element_tree_face (const t8_element_t * elem,
                                            int face) = 0;

  /** Suppose we have two trees that share a common face f.
   *  Given an element e that is a subface of f in one of the trees
   *  and given the orientation of the tree connection, construct the face
   *  element of the respective tree neighbor that logically coincides with e
   *  but lies in the coordinate system of the neighbor tree.
   *  \param [in] elem1     The face element.
   *  \param [in,out] elem2 On return the face elment \a elem1 with respective
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
  virtual void        t8_element_transform_face (const t8_element_t * elem1,
                                                 t8_element_t * elem2,
                                                 int orientation,
                                                 int sign,
                                                 int is_smaller_face) = 0;

  /** Given a boundary face inside a root tree's face construct
   *  the element inside the root tree that has the given face as a
   *  face.
   * \param [in] face     A face element.
   * \param [in,out] elem An allocated element. The entries will be filled with
   *                      the data of the element that has \a face as a face and
   *                      lies within the root tree.
   * \param [in] root_face The index of the face of the root tree in which \a face
   *                      lies.
   * \return              The face number of the face of \a elem that coincides
   *                      with \a face.
   */
  /* TODO: update documentation with face_scheme */
  virtual int         t8_element_extrude_face (const t8_element_t * face,
                                               const t8_eclass_scheme_c *
                                               face_scheme,
                                               t8_element_t * elem,
                                               int root_face) = 0;

  /** Construct the boundary element at a specific face.
   * \param [in] elem     The input element.
   * \param [in] face     The index of the face of which to construct the
   *                      boundary element.
   * \param [in,out] boundary An allocated element of dimension of \a element
   *                      minus 1. The entries will be filled with the entries
   *                      of the face of \a element.
   * If \a elem is of class T8_ECLASS_VERTEX, then \a boundary must be NULL
   * and will not be modified.
   */
  /* TODO: update documentation with boundary_scheme */
  virtual void        t8_element_boundary_face (const t8_element_t * elem,
                                                int face,
                                                t8_element_t * boundary,
                                                const t8_eclass_scheme_c *
                                                boundary_scheme) = 0;

  /** Construct the first descendant of an element that touches a given face.
   * \param [in] elem      The input element.
   * \param [in] face      A face of \a elem.
   * \param [in, out] first_desc An allocated element. This element's data will be
   *                       filled with the data of the first descendant of \a elem
   *                       that shares a face with \a face.
   */
  /* TODO: Add a level and call with forest->maxlevel */
  virtual void        t8_element_first_descendant_face (const t8_element_t *
                                                        elem, int face,
                                                        t8_element_t *
                                                        first_desc,
                                                        int level) = 0;

  /** Construct the last descendant of an element that touches a given face.
   * \param [in] elem      The input element.
   * \param [in] face      A face of \a elem.
   * \param [in, out] last_desc An allocated element. This element's data will be
   *                       filled with the data of the last descendant of \a elem
   *                       that shares a face with \a face.
   */
  /* TODO: Add a level and call with forest->maxlevel */
  virtual void        t8_element_last_descendant_face (const t8_element_t *
                                                       elem, int face,
                                                       t8_element_t *
                                                       last_desc, int level) =
    0;

  /* TODO: document better */
  /* TODO: document better.
   *        Do we need this functino at all?
   *        If not remove it. If so, what to do with prisms and pyramids?
   *        Here the boundary elements are of different eclasses, so we cannot
   *        store them in an array...
   */
/** Construct all codimension-one boundary elements of a given element. */
  virtual void        t8_element_boundary (const t8_element_t * elem,
                                           int min_dim, int length,
                                           t8_element_t ** boundary) = 0;

  /** Compute whether a given element shares a given face with its root tree.
   * \param [in] elem     The input element.
   * \param [in] face     A face of \a elem.
   * \return              True if \a face is a subface of the element's root element.
   */
  virtual int         t8_element_is_root_boundary (const t8_element_t * elem,
                                                   int face) = 0;

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
  virtual int         t8_element_face_neighbor_inside (const t8_element_t *
                                                       elem,
                                                       t8_element_t * neigh,
                                                       int face,
                                                       int *neigh_face) = 0;

  /** Initialize the entries of an allocated element according to a
   *  given linear id in a uniform refinement.
   * \param [in,out] elem The element whose entries will be set.
   * \param [in] level    The level of the uniform refinement to consider.
   * \param [in] id       The linear id.
   *                      id must fulfil 0 <= id < 'number of leafs in the uniform refinement'
   */
  virtual void        t8_element_set_linear_id (t8_element_t * elem,
                                                int level,
                                                t8_linearidx_t id) = 0;

  /** Compute the linear id of a given element in a hypothetical uniform
   * refinement of a given level.
   * \param [in] elem     The element whose id we compute.
   * \param [in] level    The level of the uniform refinement to consider.
   * \return              The linear id of the element.
   */
  virtual t8_linearidx_t t8_element_get_linear_id (const
                                                   t8_element_t *
                                                   elem, int level) = 0;

  /** Compute the first descendant of a given element.
   * \param [in] elem     The element whose descendant is computed.
   * \param [out] desc    The first element in a uniform refinement of \a elem
   *                      of the maximum possible level.
   */
  /* TODO: Add a level and call with forest->maxlevel */
  virtual void        t8_element_first_descendant (const t8_element_t *
                                                   elem,
                                                   t8_element_t * desc,
                                                   int level) = 0;

  /** Compute the last descendant of a given element.
   * \param [in] elem     The element whose descendant is computed.
   * \param [out] desc    The last element in a uniform refinement of \a elem
   *                      of the maximum possible level.
   */
  virtual void        t8_element_last_descendant (const t8_element_t *
                                                  elem,
                                                  t8_element_t * desc,
                                                  int level) = 0;

  /** Construct the successor in a uniform refinement of a given element.
   * \param [in] elem1    The element whose successor should be constructed.
   * \param [in,out] elem2  The element whose entries will be set.
   * \param [in] level    The level of the uniform refinement to consider.
   */
  virtual void        t8_element_successor (const t8_element_t * t,
                                            t8_element_t * s, int level) = 0;

/** Get the integer coordinates of the anchor node of an element */
  /* TODO: better document this */
  virtual void        t8_element_anchor (const t8_element_t * elem,
                                         int anchor[3]) = 0;

  /** Compute the root length of a given element, that is the length of
   * its level 0 ancestor.
   * \param [in] elem     The element whose root length should be computed.
   * \return              The root length of \a elem
   */
  virtual int         t8_element_root_len (const t8_element_t * elem) = 0;

  /** Compute the integer coordinates of a given element vertex.
   *   \param [in] ts     The virtual table for this element class.
   *   \param [in] t      The element to be considered.
   *   \param [in] vertex The id of the vertex whose coordinates shall be computed.
   *   \param [out] coords An array of at least as many integers as the element's dimension
   *                      whose entries will be filled with the coordinates of \a vertex.
   */
  virtual void        t8_element_vertex_coords (const t8_element_t * t,
                                                int vertex, int coords[]) = 0;

  /* TODO: deactivate */
  /** Return a pointer to a t8_element in an array indexed by a size_t.
   * \param [in] array    The \ref sc_array storing \t t8_element_t pointers.
   * \param [in] it       The index of the element that should be returned.
   * \return              A pointer to the it-th element in \b array.
   * We provide a default implementation of this routine that should suffice
   * for most use cases.
   */
  virtual t8_element_t *t8_element_array_index (sc_array_t * array,
                                                size_t it);

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
   * \note            We recommend to use the assertion T8_ASSERT (t8_element_is_valid (elem))
   *                  in the implementation of each of the functions in this file.
   */
  virtual int         t8_element_is_valid (const t8_element_t * elem) const =
    0;
#endif

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
  virtual void        t8_element_new (int length, t8_element_t ** elem) = 0;

 /** Initialize an array of allocated elements.
   * \param [in] length   The number of elements to be allocated.
   * \param [in,out] elems On input an array of \b length many allocated
   *                       elements.
   * \param [in] called_new True if the elements in \a elem were created by a call
   *                       to \ref t8_element_new. False if no element in \a elem
   *                       was created in this way. The case that only some elements
   *                       were created by \ref t8_element_new should never occur.
   * \note In debugging mode, an element that was passed to \ref t8_element_init
   * must pass \ref t8_element_is_valid.
   * \note If an element was created by \ref t8_element_new then \ref t8_element_init
   * may not be called for it. Thus, \ref t8_element_new should initialize an element
   * in the same way as a call to \ref t8_element_init would.
   * Thus, if \a called_new is true this function should usually do nothing.
   * \see t8_element_new
   * \see t8_element_is_valid
   */
  virtual void        t8_element_init (int length, t8_element_t * elem,
                                       int called_new) = 0;

  /** Deallocate an array of elements.
   * \param [in] ts       The virtual table for this element class.
   * \param [in] length   The number of elements in the array.
   * \param [in,out] elems On input an array of \b length many allocated
   *                      element pointers.
   *                      On output all these pointers will be freed.
   *                      \b elem itself will not be freed by this function.
   */
  virtual void        t8_element_destroy (int length,
                                          t8_element_t ** elem) = 0;
};

/* TODO: document */
void                t8_scheme_cxx_destroy (t8_scheme_cxx_t * s);

/* TODO: Copy the doxygen comments to the class definition above,
 * then delete all the functions below */
#if 0
/** Destroy an implementation of a particular element class. */
void                t8_eclass_scheme_destroy (t8_eclass_scheme_t * ts);

/* TODO: This function does not exist yet in the eclass_scheme class */
/** Allocate a set of elements suitable for the boundary of a given class.
 * \param [in] scheme           Defines the implementation of the element class.
 * \param [in] theclass         The element class whose boundary we want.
 * \param [in] min_dim          Ignore boundary points of lesser dimension.
 * \param [in] length           Must be equal to the return value
 *                              of \ref t8_eclass_count_boundary.
 * \param [in,out] boundary     On input, array of element pointers of at
 *                              least length \b length.  Filled on output.
 */
void                t8_eclass_boundary_new (t8_scheme_t * scheme,
                                            t8_eclass_t theclass, int min_dim,
                                            int length,
                                            t8_element_t ** boundary);

/* TODO: This function does not exist yet in the eclass_scheme class */
/** Destroy a set of elements suitable for the boundary of a given class.
 * \param [in] scheme           Defines the implementation of the element class.
 * \param [in] theclass         The element class whose boundary we have.
 * \param [in] min_dim          Ignore boundary points of lesser dimension.
 * \param [in] length           Must be equal to the return value
 *                              of \ref t8_eclass_count_boundary.
 * \param [in,out] boundary     Array of element pointers holding elements
 *                              as created by \ref t8_eclass_boundary_new.
 *                              The elements are destroyed by this function.
 */
void                t8_eclass_boundary_destroy (t8_scheme_t * scheme,
                                                t8_eclass_t theclass,
                                                int min_dim, int length,
                                                t8_element_t ** boundary);
#endif /* if 0 */

T8_EXTERN_C_END ();

#endif /* !T8_ELEMENT_CXX_HXX */
