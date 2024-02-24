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

/** \file t8_element_c_interface.h
 * This file defines the c interface to (some of) the member functions of the 
 * t8_eclass_scheme_c class.
 *
 * We recommend to use the C++ functions directly and only use this
 * interface when you really need to use C.
 */

#ifndef T8_ELEMENT_C_INTERFACE_H
#define T8_ELEMENT_C_INTERFACE_H

#include <t8_element.h>

T8_EXTERN_C_BEGIN ();

/** Return the size of any element of a given class.
 * \return                      The size of an element of class \b ts.
 * We provide a default implementation of this routine that should suffice
 * for most use cases.
 */
size_t
t8_element_size (const t8_eclass_scheme_c *ts);

/** Returns true, if there is one element in the tree, that does not refine into 2^dim children.
 * Returns false otherwise.
 */
int
t8_element_refines_irregular (const t8_eclass_scheme_c *ts);

/** Return the maximum allowed level for any element of a given class.
 * \param [in] ts     Implementation of a class scheme.
 * \return            The maximum allowed level for elements of class \b ts.
 */
int
t8_element_maxlevel (const t8_eclass_scheme_c *ts);

/** Return the type of each child in the ordering of the implementation.
 * \param [in] ts       Implementation of a class scheme.
 * \param [in] childid  Must be between 0 and the number of children (exclusive).
 *                      The number of children is defined in \a t8_element_num_children.
 * \return              The type for the given child.
 */
t8_eclass_t
t8_element_child_eclass (const t8_eclass_scheme_c *ts, int childid);
/** Return the level of a particular element.
 * \param [in] ts      Implementation of a class scheme.
 * \param [in] elem    The element whose level should be returned.
 * \return             The level of \b elem.
 */
int
t8_element_level (const t8_eclass_scheme_c *ts, const t8_element_t *elem);

/** Copy all entries of \b source to \b dest. \b dest must be an existing
 *  element. No memory is allocated by this function.
* \param [in] ts          Implementation of a class scheme.
 * \param [in] source     The element whose entries will be copied to \b dest.
 * \param [in,out] dest   This element's entries will be overwritten with the
 *                        entries of \b source.
 * \note \a source and \a dest may point to the same element.
 */
void
t8_element_copy (const t8_eclass_scheme_c *ts, const t8_element_t *source, t8_element_t *dest);

/** Compare two elements with respect to the scheme.
 * \param [in] ts     Implementation of a class scheme.
 * \param [in] elem1  The first element.
 * \param [in] elem2  The second element.
 * \return       negative if elem1 < elem2, zero if elem1 equals elem2
 *               and positive if elem1 > elem2.
 *  If elem2 is a copy of elem1 then the elements are equal.
 */
int
t8_element_compare (const t8_eclass_scheme_c *ts, const t8_element_t *elem1, const t8_element_t *elem2);

/** Check if two elements are equal.
 * \param [in] ts     Implementation of a class scheme.
 * \param [in] elem1  The first element.
 * \param [in] elem2  The second element.
 * \return            1 if the elements are equal, 0 if they are not equal
 */
int
t8_element_equal (const t8_eclass_scheme_c *ts, const t8_element_t *elem1, const t8_element_t *elem2);

/** Compute the parent of a given element \b elem and store it in \b parent.
 *  \b parent needs to be an existing element. No memory is allocated by this function.
 *  \b elem and \b parent can point to the same element, then the entries of
 *  \b elem are overwritten by the ones of its parent.
 * \param [in] ts     Implementation of a class scheme.
 * \param [in] elem   The element whose parent will be computed.
 * \param [in,out] parent This element's entries will be overwritten by those
 *                    of \b elem's parent.
 *                    The storage for this element must exist
 *                    and match the element class of the parent.
 */
void
t8_element_parent (const t8_eclass_scheme_c *ts, const t8_element_t *elem, t8_element_t *parent);

/** Compute the number of siblings of an element. That is the number of 
 * Children of its parent.
 * \param [in] ts     Implementation of a class scheme.
 * \param [in] elem   The element.
 * \return            The number of siblings of \a element.
 * Note that this number is >= 1, since we count the element itself as a sibling.
 */
int
t8_element_num_siblings (const t8_eclass_scheme_c *ts, const t8_element_t *elem);

/** Compute a specific sibling of a given element \b elem and store it in \b sibling.
 *  \b sibling needs to be an existing element. No memory is allocated by this function.
 *  \b elem and \b sibling can point to the same element, then the entries of
 *  \b elem are overwritten by the ones of its i-th sibling.
 * \param [in] ts     Implementation of a class scheme.
 * \param [in] elem   The element whose sibling will be computed.
 * \param [in] sibid  The id of the sibling computed.
 * \param [in,out] sibling This element's entries will be overwritten by those
 *                    of \b elem's sibid-th sibling.
 *                    The storage for this element must exist
 *                    and match the element class of the sibling.
 */
void
t8_element_sibling (const t8_eclass_scheme_c *ts, const t8_element_t *elem, int sibid, t8_element_t *sibling);

/** Compute the number of corners of an element.
 * \param [in] ts     Implementation of a class scheme.
 * \param [in] elem   The element.
 * \return            The number of corners of \a element.
 */
int
t8_element_num_corners (const t8_eclass_scheme_c *ts, const t8_element_t *elem);

/** Compute the number of faces of an element.
 * \param [in] ts     Implementation of a class scheme.
 * \param [in] elem   The element.
 * \return            The number of faces of \a element.
 */
int
t8_element_num_faces (const t8_eclass_scheme_c *ts, const t8_element_t *elem);

/** Compute the maximum number of faces of a given element and all of its
 *  descendants.
 * \param [in] ts     Implementation of a class scheme.
 * \param [in] elem   The element.
 * \return            The number of faces of \a element.
 */
int
t8_element_max_num_faces (const t8_eclass_scheme_c *ts, const t8_element_t *elem);

/** Compute the number of children of an element when it is refined.
 * \param [in] ts     Implementation of a class scheme.
 * \param [in] elem   The element.
 * \return            The number of children of \a element.
 */
int
t8_element_num_children (const t8_eclass_scheme_c *ts, const t8_element_t *elem);

/** Compute the number of children of an element's face when the element is refined.
 * \param [in] ts     Implementation of a class scheme.
 * \param [in] elem   The element.
 * \param [in] face   A face of \a elem.
 * \return            The number of children of \a face if \a elem is to be refined.
 */
int
t8_element_num_face_children (const t8_eclass_scheme_c *ts, const t8_element_t *elem, int face);

/** Return the corner number of an element's face corner.
 * Example quad: 2 x --- x 3
 *                 |     |
 *                 |     |   face 1
 *               0 x --- x 1
 *      Thus for face = 1 the output is: corner=0 : 1, corner=1: 3
 *
 * \param [in] ts       Implementation of a class scheme.
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
t8_element_get_face_corner (const t8_eclass_scheme_c *ts, const t8_element_t *elem, int face, int corner);

/** Compute the face numbers of the faces sharing an element's corner.
 * Example quad: 2 x --- x 3
 *                 |     |
 *                 |     |   face 1
 *               0 x --- x 1
 *                  face 2
 *      Thus for corner = 1 the output is: face=0 : 2, face=1: 1
 * \param [in] ts       Implementation of a class scheme.
 * \param [in] element  The element.
 * \param [in] corner   A corner index for the face.
 * \param [in] face     A face index for \a corner.
 * \return              The face number of the \a face-th face at \a corner.
 */
int
t8_element_get_corner_face (const t8_eclass_scheme_c *ts, const t8_element_t *elem, int corner, int face);

/** Construct the child element of a given number.
 * \param [in] ts             Implementation of a class scheme.
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
void
t8_element_child (const t8_eclass_scheme_c *ts, const t8_element_t *elem, int childid, t8_element_t *child);

/** Construct all children of a given element.
 * \param [in] ts       Implementation of a class scheme.
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
void
t8_element_children (const t8_eclass_scheme_c *ts, const t8_element_t *elem, int length, t8_element_t *c[]);

/** Compute the child id of an element.
 * \param [in] ts       Implementation of a class scheme.
 * \param [in] elem     This must be a valid element.
 * \return              The child id of elem.
 */
int
t8_element_child_id (const t8_eclass_scheme_c *ts, const t8_element_t *elem);

/** Compute the ancestor id of an element, that is the child id
 * at a given level.
 * \param [in] ts       Implementation of a class scheme.
 * \param [in] elem     This must be a valid element.
 * \param [in] level    A refinement level. Must satisfy \a level < elem.level
 * \return              The child_id of \a elem in regard to its \a level ancestor.
 */
int
t8_element_ancestor_id (const t8_eclass_scheme_c *ts, const t8_element_t *elem, int level);

/** Query whether a given set of elements is a family or not.
 * \param [in] ts       Implementation of a class scheme.
 * \param [in] fam      An array of as many elements as an element of class
 *                      \b ts has children.
 * \return              Zero if \b fam is not a family, nonzero if it is.
 */
int
t8_element_is_family (const t8_eclass_scheme_c *ts, t8_element_t **fam);

/** Compute the nearest common ancestor of two elements. That is,
 * the element with highest level that still has both given elements as
 * descendants.
 * \param [in] ts       Implementation of a class scheme.
 * \param [in] elem1    The first of the two input elements.
 * \param [in] elem2    The second of the two input elements.
 * \param [in,out] nca  The storage for this element must exist
 *                      and match the element class of the child.
 *                      On output the unique nearest common ancestor of
 *                      \b elem1 and \b elem2.
 */
void
t8_element_nca (const t8_eclass_scheme_c *ts, const t8_element_t *elem1, const t8_element_t *elem2, t8_element_t *nca);

/** Compute the shape of the face of an element.
 * \param [in] ts       Implementation of a class scheme.
 * \param [in] elem     The element.
 * \param [in] face     A face of \a elem.
 * \return              The element shape of the face.
 * I.e. T8_ECLASS_LINE for quads, T8_ECLASS_TRIANGLE for tets
 *      and depending on the face number either T8_ECLASS_QUAD or
 *      T8_ECLASS_TRIANGLE for prisms.
 */
t8_element_shape_t
t8_element_face_shape (const t8_eclass_scheme_c *ts, const t8_element_t *elem, int face);

/** Given an element and a face of the element, compute all children of
 * the element that touch the face.
 * \param [in] ts       Implementation of a class scheme.
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
void
t8_element_children_at_face (const t8_eclass_scheme_c *ts, const t8_element_t *elem, int face, t8_element_t *children[],
                             int num_children, int *child_indices);

/** Given a face of an element and a child number of a child of that face, return the face number
 * of the child of the element that matches the child face.
 * \verbatim
 *  x ---- x   x      x           x ---- x
 *  |      |   |      |           |   |  | <-- f
 *  |      |   |      x           |   x--x
 *  |      |   |                  |      |
 *  x ---- x   x                  x ---- x
 *   elem    face  face_child    Returns the face number f
 * \endverbatim
 *
 * \param [in] ts         Implementation of a class scheme.
 * \param [in] elem       The element.
 * \param [in] face       Then number of the face.
 * \param [in] face_child A number 0 <= \a face_child < num_face_children,
 *                        specifying a child of \a elem that shares a face with \a face.
 *                        These children are counted in linear order. This coincides with
 *                        the order of children from a call to \ref t8_element_children_at_face.
 * \return                The face number of the face of a child of \a elem
 *                        that coincides with \a face_child.
 */
int
t8_element_face_child_face (const t8_eclass_scheme_c *ts, const t8_element_t *elem, int face, int face_child);

/** Given a face of an element return the face number
 * of the parent of the element that matches the element's face. Or return -1 if
 * no face of the parent matches the face.
 *
 * \param [in] ts     Implementation of a class scheme.
 * \param [in] elem   The element.
 * \param [in] face   Then number of the face.
 * \return            If \a face of \a elem is also a face of \a elem's parent,
 *                    the face number of this face. Otherwise -1.
 * \note For the root element this function always returns \a face.
 */
int
t8_element_face_parent_face (const t8_eclass_scheme_c *ts, const t8_element_t *elem, int face);

/** Given an element and a face of this element. If the face lies on the
 *  tree boundary, return the face number of the tree face.
 *  If not the return value is arbitrary.
 * \param [in] ts       Implementation of a class scheme.
 * \param [in] elem     The element.
 * \param [in] face     The index of a face of \a elem.
 * \return The index of the tree face that \a face is a subface of, if
 *         \a face is on a tree boundary.
 *         Any arbitrary integer if \a is not at a tree boundary.
 */
int
t8_element_tree_face (const t8_eclass_scheme_c *ts, const t8_element_t *elem, int face);

/** Suppose we have two trees that share a common face f.
 *  Given an element e that is a subface of f in one of the trees
 *  and given the orientation of the tree connection, construct the face
 *  element of the respective tree neighbor that logically coincides with e
 *  but lies in the coordinate system of the neighbor tree.
 * \param [in] ts        Implementation of a class scheme.
 * \param [in] elem1     The face element.
 * \param [in,out] elem2 On return the face element \a elem1 with respect
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
void
t8_element_transform_face (const t8_eclass_scheme_c *ts, const t8_element_t *elem1, t8_element_t *elem2,
                           int orientation, int sign, int is_smaller_face);

/** Given a boundary face inside a root tree's face construct
 *  the element inside the root tree that has the given face as a
 *  face.
 * \param [in] ts       Implementation of a class scheme.
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
int
t8_element_extrude_face (const t8_eclass_scheme_c *ts, const t8_element_t *face, const t8_eclass_scheme_c *face_scheme,
                         t8_element_t *elem, int root_face);

/** Construct the boundary element at a specific face.
 * \param [in] ts       Implementation of a class scheme.
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
void
t8_element_boundary_face (const t8_eclass_scheme_c *ts, const t8_element_t *elem, int face, t8_element_t *boundary,
                          const t8_eclass_scheme_c *boundary_scheme);

/** Construct the first descendant of an element at a given level that touches a given face.
 * \param [in] ts        Implementation of a class scheme.
 * \param [in] elem      The input element.
 * \param [in] face      A face of \a elem.
 * \param [in, out] first_desc An allocated element. This element's data will be
 *                       filled with the data of the first descendant of \a elem
 *                       that shares a face with \a face.
 * \param [in] level     The level, at which the first descendant is constructed
 */
void
t8_element_first_descendant_face (const t8_eclass_scheme_c *ts, const t8_element_t *elem, int face,
                                  t8_element_t *first_desc, int level);

/** Construct the last descendant of an element at a given level that touches a given face.
 * \param [in] ts        Implementation of a class scheme.
 * \param [in] elem      The input element.
 * \param [in] face      A face of \a elem.
 * \param [in, out] last_desc An allocated element. This element's data will be
 *                       filled with the data of the last descendant of \a elem
 *                       that shares a face with \a face.
 * \param [in] level     The level, at which the last descendant is constructed
 */
void
t8_element_last_descendant_face (const t8_eclass_scheme_c *ts, const t8_element_t *elem, int face,
                                 t8_element_t *last_desc, int level);

/** Compute whether a given element shares a given face with its root tree.
 * \param [in] ts       Implementation of a class scheme.
 * \param [in] elem     The input element.
 * \param [in] face     A face of \a elem.
 * \return              True if \a face is a subface of the element's root element.
 */
int
t8_element_is_root_boundary (const t8_eclass_scheme_c *ts, const t8_element_t *elem, int face);

/** Construct the face neighbor of a given element if this face neighbor
 * is inside the root tree. Return 0 otherwise.
 * \param [in] ts     Implementation of a class scheme.
 * \param [in] elem  The element to be considered.
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
int
t8_element_face_neighbor_inside (const t8_eclass_scheme_c *ts, const t8_element_t *elem, t8_element_t *neigh, int face,
                                 int *neigh_face);

/** Return the shape of an allocated element according its type.
*  For example, a child of an element can be an element of a different shape
*  and has to be handled differently - according to its shape.
 * \param [in] ts     Implementation of a class scheme.
*  \param [in] elem   The element to be considered
*  \return            The shape of the element as an eclass
*/
t8_element_shape_t
t8_element_shape (const t8_eclass_scheme_c *ts, const t8_element_t *elem);

/** Initialize the entries of an allocated element according to a
 *  given linear id in a uniform refinement.
 * \param [in] ts       Implementation of a class scheme.
 * \param [in,out] elem The element whose entries will be set.
 * \param [in] level    The level of the uniform refinement to consider.
 * \param [in] id       The linear id.
 *                      id must fulfil 0 <= id < 'number of leaves in the uniform refinement'
 */
void
t8_element_set_linear_id (const t8_eclass_scheme_c *ts, t8_element_t *elem, int level, t8_linearidx_t id);

/** Compute the linear id of a given element in a hypothetical uniform
 * refinement of a given level.
 * \param [in] ts       Implementation of a class scheme.
 * \param [in] elem     The element whose id we compute.
 * \param [in] level    The level of the uniform refinement to consider.
 * \return              The linear id of the element.
 */
t8_linearidx_t
t8_element_get_linear_id (const t8_eclass_scheme_c *ts, const t8_element_t *elem, int level);

/** Compute the first descendant of a given element.
 * \param [in] ts       Implementation of a class scheme.
 * \param [in] elem     The element whose descendant is computed.
 * \param [out] desc    The first element in a uniform refinement of \a elem
 *                      of the maximum possible level.
 */
void
t8_element_first_descendant (const t8_eclass_scheme_c *ts, const t8_element_t *elem, t8_element_t *desc, int level);

/** Compute the last descendant of a given element.
 * \param [in] ts       Implementation of a class scheme.
 * \param [in] elem     The element whose descendant is computed.
 * \param [out] desc    The last element in a uniform refinement of \a elem
 *                      of the maximum possible level.
 */
void
t8_element_last_descendant (const t8_eclass_scheme_c *ts, const t8_element_t *elem, t8_element_t *desc, int level);

/** Construct the successor in a uniform refinement of a given element.
 * \param [in] ts         Implementation of a class scheme.
 * \param [in] elem1      The element whose successor should be constructed.
 * \param [in,out] elem2  The element whose entries will be set.
 * \param [in] level      The level of the uniform refinement to consider.
 */
void
t8_element_successor (const t8_eclass_scheme_c *ts, const t8_element_t *elem1, t8_element_t *elem2, int level);

/** Compute the coordinates of a given element vertex inside a reference tree
 *  that is embedded into [0,1]^d (d = dimension).
 * \param [in] ts       Implementation of a class scheme.
 * \param [in] t        The element to be considered.
 * \param [in] vertex   The id of the vertex whose coordinates shall be computed.
 * \param [out] coords  An array of at least as many doubles as the element's dimension
 *                      whose entries will be filled with the coordinates of \a vertex.
 */
void
t8_element_vertex_reference_coords (const t8_eclass_scheme_c *ts, const t8_element_t *t, const int vertex,
                                    double coords[]);

/* TODO: deactivate */
/** Return a pointer to a t8_element in an array indexed by a size_t.
 * \param [in] array    The \ref sc_array storing \t t8_element_t pointers.
 * \param [in] it       The index of the element that should be returned.
 * \return              A pointer to the it-th element in \b array.
 * We provide a default implementation of this routine that should suffice
 * for most use cases.
 */
/* t8_element_t *t8_element_array_index (sc_array_t *array, size_t it); */

/** Count how many leaf descendants of a given uniform level an element would produce.
 * \param [in] ts     Implementation of a class scheme.
 * \param [in] t      The element to be checked.
 * \param [in] level  A refinement level.
 * \return Suppose \a t is uniformly refined up to level \a level. The return value
 * is the resulting number of elements (of the given level).
 * If \a level < t8_element_level(t), the return value should be 0.
 *
 * Example: If \a t is a line element that refines into 2 line elements on each level,
 *  then the return value is max(0, 2^{\a level - level(\a t)}).
 *  Thus, if \a t's level is 0, and \a level = 3, the return value is 2^3 = 8.
 */
t8_gloidx_t
t8_element_count_leaves (const t8_eclass_scheme_c *ts, const t8_element_t *t, int level);

/** Count how many leaf descendants of a given uniform level the root element will produce.
 * \param [in] ts    Implementation of a class scheme.
 * \param [in] level A refinement level.
 * \return The value of \ref t8_element_count_leaves if the input element
 *      is the root (level 0) element.
 *
 * This is a convenience function, and can be implemented via
 * \ref t8_element_count_leaves.
 */
t8_gloidx_t
t8_element_count_leaves_from_root (const t8_eclass_scheme_c *ts, int level);

/** This function has no defined effect but each implementation is free to
 *  provide its own meaning of it. Thus this function can be used to compute or
 *  lookup very scheme implementation specific data.
 *  \param [in] ts        Implementation of a class scheme.
 *  \param [in] elem      An valid element
 *  \param [in] indata    Pointer to input data
 *  \param [out] outdata  Pointer to output data.
 *  For the correct usage of \a indata and \a outdata see the specific implementations
 *  of the scheme.
 *  For example the default scheme triangle and tetrahedron implementations use 
 *  this function to return the type of a tri/tet to the caller.
 */
void
t8_element_general_function (const t8_eclass_scheme_c *ts, const t8_element_t *elem, const void *indata, void *outdata);

#ifdef T8_ENABLE_DEBUG
/** Query whether a given element can be considered as 'valid' and it is
   *  safe to perform any of the above algorithms on it.
   *  For example this could mean that all coordinates are in valid ranges
   *  and other membervariables do have meaningful values.
   * \param [in] ts   Implementation of a class scheme.
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
int
t8_element_is_valid (const t8_eclass_scheme_c *ts, const t8_element_t *elem);

/**
 * Print a given element. For a example for a triangle print the coordinates
 * and the level of the triangle. This function is only available in the
 * debugging configuration. 
 * 
 * \param [in] ts     Implementation of a class scheme.
 * \param [in] elem   The element to print
 */
void
t8_element_debug_print (const t8_eclass_scheme_c *ts, const t8_element_t *elem);

/**
 * \brief Fill a string with readable information about the element
 * 
 * \param[in] elem The element to translate into human-readable information
 * \param[in, out] debug_string The string to fill. 
 */
void
t8_element_to_string (const t8_eclass_scheme_c *ts, const t8_element_t *elem, char *debug_string,
                      const int string_size);
#endif

/** Allocate memory for an array of elements of a given class and initialize them.
 * \param [in] ts         Implementation of a class scheme.
 * \param [in] length     The number of elements to be allocated.
 * \param [in,out] elems  On input an array of \b length many unallocated element pointers.
 *                        On output all these pointers will point to an allocated and initialized element.
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
void
t8_element_new (const t8_eclass_scheme_c *ts, int length, t8_element_t **elems);

/** Deallocate an array of elements.
 * \param [in] ts         Implementation of a class scheme.
 * \param [in] length     The number of elements in the array.
 * \param [in,out] elems  On input an array of \b length many allocated element pointers. On output all these pointers
 *                        will be freed. \b elem itself will not be freed by this function.
 */
void
t8_element_destroy (const t8_eclass_scheme_c *ts, int length, t8_element_t **elems);

T8_EXTERN_C_END ();

#endif /* !T8_ELEMENT_C_INTERFACE_H */
