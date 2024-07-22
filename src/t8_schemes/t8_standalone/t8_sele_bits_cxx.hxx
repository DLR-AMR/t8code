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

#ifndef T8_SELE_BITS_CXX_HXX
#define T8_SELE_BITS_CXX_HXX

#include "t8_element.h"
#include "t8_sele_cxx.hxx"
#include <sc_functions.h>

/** Initialize a pyramid as the pyramid with a given global id in a uniform
 *  refinement of a given level. *
 * \param [in,out] p  Existing pyramid whose data will be filled.
 * \param [in] level  level of uniform grid to be considered.
 * \param [in] id     Index to be considered.
 */
template <t8_eclass_t eclass_T>
void
t8_sele_init_linear_id (t8_standalone_element_t<eclass_T> *p, const int level, t8_linearidx_t id);

/** Copy the data from source to dest
 * \param[in] source    The source-pyramid
 * \param[in,out] dest  The destination
 */
template <t8_eclass_t eclass_T>
void
t8_sele_copy (const t8_standalone_element_t<eclass_T> *source, t8_standalone_element_t<eclass_T> *dest);

/** Computes the linear position of a pyramid in an uniform grid.
 * \param [in] p          pyramid whose id will be computed.
 * \param [in] level      The level on which the linear-id should be computed.
 * \return                Returns the linear position of this pyramid on a grid.
 */
template <t8_eclass_t eclass_T>
t8_linearidx_t
t8_sele_linear_id (const t8_standalone_element_t<eclass_T> *p, const int level);

/** Compute the childid-th child in Morton order of a pyramid.
 * \param [in] elem         Input pyramid.
 * \param [in,out] childid  The id of the child, 0..7 in Morton order.
 * \param [out] child       Existing pyramid whose data will be filled
 * 		                      with the date of t's childid-th child.
 */
template <t8_eclass_t eclass_T>
void
t8_sele_child (const t8_standalone_element_t<eclass_T> *elem, const int child_id,
               t8_standalone_element_t<eclass_T> *child);

/** Compute the children of a pyramid, array version
 * \param [in] p        Input pyramid
 * \param [in, out] c   Pointers to the computed children in Morton order
 */
template <t8_eclass_t eclass_T>
void
t8_sele_children (const t8_standalone_element_t<eclass_T> *p, t8_standalone_element_t<eclass_T> **c);

/** Given a pyramid and a face, compute all children touching this face
 * \param [in] p                  Input pyramid
 * \param [in] face               The face to compute the childran at
 * \param [in, out] children      The children of \a p at \a face
 * \param [in] num_childrem       The number of children at this face
 * \param [in, out] child_indices An array to be filled with the local-ids of the children.
 */
template <t8_eclass_t eclass_T>
void
t8_sele_children_at_face (const t8_standalone_element_t<eclass_T> *p, const int face,
                          t8_standalone_element_t<eclass_T> *children[], const int num_children, int *child_indices);

/** Given a face of a pyramid and a child number of a child of that face,
 * return the face number of the child of the pyramid that matches the child
 * face.
 * \param[in] p               Input pyramid
 * \param[in] face            A face of \a p
 * \param [in] face_child     A number specifying a child on the \a face
 * \returns                   The number of the face of the child \a face_child
 */
template <t8_eclass_t eclass_T>
int
t8_sele_face_child_face (const t8_standalone_element_t<eclass_T> *p, const int face, const int face_child);

/** Given the facenumber of a pyramid, return the shape of the face
 * \param[in] pyra       Input pyramid
 * \param[in] face       The facenumber
 * \return               the shape of the face  
 */
template <t8_eclass_t eclass_T>
t8_element_shape_t
t8_sele_face_shape (const t8_standalone_element_t<eclass_T> *pyra, int face);

/** Returns the corner number of a pyramid given a face of a pyramid and
 *  a corner number regarding that face.
 * \param [in]  pyra    Input pyramid
 * \param [in]  face    The facenumber of a face of \a pyra
 * \param [in]  corner  The cornernumber of a corner of \a face
 * \return              The cornernumber of \a pyra
 */
template <t8_eclass_t eclass_T>
int
t8_sele_get_face_corner (const t8_standalone_element_t<eclass_T> *pyra, int face, int corner);

/** Compare two elements. returns negative if p1 < p2, zero if p1 equals p2
 *  and positive if p1 > p2.
 *  If p2 is a copy of p1 then the elements are equal.
 * \param[in] p1    A pyramid
 * \param[in] p2    Another pyramid
 * \returns         an integer describing which pyramid is larger.
 */
template <t8_eclass_t eclass_T>
int
t8_sele_compare (const t8_standalone_element_t<eclass_T> *p1, const t8_standalone_element_t<eclass_T> *p2);

template <t8_eclass_t eclass_T>
int
t8_sele_equal (const t8_standalone_element_t<eclass_T> *p1, const t8_standalone_element_t<eclass_T> *p2);

/** Check whether a collection of 10 pyramids is a family in Morton order.
 * \param [in]  fam A collection of pyramids
 * \return      Nonzero if \a fam is a family of pyramids
 */
template <t8_eclass_t eclass_T>
int
t8_sele_is_family (t8_standalone_element_t<eclass_T> **fam);

/** Compute whether a given pyramid shares a given face with its root tree.
 * \param [in]  p       The input pyramid
 * \param [in]  face    A face of \a p
 * \return              True, if \a is a subface of the pyramid root element.
 */
template <t8_eclass_t eclass_T>
int
t8_sele_is_root_boundary (const t8_standalone_element_t<eclass_T> *p, const int face);

/** Compute the neighbor of p along a given face and the number of the dual face if
 * the neighbor is inside the root pyramid.
 * \param[in] p                 Input pyramid
 * \param[in, out] neigh        The neighbor of \a p
 * \param[in] face              The face of \a p along which \a neigh is computed
 * \param [in, out] neigh_face  The dual face
 */
template <t8_eclass_t eclass_T>
int
t8_sele_face_neighbor (const t8_standalone_element_t<eclass_T> *p, t8_standalone_element_t<eclass_T> *neigh,
                       const int face, int *neigh_face);

template <t8_eclass_t eclass_T>
void
t8_sele_children_at_face (const t8_standalone_element_t<eclass_T> elem, int face,
                          t8_standalone_element_t<eclass_T> **children, int num_children, int *child_indices);

/** Compute the position of the ancestor of this child at level \a level within
 * its siblings.
 * \param [in] p  pyramid to be considered.
 * \return Returns its child id in 0..9
 */
template <t8_eclass_t eclass_T>
int
t8_sele_child_id (const t8_standalone_element_t<eclass_T> *p);

/** Returns zero if p is not inside root, 1 ow
 * \param[in] p       Pyramid to check
 * \returns 0 if p is inside root, 1, ow
 */
template <t8_eclass_t eclass_T>
int
t8_sele_is_inside_root (const t8_standalone_element_t<eclass_T> *p);

/** Check, if a tet of type 0 or 3 has a common face with its pyramid-ancestor
 * \param [in] p      input pyramid
 * \param [in] face   A face of \a p.
 * \return            false if they don't share a face, true otherwise
 */
template <t8_eclass_t eclass_T>
int
t8_sele_tet_boundary (const t8_standalone_element_t<eclass_T> *p, const int face);

/** compute if a given element lies on the tree boundary and return the face number
 * of the tree face. If not the return value is arbitrary
 * \param [in] elem     pyramid
 * \param [in] face     a face of \a elem
 * \return              See description
 */
template <t8_eclass_t eclass_T>
int
t8_sele_tree_face (const t8_standalone_element_t<eclass_T> *p, const int face);

/** Compute the first descendant of a pyramid at a given level. This is the descendant of
 * the pyramid in a uniform level refinement that has the smallest id.
 * \param [in] p        pyramid whose descendant is computed.
 * \param [out] desc    Existing pyramid whose data will be filled with the data
 *                      of \a p's first descendant on level \a level.
 * \param [in] level    The refinement level. Must be greater than \a p's refinement
 *                      level.
 */
template <t8_eclass_t eclass_T>
void
t8_sele_first_descendant (const t8_standalone_element_t<eclass_T> *p, t8_standalone_element_t<eclass_T> *desc,
                          const int level);

/** Construct the first descendant of a pyramid touching a given face
 * \param [in] p        pyramid whose descendant is computed.
 * \param [in] face     The face at which the descendant is computed
 * \param [out] first_desc       Existing pyramid whose data will be filled with the data
 *                      of \a p's first descendant on level \a level.
 * \param [in] level    The refinement level. Must be greater than \a p's refinement
 *                      level
 */
template <t8_eclass_t eclass_T>
void
t8_sele_first_descendant_face (const t8_standalone_element_t<eclass_T> *p, const int face,
                               t8_standalone_element_t<eclass_T> *first_desc, const int level);

/** Compute the last descendant of a pyramid at a given level. This is the descendant of
 * the pyramid in a uniform level refinement that has the largest id.
 * \param [in] p        pyramid whose descendant is computed.
 * \param [out] desc    Existing pyramid whose data will be filled with the data
 *                      of \a p's last descendant on level \a level.
 * \param [in] level    The refinement level. Must be greater than \a p's refinement
 *                      level.
 */
template <t8_eclass_t eclass_T>
void
t8_sele_last_descendant (const t8_standalone_element_t<eclass_T> *p, t8_standalone_element_t<eclass_T> *desc,
                         int level);

/** Construct the last descendant of a pyramid touching a given face
 * \param [in] p            pyramid whose descendant is computed.
 * \param [in] face         The face at which the descendant is computed
 * \param [out] last_desc   Existing pyramid whose data will be filled with the data
 *                          of \a p's first descendant on level \a level.
 * \param [in] level        The refinement level. Must be greater than \a p's refinement
 *                          level.
 */
template <t8_eclass_t eclass_T>
void
t8_sele_last_descendant_face (const t8_standalone_element_t<eclass_T> *p, const int face,
                              t8_standalone_element_t<eclass_T> *last_desc, const int level);

/** Compute the coordinates of a vertex of a element.
 * \param [in] p    Input pyramid.
 * \param [in] vertex The number of the vertex.
 * \param [out] coords An array of 3 t8_element_coord_t that
 * 		     will be filled with the coordinates of the vertex.
 */
template <t8_eclass_t eclass_T>
void
t8_sele_compute_coords (const t8_standalone_element_t<eclass_T> *p, const int vertex, int coords[]);

/** Compute the parent of a given pyramid
 * \param [in] p        Input pyramid.
 * \param [out] parent  The parent of \a p.
 */
template <t8_eclass_t eclass_T>
void
t8_sele_parent (const t8_standalone_element_t<eclass_T> *p, t8_standalone_element_t<eclass_T> *parent);

/**
 * Compute the number of corners of a pyramid. If pyramid has type less than 6,
 * it is actually a tetrahedron.
 * \param [in] p    Input pyramid.
 * \return          The number of corners of p.
 */
template <t8_eclass_t eclass_T>
int
t8_sele_num_corners (const t8_standalone_element_t<eclass_T> *p);

/** Compute the number of children of p
 * \param [in] p    Input pyramid.
 * \return          The number of children of p.
 */
template <t8_eclass_t eclass_T>
int
t8_sele_num_children (const t8_standalone_element_t<eclass_T> *p);

/** Compute the number of siblings of p
 * \param [in] p    Input pyramid
 * \return          The number of siblings of p.
 */
template <t8_eclass_t eclass_T>
int
t8_sele_num_siblings (const t8_standalone_element_t<eclass_T> *p);

/** Return the number of faces of p
 * \param [in] p    Input pyramid
 * \return          The number of faces of p
 */
template <t8_eclass_t eclass_T>
int
t8_sele_num_faces (const t8_standalone_element_t<eclass_T> *p);

/** Return the maximal number of faces of an element p
 * \param [in] p    Input pyramid
 * \return          The maximal number of faces of p
 */
template <t8_eclass_t eclass_T>
int
t8_sele_max_num_faces (const t8_standalone_element_t<eclass_T> *p);

/** Given a face of an element return the face number
 * of the parent of the element that matches the element's face. Or return -1 if
 * no face of the parent matches the face.
 * \param [in] elem Input pyramid
 * \param [in] face a face of \a elem
 * \return          the facenumber of the parent of \a elem matching \a face or -1
 */
template <t8_eclass_t eclass_T>
int
t8_sele_face_parent_face (const t8_standalone_element_t<eclass_T> *elem, const int face);

template <t8_eclass_t eclass_T>
void
t8_sele_transform_face (const t8_standalone_element_t<eclass_T> *elem1, t8_standalone_element_t<eclass_T> *elem2,
                        int orientation, int sign, int smaller_faces);

template <t8_eclass_t eclass_T, t8_eclass_t face_eclass_T>
void
t8_sele_boundary_face (const t8_standalone_element_t<eclass_T> *elem, const int face,
                       t8_standalone_element_t<face_eclass_T> *boundary);

template <t8_eclass_t eclass_T, t8_eclass_t face_eclass_T>
int
t8_sele_extrude_face (const t8_standalone_element_t<face_eclass_T> *face, t8_standalone_element_t<eclass_T> *elem,
                      int root_face);

/** Return the child-id of the ancestor of p at level level
 * \param [in] p    Input pyramid
 * \param [in] level The ancestor-level
 * \return          The child-id of the ancestor
 */
template <t8_eclass_t eclass_T>
int
t8_sele_ancestor_id (const t8_standalone_element_t<eclass_T> *p, const int level);

/**
 * Compute the ancestor of \a el at a given level via a loop over the parent function.
 * Used in debug mode to compare against the equation based calculation.
 * 
 * \param[in] el        Input element
 * \param[in] level     Level of the ancestor to compute
 * \param[in, out] anc  Allocated element that will be filled with the data of the ancestor.
 */
template <t8_eclass_t eclass_T>
void
t8_sele_ancestor_loop (const t8_standalone_element_t<eclass_T> *el, const int level,
                       t8_standalone_element_t<eclass_T> *anc);
/**
 * Compute the ancestor of \a el at a given level via the equation properties
 * 
 * \param[in] pyra      Input pyramid
 * \param[in] level     Level of the ancestor to compute
 * \param[in, out] anc  Allocated element that will be filled with the data of the ancestor.
 */
template <t8_eclass_t eclass_T>
void
t8_sele_ancestor_equation (const t8_standalone_element_t<eclass_T> *el, const int level,
                           t8_standalone_element_t<eclass_T> *anc);

/** Compute the type of a pyramid at a given level. Starting from its own level,
 * we iterate over the levels and compute the type of this level. If p is a tetrahedron,
 * we compute it in a tetrahedral fashion up unto the last level where p is a tet and
 * continue in a pyramidal fashion 
 * \param [in] p      Input pyramid
 * \param [in] level  The level at which the type is computed
 * \return            The type of \a p at level \a level.
 */
template <t8_eclass_t eclass_T>
t8_element_type_t<eclass_T>
t8_sele_compute_type_at_level (const t8_standalone_element_t<eclass_T> *p, const int level);

/** Returns the shape of the pyramid (pyramid or tetrahedron)
 * \param [in] p    Input pyramid.
 * \return          The eclass of the element
 */
template <t8_eclass_t eclass_T>
t8_element_shape_t
t8_sele_shape (const t8_standalone_element_t<eclass_T> *p);

/** Computes the successor of a pyramid in a uniform grid of level \a level.
 * \param [in] elem  pyramid whose id will be computed.
 * \param [in,out] s Existing pyramid whose data will be filled with the
 *                data of \a l's successor on level \a level.
 */
template <t8_eclass_t eclass_T>
void
t8_sele_successor (const t8_standalone_element_t<eclass_T> *elem, t8_standalone_element_t<eclass_T> *s,
                   const int level);

/** Compute the reference coordinates of a vertex of a pyramid when the
 * tree (level 0 triangle) is embedded in [0,1]^3.
 * \param [in] elem    Input pyramid.
 * \param [in] vertex The number of the vertex.
 * \param [out] coordinates An array of 3 double that
 * 		     will be filled with the reference coordinates of the vertex.
 */
template <t8_eclass_t eclass_T>
void
t8_sele_vertex_reference_coords (const t8_standalone_element_t<eclass_T> *elem, const int vertex, double coords[]);

/** Compute the nearest common ancestor of two elements
 * \param [in]      el1       The first pyramid
 * \param [in]      el2       The second pyramid
 * \param [in,out]  nca         Existing pyramid whose data will be filled with
 *                              the data of \a el1 and \a el2 nearest common ancestor.
 */
template <t8_eclass_t eclass_T>
void
t8_sele_nearest_common_ancestor (const t8_standalone_element_t<eclass_T> *el1,
                                 const t8_standalone_element_t<eclass_T> *el2, t8_standalone_element_t<eclass_T> *nca);

/** Query whether all entries of a pyramid are in valid ranges.
 * A pyramid is valid if and only if its triangle and line member are valid.
 * \param [in] p  pyramid to be considered.
 * \return        True, if \a p is a valid pyramid and it is safe to call any
 *                function in this file on \a p.
 *                False otherwise.
 */
template <t8_eclass_t eclass_T>
int
t8_sele_is_valid (const t8_standalone_element_t<eclass_T> *p);

/** Print the coordinates, the level and the type of a pyramid
 * \param[in] p        The pyramid to print
 * 
 */
template <t8_eclass_t eclass_T>
void
t8_sele_debug_print (const t8_standalone_element_t<eclass_T> *p);

template <t8_eclass_t eclass_T>
void
t8_sele_global_print (const t8_standalone_element_t<eclass_T> *p);

#include "t8_sele_bits_cxx.txx"

#endif /* T8_SELE_BITS_CXX_HXX */
