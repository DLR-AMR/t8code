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

#ifndef T8_DPYRAMID_BITS_H
#define T8_DPYRAMID_BITS_H

#include "t8_element.h"
#include "t8_dpyramid.h"

T8_EXTERN_C_BEGIN ();

/** Initialize a pyramid as the pyramid with a given global id in a uniform
 *  refinement of a given level. *
 * \param [in,out] p  Existing pyramid whose data will be filled.
 * \param [in] id     Index to be considered.
 * \param [in] level  level of uniform grid to be considered.
 */
void                t8_dpyramid_init_linear_id (t8_dpyramid_t * p, int level,
                                                uint64_t id);

/** Compute the level of a pyramid.
 * \param [in] p    Line whose pyramid is computed.
 * \return          The level of \a p.
 */
int                 t8_dpyramid_get_level (const t8_dpyramid_t * p);

/** Copy the data from source to dest
 * \param[in] source    The source-pyramid
 * \param[in/out] dest  The destination
 */
void                t8_dpyramid_copy (const t8_dpyramid_t * source,
                                      t8_dpyramid_t * dest);

/** Computes the linear position of a pyramid in an uniform grid.
 * \param [in] p  pyramid whose id will be computed.
 * \return Returns the linear position of this pyramid on a grid.
 */
t8_linearidx_t      t8_dpyramid_linear_id (const t8_dpyramid_t * p,
                                           int level);

/** Compute the childid-th child in Morton order of a pyramid.
 * \param [in] t    Input pyramid.
 * \param [in,out] childid The id of the child, 0..7 in Morton order.
 * \param [out] child  Existing pyramid whose data will be filled
 * 		    with the date of t's childid-th child.
 */
void                t8_dpyramid_child (const t8_dpyramid_t * elem,
                                       int child_id, t8_dpyramid_t * child);

/** Compute the children of a pyramid, array version
 * \param [in] p        Input pyramid
 * \param [in, out] c   Pointers to the computed children in Morton order
 */
void                t8_dpyramid_children (const t8_dpyramid_t * p, int length,
                                          t8_dpyramid_t ** c);

void                t8_dpyramid_children_at_face(const t8_dpyramid_t * p, int face,
                             t8_dpyramid_t * children[], int num_children,
                             int *child_indices);

int                 t8_dpyramid_face_child_face(const t8_dpyramid_t * p,
                                       int face, int face_child);

/** Compare two elements. returns negativ if p1 < p2, zero if p1 equals p2
 *  and positiv if p1 > p2.
 *  If p2 is a copy of p1 then the elements are equal.
 */
int                 t8_dpyramid_compare (const t8_dpyramid_t * p1,
                                         const t8_dpyramid_t * p2);

/** Checks if two pyramids have the same coordinates, same type and same level.
 * \param [in]  p   first input pyramid
 * \param [in]  q   second input pyramid
 * \return          0, if they are equal, 1 ow.
 */
int                 t8_dpyramid_is_equal (const t8_dpyramid_t * p,
                                          const t8_dpyramid_t * q);

/** Check wether a collection of 10 pyramids is a family in Morton order.
 * \param [in]  fam A collection of pyramids
 * \return      Nonzero if \a fam is a family of pyramids
 */
int                 t8_dpyramid_is_family (const t8_dpyramid_t ** fam);

/** Copmute whether a given pyramid shares a given face with its root tree.
 * \param [in]  p       The input pyramid
 * \param [in]  face    A face of \a p
 * \return              True, if \a is a subface of the pyramid root element.
 */
int                 t8_dpyramid_is_root_boundary (const t8_dpyramid_t * p,
                                                  int face);

int                 t8_dpyramid_face_neighbor_inside (const t8_dpyramid_t * p,
                                                      t8_dpyramid_t * neigh,
                                                      int face,
                                                      int *neigh_face);

/** Compute the position of the ancestor of this child at level \a level within
 * its siblings.
 * \param [in] p  pyramid to be considered.
 * \return Returns its child id in 0..9
 */
int                 t8_dpyramid_child_id (const t8_dpyramid_t * p);

/** Compute the position of the ancestor of this child at level \a level within
 * its siblings. The parent of p is known.
 * \param [in] p  pyramid to be considered.
 * \param [in] parent The parent of \a p.
 * \return Returns its child id in 0..9
 */
int                 t8_dpyramid_child_id_known_parent (const t8_dpyramid_t *
                                                       p,
                                                       t8_dpyramid_t *
                                                       parent);

int
t8_dpyramid_is_inside_pyra(const t8_dpyramid_t *p, const t8_dpyramid_t *check);

/** Check, if the input pyramid with shape of a tet is inside of a tetrahedron
 * \param [in] p  pyramid with tet shape
 * \return      0, if the pyramid is insed of a tetrahedron
 */
int                 t8_dpyramid_is_inside_tet (const t8_dpyramid_t * p,
                                               int level,
                                               t8_dpyramid_t * anc);

/** Check, if a tet of type 0 or 3 has a common face with its ancestor at a given level
 * \param [in]p  input pyramid
 * \level [in]  The level of the ancestor
 * \face [in] face  A face of \a p.
 * \return      false if they don't share a face, true otherwise
 */
int                 t8_dpyramid_tet_boundary (const t8_dpyramid_t * p, int face);

/** compute if a given element lies on the tree boundary and return the face number
 * of the tree face. If not the return value is arbitrary
 * \param [in] elem     pyramid
 * \param [in] face     a face of \a elem
 * \return              See discription
 */
int                 t8_dpyramid_tree_face(const t8_dpyramid_t * p , int face);

/** Compute the first descendant of a pyramid at a given level. This is the descendant of
 * the pyramid in a uniform level refinement that has the smallest id.
 * \param [in] p        pyramid whose descendant is computed.
 * \param [out] s       Existing pyramid whose data will be filled with the data
 *                      of \a p's first descendant on level \a level.
 * \param [in] level    The refinement level. Must be greater than \a p's refinement
 *                      level.
 */
void                t8_dpyramid_first_descendant (const t8_dpyramid_t * p,
                                                  t8_dpyramid_t * desc,
                                                  int level);

void                t8_dpyramid_first_descendant_face(const t8_dpyramid_t * p,
                                                      int face,
                                                      t8_dpyramid_t * first_desc,
                                                      int level);

/** Compute the last descendant of a pyramid at a given level. This is the descendant of
 * the pyramid in a uniform level refinement that has the largest id.
 * \param [in] p        pyramid whose descendant is computed.
 * \param [out] s       Existing pyramid whose data will be filled with the data
 *                      of \a p's last descendant on level \a level.
 * \param [in] level    The refinement level. Must be greater than \a p's refinement
 *                      level.
 */
void                t8_dpyramid_last_descendant (const t8_dpyramid_t * p,
                                                 t8_dpyramid_t * desc,
                                                 int level);

void                t8_dpyramid_last_descendant_face(const t8_dpyramid_t * p,
                                 int face, t8_dpyramid_t * last_desc,
                                 int level);

/** Compute the coordinates of a vertex of a pyramid.
 * \param [in] p    Input pyramid.
 * \param [out] coordinates An array of 3 t8_dpyramid_coord_t that
 * 		     will be filled with the coordinates of the vertex.
 * \param [in] vertex The number of the vertex.
 */
void                t8_dpyramid_compute_coords (const t8_dpyramid_t * p,
                                                int vertex, int coords[]);

/** Compute the pyramid-parent-type of a tetrahedron
 * \param [in] p        Input pyramid
 * \param [out] parent  The parent whose type has to be computed.
 */
void                t8_dpyramid_tetparent_type (const t8_dpyramid_t * p,
                                                t8_dpyramid_t * parent);

/** Compute the parent of a given pyramid
 * \param [in] p        Input pyramid.
 * \param [out] parent  The parent of \a p.
 */
void                t8_dpyramid_parent (const t8_dpyramid_t * p,
                                        t8_dpyramid_t * parent);

/**
 * Compute the number of corners of a pyramid. If pyramid has type less than 6,
 * it is actually a tetrahedron.
 * \param [in] p    Input pyramid.
 * \return          The number of corners of p.
 */
int                 t8_dpyramid_num_vertices (const t8_dpyramid_t * p);

/** Compute the number of children of p
 * \param [in] p    Input pyramid.
 * \return          The number of children of p.
 */
int                 t8_dpyramid_num_children (const t8_dpyramid_t * p);

/** Compute the number of siblings of p
 * \param [in] p    Input pyramid
 * \return          The number of siblings of p.
 */
int                 t8_dpyramid_num_siblings (const t8_dpyramid_t * p);

/** Return the number of faces of p
 * \param [in] p    Input pyramid
 * \return          The number of faces of p
 */
int                 t8_dpyramid_num_faces(const t8_dpyramid_t *p);


/** Returns the shape of the pyramid (pyramid or tetrahedron)
 * \param [in] p    Input pyramid.
 * \return          The eclass of the element
 */
t8_eclass_t         t8_dpyramid_shape (const t8_dpyramid_t * p);

/** Computes the successor of a pyramid in a uniform grid of level \a level.
 * \param [in] elem  pyramid whose id will be computed.
 * \param [in,out] s Existing pyramid whose data will be filled with the
 *                data of \a l's successor on level \a level.
 * \param [in] level level of uniform grid to be considered.
 */
void                t8_dpyramid_successor (const t8_dpyramid_t * elem,
                                           t8_dpyramid_t * s, int level);

T8_EXTERN_C_END ();

#endif /* T8_DPYRAMID_BITS_H */
