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

/** \file t8_dline_bits.h
 * Definitions of line-specific functions.
 */

#ifndef T8_DLINE_BITS_H
#define T8_DLINE_BITS_H

#include <t8_element.h>
#include <t8_schemes/t8_default/t8_default_line/t8_dline.h>
#include <t8_schemes/t8_default/t8_default_vertex/t8_dvertex.h>

T8_EXTERN_C_BEGIN ();

/** Compute the level of a line.
 * \param [in] line    Line whose level is computed.
 * \return          The level of \a line.
 */
int
t8_dline_get_level (const t8_dline_t *line);

/** Copy all values from one line to another.
 * \param [in] line    The line to be copied.
 * \param [in,out] dest Existing line whose data will be filled with the data
 *                   of \a line.
 */
void
t8_dline_copy (const t8_dline_t *line, t8_dline_t *dest);

/** Compare two elements. returns negative if line1 < line2, zero if line1 equals line2
 *  and positive if line1 > line2.
 *  If line2 is a copy of line1 then the elements are equal.
 */
int
t8_dline_compare (const t8_dline_t *line1, const t8_dline_t *line2);

/** Check if two elements are equal.
* \param [in] line1  The first element.
* \param [in] line2  The second element.
* \return            1 if the elements are equal, 0 if they are not equal
*/
int
t8_dline_equal (const t8_dline_t *line1, const t8_dline_t *line2);

/** Compute the parent of a line.
 * \param [in]  line   The input line.
 * \param [in,out] parent Existing line whose data will be filled with the parent
 *                  data of \a line.
 */
void
t8_dline_parent (const t8_dline_t *line, t8_dline_t *parent);

/** Compute the ancestor of a line at a given level.
 * \param [in]  line   Input line.
 * \param [in]  level A smaller level than \a line.
 * \param [in,out] ancestor Existing line whose data will
 *                  be filled with the data of \a line's ancestor on
 *                  level \a level.
 * \note The line \a ancestor may point to the same line as \a line.
 */
void
t8_dline_ancestor (const t8_dline_t *line, int level, t8_dline_t *ancestor);

/** Compute the childid-th child in Morton order of a line.
 * \param [in] line    Input Line.
 * \param [in] childid The id of the child, 0 or 1, in Morton order.
 * \param [in,out] child  Existing Line whose data will be filled
 * 		    with the date of l's childid-th child.
 */
void
t8_dline_child (const t8_dline_t *line, int childid, t8_dline_t *child);

/** Compute the face neighbor of a line.
 * \param [in]     line      Input line.
 * \param [in,out] neigh  Existing line whose data will be filled.
 * \param [in]     face   The face across which to generate the neighbor.
 * \param [out]    dual_face If not NULL, the face number as seen from \a neigh
 *                        is stored.
 * \note \a line may point to the same line as \a neigh.
 */
void
t8_dline_face_neighbour (const t8_dline_t *line, t8_dline_t *neigh, int face, int *dual_face);

/** Computes the nearest common ancestor of two lines in the same tree.
 * \param [in]     line1 First input line.
 * \param [in]     line2 Second input line.
 * \param [in,out] nca Existing line whose data will be filled.
 * \note \a line1, \a line2, \a nca may point to the same line.
 */
void
t8_dline_nearest_common_ancestor (const t8_dline_t *line1, const t8_dline_t *line2, t8_dline_t *nca);

/** Compute the position of the ancestor of this child at level \a level within
 * its siblings.
 * \param [in] line  line to be considered.
 * \param [in] level level to be considered.
 * \return Returns its child id 0 or 1.
 */
int
t8_dline_ancestor_id (const t8_dline_t *line, int level);

/** Given a face of a line return the face number
 * of the parent of the line that matches the line's face. Or return -1 if
 * no face of the parent matches the face.

 * \param [in]  line       The line.
 * \param [in]  face    The number of the face.
 * \return              If \a face of \a line is also a face of \a line's parent,
 *                      the face number of this face. Otherwise -1.
 */
int
t8_dline_face_parent_face (const t8_dline_t *line, int face);

/** Compute the position of the ancestor of this child at level \a level within
 * its siblings.
 * \param [in] line  line to be considered.
 * \return Returns its child id in 0,1
 */
int
t8_dline_child_id (const t8_dline_t *line);

/** Compute the 2 children of a line, array version.
 * \param [in]     line  Input line.
 * \param [in,out] c  Pointers to the 2 computed children in Morton order.
 *                    t may point to the same quadrant as c[0].
 */
void
t8_dline_childrenpv (const t8_dline_t *line, t8_dline_t *c[T8_DLINE_CHILDREN]);

/** Check whether a collection of two lines is a family in Morton order.
 * \param [in]     f  An array of two lines.
 * \return            Nonzero if \a f is a family of lines.
 */
int
t8_dline_is_familypv (const t8_dline_t *f[]);

/** Compute whether a given line shares a given face with its root tree.
 * \param [in] line       The input line.
 * \param [in] face     A face of \a line.
 * \return              True if \a face is a subface of the line's root element.
 */
int
t8_dline_is_root_boundary (const t8_dline_t *line, int face);

/** Test if a line lies inside of the root line,
 *  that is the line of level 0, anchor node (0,0)
 *  \param [in]     line Input line.
 *  \return true    If \a line lies inside of the root line.
 */
int
t8_dline_is_inside_root (const t8_dline_t *line);

/** Initialize a line as the line with a given global id in a uniform
 *  refinement of a given level. *
 * \param [in,out] line  Existing line whose data will be filled.
 * \param [in] id     Index to be considered.
 * \param [in] level  level of uniform grid to be considered.
 */
void
t8_dline_init_linear_id (t8_dline_t *line, int level, t8_linearidx_t id);

/** Computes the successor of a line in a uniform grid of level \a level.
 * \param [in] line  line whose id will be computed.
 * \param [in,out] succ Existing line whose data will be filled with the
 *                data of \a line's successor on level \a level.
 * \param [in] level level of uniform grid to be considered.
 */
void
t8_dline_successor (const t8_dline_t *line, t8_dline_t *succ, int level);

/** Suppose we have two trees that share a common face f.
 *  Given a Line e that is a subface of f in one of the trees
 *  and given the orientation of the tree connection, construct the face
 *  Line of the respective tree neighbor that logically coincides with e
 *  but lies in the coordinate system of the neighbor tree.
 *  \param [in] line1     The face element.
 *  \param [in,out] line2 On return the face element \a line1 with respect
 *                        to the coordinate system of the other tree.
 *  \param [in] orientation The orientation of the tree-tree connection.
 *                        0 if vertex 0 of face 0 coincides with vertex 0 of face 1.
 *                        1 if vertex 0 of face 0 coincides with vertex 1 of face 1.
 */
void
t8_dline_transform_face (const t8_dline_t *line1, t8_dline_t *line2, int orientation);

/** Given a vertex at the boundary of a line at a root tree boundary,
 *  construct the line from it.
 * \param [in] face     The face element (vertex).
 * \param [in] root_face  The index of the face of the tree.
 * \param [out] line    The line that has \a face as face element at face \a root_face
 * \return              The face number pf \a line that coincides with \a face,
 *                      thus \a root_face is returned.
 */
int
t8_dline_extrude_face (const t8_dvertex_t *face, int root_face, t8_dline_t *line);

/** Compute the first descendant of a line at a given level. This is the descendant of
 * the line in a uniform level refinement that has the smallest id.
 * \param [in] line       Line whose descendant is computed.
 * \param [out] desc       Existing line whose data will be filled with the data
 *                      of \a line's first descendant on level \a level.
 * \param [in] level    The refinement level. Must be greater than \a line's refinement
 *                      level.
 */
void
t8_dline_first_descendant (const t8_dline_t *line, t8_dline_t *desc, int level);

/** Compute the last descendant of a line at a given level. This is the descendant of
 * the line in a uniform level refinement that has the largest id.
 * \param [in] line        Line whose descendant is computed.
 * \param [out] desc       Existing line whose data will be filled with the data
 *                      of \a line's last descendant on level \a level.
 * \param [in] level    The refinement level. Must be greater than \a line's refinement
 *                      level.
 */
void
t8_dline_last_descendant (const t8_dline_t *line, t8_dline_t *desc, int level);

/** Compute the first or second vertex of a line.
 * \param [in] line     Line whose vertex is computed.
 * \param [in] vertex   The number of the vertex of \a line
 * \param [out] coords   The coordinates of the computed vertex
 */
void
t8_dline_vertex_integer_coords (const t8_dline_t *line, const int vertex, int coords[]);

/** Compute the coordinates of a vertex of a line when the 
 * tree (level 0 line) is embedded in [0,1]^1.
 * \param [in] line    Input line.
 * \param [in] vertex The number of the vertex.
 * \param [out] coordinates An array of 1 double that
 * 		     will be filled with the reference coordinates of the vertex.
 */
void
t8_dline_vertex_ref_coords (const t8_dline_t *line, const int vertex, double coordinates[1]);

/** Convert points in the reference space of a line element to points in the
 *  reference space of the tree (level 0) embedded in [0,1]^1.
 * \param [in]  line       Input line.
 * \param [in]  ref_coords The reference coordinates in the line
 *                         (\a num_coords times \f$ [0,1]^1 \f$)
 * \param [in]  num_coords Number of coordinates to evaluate
 * \param [in]  skip_coords Only used for batch computation of prisms.
 *                          In all other cases 0.
 *                          Skip coordinates in the \a ref_coords and
 *                          \a out_coords array.
 * \param [out] out_coords An array of \a num_coords x 1 x double that
 * 		                     will be filled with the reference coordinates
 *                         of the points on the line.
 */
void
t8_dline_compute_reference_coords (const t8_dline_t *line, const double *ref_coords, const size_t num_coords,
                                   const size_t skip_coords, double *out_coords);

/** Computes the linear position of a line in an uniform grid.
 * \param [in] line  Pointer to a line element whose id will be computed.
 * \param [in] level Refinement level of the line element.
 * \return Returns the linear position of this line on a grid.
 */
t8_linearidx_t
t8_dline_linear_id (const t8_dline_t *line, int level);

/** Query whether all entries of a line are in valid ranges.
 * \param [in] line  line to be considered.
 * \return        True, if \a line is a valid line and it is safe to call any
 *                function in this file on \a line.
 *                False otherwise.
 */
int
t8_dline_is_valid (const t8_dline_t *line);

/** Set default values for a line, such that it passes \ref t8_dline_is_valid.
 * \param [in] line  line to be initialized
 */
void
t8_dline_init (t8_dline_t *line);

T8_EXTERN_C_END ();

#endif /* T8_DLINE_BITS_H */
