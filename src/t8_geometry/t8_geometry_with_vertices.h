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

/** \file t8_geometry_with_vertices.h
 * This header file provides a C interface for functions for the 
 * t8_geometry_with_vertices class.
 */

#ifndef T8_GEOMETRY_WITH_VERTICES_H
#define T8_GEOMETRY_WITH_VERTICES_H

#include <t8_cmesh.h>

T8_EXTERN_C_BEGIN ();

/** Set the vertex coordinates of a tree in the cmesh.
 * This is currently inefficient, since the vertices are duplicated for each
 * tree. Eventually this function will be replaced by a more efficient one.
 * It is not allowed to call this function after \ref t8_cmesh_commit.
 * The eclass of the tree has to be set before calling this function.
 * \param [in,out] cmesh        The cmesh to be updated.
 * \param [in]     gtree_id     The global number of the tree.
 * \param [in]     vertices     An array of 3 doubles per tree vertex.
 * \param [in]     num_vertices The number of verticess in \a vertices. Must
 *                              match the number of corners of the tree.
 */
void
t8_cmesh_set_tree_vertices (t8_cmesh_t cmesh, const t8_gloidx_t gtree_id, const double *vertices,
                            const int num_vertices);

T8_EXTERN_C_END ();

#endif /* !T8_GEOMETRY_WITH_VERTICES_H */
