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

/** \file t8_eclass.c
 * Definition of C const integers using the value macros from \ref t8_eclass.h.
 */

#define KEEP_ECLASS_VALUE_DEFINITIONS
#include <t8_element/t8_eclass.h>
#undef KEEP_ECLASS_VALUE_DEFINITIONS

const int t8_eclass_to_dimension[T8_ECLASS_COUNT] = T8_ECLASS_TO_DIMENSION_VALUES;

const int t8_eclass_num_faces[T8_ECLASS_COUNT] = T8_ECLASS_NUM_FACES_VALUES;

const int t8_eclass_max_num_faces[T8_ECLASS_MAX_DIM + 1] = T8_ECLASS_MAX_NUM_FACES_VALUES;

const int t8_eclass_max_num_children[T8_ECLASS_COUNT] = T8_ECLASS_MAX_NUM_CHILDREN_VALUES;

const int t8_face_vertex_to_tree_vertex[T8_ECLASS_COUNT][T8_ECLASS_MAX_FACES][T8_ECLASS_MAX_CORNERS_2D]
  = T8_FACE_VERTEX_TO_TREE_VERTEX_VALUES;

const int t8_face_edge_to_tree_edge[T8_ECLASS_COUNT][T8_ECLASS_MAX_FACES][T8_ECLASS_MAX_EDGES_2D]
  = T8_FACE_EDGE_TO_TREE_EDGE_VALUES;

const int t8_face_to_edge_neighbor[T8_ECLASS_COUNT][T8_ECLASS_MAX_FACES][T8_ECLASS_MAX_CORNERS_2D]
  = T8_FACE_TO_EDGE_NEIGHBOR_VALUES;

const int t8_edge_vertex_to_tree_vertex[T8_ECLASS_COUNT][T8_ECLASS_MAX_EDGES][2] = T8_EDGE_VERTEX_TO_TREE_VERTEX_VALUES;

const int t8_edge_to_face[T8_ECLASS_COUNT][T8_ECLASS_MAX_EDGES][2] = T8_EDGE_TO_FACE_VALUES;

const int t8_eclass_face_orientation[T8_ECLASS_COUNT][T8_ECLASS_MAX_FACES] = T8_ECLASS_FACE_ORIENTATION_VALUES;

const int t8_reference_face_normal_tet[T8_ECLASS_MAX_FACES][3] = T8_REFERENCE_FACE_NORMAL_TET_VALUES;

const int t8_eclass_num_vertices[T8_ECLASS_COUNT] = T8_ECLASS_NUM_VERTICES_VALUES;

const int t8_eclass_num_edges[T8_ECLASS_COUNT] = T8_ECLASS_NUM_EDGES_VALUES;

const int t8_eclass_vtk_type[T8_ECLASS_COUNT] = T8_ECLASS_VTK_TYPE_VALUES;

const int t8_eclass_vtk_to_t8_corner_number[T8_ECLASS_COUNT][T8_ECLASS_MAX_CORNERS]
  = T8_ECLASS_VTK_TO_T8_CORNER_NUMBER_VALUES;

const int t8_eclass_t8_to_vtk_corner_number[T8_ECLASS_COUNT][T8_ECLASS_MAX_CORNERS]
  = T8_ECLASS_T8_TO_VTK_CORNER_NUMBER_VALUES;

const int t8_eclass_face_types[T8_ECLASS_COUNT][T8_ECLASS_MAX_FACES] = T8_ECLASS_FACE_TYPES_VALUES;

const int t8_eclass_boundary_count[T8_ECLASS_COUNT][T8_ECLASS_COUNT] = T8_ECLASS_BOUNDARY_COUNT_VALUES;

const char *t8_eclass_to_string[T8_ECLASS_INVALID] = T8_ECLASS_TO_STRING_VALUES;

int
t8_eclass_count_boundary (t8_eclass_t theclass, int min_dim, int *per_eclass)
{
  int sum = 0;
  for (int t = T8_ECLASS_ZERO; t < T8_ECLASS_COUNT; ++t) {
    if (t8_eclass_to_dimension[t] >= min_dim) {
      sum += (per_eclass[t] = t8_eclass_boundary_count[theclass][t]);
    }
    else {
      per_eclass[t] = 0;
    }
  }

  return sum;
}

/* Compares two eclasses within the order
 * Tri < Quad
 * Tet < Hex < Prism < Pyramid
 * Eclasses of different dimension are not allowed to be compared.
 */
int
t8_eclass_compare (t8_eclass_t eclass1, t8_eclass_t eclass2)
{
  int dim = t8_eclass_to_dimension[eclass1];
  T8_ASSERT (dim == t8_eclass_to_dimension[eclass2]);

  if (eclass1 == eclass2) {
    /* If both are equal return 0.
     * This also captures the case dim <= 1. */
    return 0;
  }
  else if (dim == 2) {
    /* Either eclass1 = tri and eclass2 = quad or the other way around. */
    return eclass1 == T8_ECLASS_TRIANGLE ? -1 : 1;
  }
  else {
    T8_ASSERT (dim == 3);
    switch (eclass1) {
    case T8_ECLASS_TET:
      return -1;
    case T8_ECLASS_HEX:
      return eclass2 == T8_ECLASS_TET ? 1 : -1;
    case T8_ECLASS_PRISM:
      return eclass2 == T8_ECLASS_PYRAMID ? -1 : 1;
    default:
      T8_ASSERT (eclass1 == T8_ECLASS_PYRAMID);
      return 1;
    }
  }
}

int
t8_eclass_is_valid (t8_eclass_t eclass)
{
  /* every eclass up to T8_ECLASS_COUNT is a valid class T8_ECLASS_COUNT
   * itself is invalid, every class higher than eclass count is considered
   * invalid.*/
  return eclass < T8_ECLASS_COUNT;
}
