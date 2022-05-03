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

/**
 * \file t8_cmesh_reader_helper.cxx
 * This file provides helper-functions to use when we read a mesh
 * comming from an external mesh-generator. 
 * 
 */

#include <t8_eclass.h>
#include <t8_cmesh/t8_cmesh_reader_helper.hxx>
#include <t8_cmesh.h>

T8_EXTERN_C_BEGIN ();

void
t8_cmesh_correct_volume (double *tree_vertices, t8_eclass_t eclass)
{
  /* The volume described is negative. We need to change vertices.
   * For tets we switch 0 and 3.
   * For prisms we switch 0 and 3, 1 and 4, 2 and 5.
   * For hexahedra we switch 0 and 4, 1 and 5, 2 and 6, 3 and 7.
   * For pyramids we switch 0 and 4 */
  double              temp;
  int                 num_switches = 0;
  int                 switch_indices[4] = { 0 };
  int                 iswitch;
  T8_ASSERT (t8_eclass_to_dimension[eclass] == 3);
  t8_debugf ("Correcting negative volume.\n");
  switch (eclass) {
  case T8_ECLASS_TET:
    /* We switch vertex 0 and vertex 3 */
    num_switches = 1;
    switch_indices[0] = 3;
    break;
  case T8_ECLASS_PRISM:
    num_switches = 3;
    switch_indices[0] = 3;
    switch_indices[1] = 4;
    switch_indices[2] = 5;
    break;
  case T8_ECLASS_HEX:
    num_switches = 4;
    switch_indices[0] = 4;
    switch_indices[1] = 5;
    switch_indices[2] = 6;
    switch_indices[3] = 7;
    break;
  case T8_ECLASS_PYRAMID:
    num_switches = 1;
    switch_indices[0] = 4;
    break;
  default:
    SC_ABORT_NOT_REACHED ();
  }

  for (iswitch = 0; iswitch < num_switches; ++iswitch) {
    /* We switch vertex 0 + iswitch and vertex switch_indices[iswitch] */
    for (int i = 0; i < 3; i++) {
      temp = tree_vertices[3 * iswitch + i];
      tree_vertices[3 * iswitch + i] =
        tree_vertices[3 * switch_indices[iswitch] + i];
      tree_vertices[3 * switch_indices[iswitch] + i] = temp;
    }
  }
  T8_ASSERT (!t8_cmesh_tree_vertices_negative_volume
             (eclass, tree_vertices, t8_eclass_num_vertices[eclass]));
}

const int           t8_vtk_cell_face_to_vertex_num[T8_ECLASS_COUNT][6][4] = {
  /*Vertex */
  {{0, -1, -1, -1},
   {-1, -1, -1, -1},
   {-1, -1, -1, -1},
   {-1, -1, -1, -1},
   {-1, -1, -1, -1},
   {-1, -1, -1, -1}},
  /*line */
  {{0, -1, -1, -1},
   {0, -1, -1, -1},
   {-1, -1, -1, -1},
   {-1, -1, -1, -1},
   {-1, -1, -1, -1},
   {-1, -1, -1, -1}},
  /*triangle */
  {{0, 1, -1, -1},
   {1, 2, -1, -1},
   {2, 0, -1, -1},
   {-1, -1, -1, -1},
   {-1, -1, -1, -1},
   {-1, -1, -1, -1}},
  /*quad */
  {{0, 1, -1, -1},
   {1, 3, -1, -1},
   {3, 2, -1, -1},
   {2, 0, -1, -1},
   {-1, -1, -1, -1},
   {-1, -1, -1, -1}},
  /*tetra */
  {{0, 1, 2, -1},
   {0, 1, 3, -1},
   {1, 2, 3, -1},
   {2, 0, 3, -1},
   {-1, -1, -1, -1},
   {-1, -1, -1, -1}},
  /*hex */
  {{0, 1, 2, 3},
   {0, 1, 5, 4},
   {1, 3, 7, 5},
   {3, 2, 6, 7},
   {2, 0, 4, 6},
   {4, 5, 6, 7}},
  /*prism/vedge */
  {{0, 1, 2, -1},
   {0, 1, 3, 4},
   {1, 2, 4, 5},
   {0, 2, 3, 5},
   {3, 4, 5, -1},
   {-1, -1, -1, -1}},
  /*pyramid */
  {{0, 1, 2, 3},
   {0, 1, 4, -1},
   {1, 2, 4, -1},
   {2, 3, 4, -1},
   {3, 0, 4, -1},
   {-1, -1, -1, -1}},
};

T8_EXTERN_C_END ();
