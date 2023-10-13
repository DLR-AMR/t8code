/*
This file is part of t8code.
t8code is a C library to manage a collection (a forest) of multiple
connected adaptive space-trees of general element classes in parallel.

Copyright (C) 2023 the developers

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
 * This test program tests if the face connectivity between two trees of the same
 * eclass is correct w.r.t. orientation. We compute the faceneighbor and the coordinates
 * of the corners of the touching faces. If they match for children of a the face the
 * test is passed. 
 */

#include <gtest/gtest.h>
#include <t8_cmesh.h>
#include <t8_cmesh_vtk_writer.h>
#include <t8_eclass.h>
#include <t8_schemes/t8_default/t8_default_cxx.hxx>
#include <t8_cmesh/t8_cmesh_examples.h>

class tree_face_neigh: public testing::TestWithParam<t8_eclass_t> {
 protected:
  void
  SetUp () override
  {
    eclass = GetParam ();
    if (eclass == T8_ECLASS_PRISM) {
      GTEST_SKIP ();
    }
    scheme = t8_scheme_new_default_cxx ();
  }
  void
  TearDown () override
  {
    if (eclass != T8_ECLASS_PRISM) {
      t8_scheme_cxx_unref (&scheme);
    }
  }
  t8_eclass_t eclass;
  t8_scheme_cxx_t *scheme;
};

TEST_P (tree_face_neigh, t8_forest_neigh_face_test)
{
  for (int face = 0; face < t8_eclass_num_faces[eclass]; face++) {
    const int face_type = t8_eclass_face_types[eclass][face];
    const int num_vertices = t8_eclass_num_vertices[face_type];
    for (int orientation = 0; orientation < num_vertices; orientation++) {
      char fileprefix[BUFSIZ];
      t8_cmesh_t cmesh = t8_cmesh_new_two_trees_face_orientation (eclass, face, orientation, sc_MPI_COMM_WORLD);

      snprintf (fileprefix, BUFSIZ, "%s_face_%i_%i", t8_eclass_to_string[eclass], face, orientation);
      t8_cmesh_vtk_write_file (cmesh, fileprefix, 1.0);
      t8_cmesh_destroy (&cmesh);
    }
  }
}

INSTANTIATE_TEST_SUITE_P (t8_gtest_tree_face_neigh, tree_face_neigh, testing::Range (T8_ECLASS_ZERO, T8_ECLASS_COUNT));