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
#include <t8_cmesh/t8_cmesh_trees.h>
#include <t8_cmesh_vtk_writer.h>
#include <t8_eclass.h>
#include <t8_schemes/t8_default/t8_default_cxx.hxx>
#include <t8_cmesh/t8_cmesh_examples.h>
#include <t8_forest/t8_forest_general.h>
#include <t8_forest/t8_forest_geometrical.h>

class tree_face_neigh: public testing::TestWithParam<t8_eclass_t> {
 protected:
  void
  SetUp () override
  {
    eclass = GetParam ();
    if (eclass == T8_ECLASS_PRISM) {
      GTEST_SKIP ();
    }
  }
  void
  TearDown () override
  {
  }
  t8_eclass_t eclass;
  t8_scheme_cxx_t *scheme;
  sc_MPI_Comm comm = sc_MPI_COMM_WORLD;
};

TEST_P (tree_face_neigh, t8_forest_neigh_face_test)
{
  /*Iterate over all faces of the tree */
  for (int face = 0; face < t8_eclass_num_faces[eclass]; face++) {
    const int face_type = t8_eclass_face_types[eclass][face];
    const int num_vertices = t8_eclass_num_vertices[face_type];
    /* Check all possible orientations */
    for (int orientation = 0; orientation < num_vertices; orientation++) {
      t8_debugf ("[D] ------------------Test: %s face %i orientation %i -----------------\n",
                 t8_eclass_to_string[eclass], face, orientation);
      char fileprefix[BUFSIZ];
      t8_cmesh_t cmesh = t8_cmesh_new_two_trees_face_orientation (eclass, face, orientation, sc_MPI_COMM_WORLD);
      scheme = t8_scheme_new_default_cxx ();
      t8_forest_t forest = t8_forest_new_uniform (cmesh, scheme, 0, 0, comm);
      t8_element_t *elem = t8_forest_get_element_in_tree (forest, 0, 0);
      /**/
      const int face_class = t8_eclass_face_types[eclass][face];
      t8_eclass_scheme_c *ts = scheme->eclass_schemes[eclass];
      t8_eclass_scheme_c *face_scheme = scheme->eclass_schemes[face_class];
      t8_element_t *t0_face;
      t8_element_t *t1_face;
      const int num_face_children = ts->t8_element_num_face_children (elem, face);
      face_scheme->t8_element_new (1, &t0_face);
      face_scheme->t8_element_new (1, &t1_face);
      t8_element_t **children = T8_ALLOC (t8_element_t *, num_face_children);
      t8_element_t *face_elem;
      ts->t8_element_new (num_face_children, children);
      ts->t8_element_new (1, &face_elem);
      ts->t8_element_children_at_face (elem, face, children, num_face_children, NULL);
      const int num_face_vertices = t8_eclass_num_vertices[face_class];

      for (int ichild = 0; ichild < num_face_children; ichild++) {
        /*Get the child at the face. */
        const int child_face = ts->t8_element_face_child_face (elem, face, ichild);
        /*Construct the face of the child touching the boundary*/
        ts->t8_element_boundary_face (children[ichild], child_face, t0_face, face_scheme);

        /* Transfer the face onto tree 1*/
        t8_locidx_t *face_neighbor;
        int8_t *ttf;
        (void) t8_cmesh_trees_get_tree_ext (cmesh->trees, 0, &face_neighbor, &ttf);
        /* Compute the local id of the face neighbor tree. */
        const t8_locidx_t lcneigh_id = face_neighbor[face];
        EXPECT_EQ (lcneigh_id, 1);
        /* F is needed to compute the neighbor face number and the orientation.
        * tree_neigh_face = ttf % F
        * or = ttf / F
        */
        const int F = t8_eclass_max_num_faces[cmesh->dimension];
        /* compute the neighbor face */
        const int tree_neigh_face = ttf[face] % F;
        /* eclasses of the trees are the same, so it is sufficient to compare the face-numbers */
        const int is_smaller = face <= tree_neigh_face;
        const int sign
          = t8_eclass_face_orientation[eclass][face] == t8_eclass_face_orientation[eclass][tree_neigh_face];
        /* Actual face transfer */
        face_scheme->t8_element_transform_face (t0_face, t1_face, ttf[face] / F, sign, is_smaller);
        /* Extrude face onto tree 1*/
        const int neigh_face = ts->t8_element_extrude_face (t1_face, face_scheme, face_elem, tree_neigh_face);
        for (int ivertex = 0; ivertex < num_face_vertices; ivertex++) {
          /* Compute the coordinates of all vertices touching the face of both elements*/
          int t0_face_vertex = t8_face_vertex_to_tree_vertex[eclass][child_face][ivertex];
          int t1_face_vertex
            = t8_face_vertex_to_tree_vertex[eclass][neigh_face][(ivertex + orientation) % num_face_vertices];
          double t0_face_vertex_coords[3];
          double t1_face_vertex_coords[3];
          t8_forest_element_coordinate (forest, 0, children[ichild], t0_face_vertex, t0_face_vertex_coords);
          t8_forest_element_coordinate (forest, 1, face_elem, t1_face_vertex, t1_face_vertex_coords);
          /* If everything is set correct the coordinates should be equal. */
          for (int icoord = 0; icoord < 3; icoord++) {
            EXPECT_NEAR (t0_face_vertex_coords[icoord], t1_face_vertex_coords[icoord], T8_PRECISION_EPS);
          }
        }
      }
      face_scheme->t8_element_destroy (1, &t0_face);
      face_scheme->t8_element_destroy (1, &t1_face);
      ts->t8_element_destroy (num_face_children, children);
      ts->t8_element_destroy (1, &face_elem);
      T8_FREE (children);

      //snprintf (fileprefix, BUFSIZ, "%s_face_%i_%i", t8_eclass_to_string[eclass], face, orientation);
      //t8_cmesh_vtk_write_file (cmesh, fileprefix, 1.0);
      t8_forest_unref (&forest);
    }
  }
}

INSTANTIATE_TEST_SUITE_P (t8_gtest_tree_face_neigh, tree_face_neigh, testing::Range (T8_ECLASS_ZERO, T8_ECLASS_COUNT));