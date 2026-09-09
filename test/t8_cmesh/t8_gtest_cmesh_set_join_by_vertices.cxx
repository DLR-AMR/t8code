/*
  This file is part of t8code.
  t8code is a C library to manage a collection (a forest) of multiple
  connected adaptive space-trees of general element types in parallel.

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

#include <gtest/gtest.h>
#include <t8.h>
#include <t8_eclass/t8_eclass.h>
#include <t8_cmesh/t8_cmesh.h>
#include <t8_cmesh/t8_cmesh_internal/t8_cmesh_types.h>
#include <t8_cmesh/t8_cmesh_examples.h>
#include <t8_cmesh/t8_cmesh_internal/t8_cmesh_helpers.hxx>
#include <t8_cmesh/t8_cmesh_vertex_connectivity/t8_cmesh_vertex_connectivity.h>
#include <t8_geometry/t8_geometry_with_vertices.h>
#include "test/t8_cmesh_generator/t8_cmesh_example_sets.hxx"
#include <t8_geometry/t8_geometry_implementations/t8_geometry_linear.hxx>
#include <test/t8_gtest_macros.hxx>
#include "t8_test_data_dir.h"
#include <string>

#include <p8est_geometry.h>

/* In this file we test `t8_cmesh_join_by_vertices` with a lot of example meshes
 * provided by t8code and p4est. The general idea is as follows: We build a second, fresh,
 * uncommitted cmesh with the same tree classes and vertices as a given (already committed)
 * cmesh, but with no explicit joins, so that `t8_cmesh_join_by_vertices` has to
 * derive all face connectivity itself (the coordinate-matching path, since no global vertex
 * ids were ever set on the fresh cmesh). We then compare the result against the information
 * given by `t8_cmesh_get_face_neighbor` on the original cmesh.
 */

static void
test_with_cmesh (t8_cmesh_t cmesh)
{
  const t8_locidx_t ntrees = t8_cmesh_get_num_local_trees (cmesh);

  /* Build a fresh, uncommitted cmesh with the same tree classes and vertices as `cmesh`, but
   * no explicit joins. */
  t8_cmesh_t cmesh2;
  t8_cmesh_init (&cmesh2);
  t8_cmesh_set_dimension (cmesh2, t8_cmesh_get_dimension (cmesh));
  for (t8_locidx_t itree = 0; itree < ntrees; itree++) {
    const t8_eclass_t eclass = t8_cmesh_get_tree_class (cmesh, itree);
    t8_cmesh_set_tree_class (cmesh2, itree, eclass);
    t8_cmesh_set_tree_vertices (cmesh2, itree, t8_cmesh_get_tree_vertices (cmesh, itree),
                                t8_eclass_num_vertices[eclass]);
  }

  /* Compute face connectivity. */
  t8_cmesh_join_by_vertices (cmesh2);
  t8_cmesh_commit (cmesh2, sc_MPI_COMM_WORLD);

  /* Compare results with `t8_cmesh_get_face_neighbor`. */
  for (int this_itree = 0; this_itree < ntrees; this_itree++) {
    const t8_eclass_t this_eclass = t8_cmesh_get_tree_class (cmesh, this_itree);
    const int this_nfaces = t8_eclass_num_faces[this_eclass];

    for (int this_iface = 0; this_iface < this_nfaces; this_iface++) {
      int conn_dual_iface;
      int conn_orientation;
      const t8_locidx_t conn_dual_itree
        = t8_cmesh_get_face_neighbor (cmesh2, this_itree, this_iface, &conn_dual_iface, &conn_orientation);

      int cmesh_dual_iface;
      int cmesh_orientation;

      t8_locidx_t cmesh_dual_itree
        = t8_cmesh_get_face_neighbor (cmesh, this_itree, this_iface, &cmesh_dual_iface, &cmesh_orientation);

      /* Here we check for a connected domain boundary (e.g. periodic
       * boundary). In this case we skip the test. */
      if (cmesh_dual_itree > -1) {
        const t8_eclass_t dual_eclass = t8_cmesh_get_tree_class (cmesh, cmesh_dual_itree);

        const int this_nface_verts = t8_eclass_num_vertices[t8_eclass_face_types[this_eclass][this_iface]];
        const int dual_nface_verts = t8_eclass_num_vertices[t8_eclass_face_types[dual_eclass][cmesh_dual_iface]];

        const double *this_vertices = t8_cmesh_get_tree_vertices (cmesh, this_itree);
        const double *dual_vertices = t8_cmesh_get_tree_vertices (cmesh, cmesh_dual_itree);

        int match_count = 0;
        for (int this_iface_vert = 0; this_iface_vert < this_nface_verts; this_iface_vert++) {
          const int this_ivert = t8_face_vertex_to_tree_vertex[this_eclass][this_iface][this_iface_vert];

          for (int dual_iface_vert = 0; dual_iface_vert < dual_nface_verts; dual_iface_vert++) {
            const int dual_ivert = t8_face_vertex_to_tree_vertex[dual_eclass][cmesh_dual_iface][dual_iface_vert];

            int match_count_per_coord = 0;
            for (int icoord = 0; icoord < T8_ECLASS_MAX_DIM; icoord++) {
              const double this_face_vert
                = this_vertices[T8_2D_TO_1D (this_nface_verts, T8_ECLASS_MAX_DIM, this_ivert, icoord)];
              const double dual_face_vert
                = dual_vertices[T8_2D_TO_1D (dual_nface_verts, T8_ECLASS_MAX_DIM, dual_ivert, icoord)];

              if (fabs (this_face_vert - dual_face_vert) < 10 * T8_PRECISION_EPS) {
                match_count_per_coord++;
              }
              else {
                break;
              }
            }
            if (match_count_per_coord == T8_ECLASS_MAX_DIM) {
              match_count++;
              continue;
            }
          }
        }

        if (match_count < this_nface_verts) {
          continue;
        }
      }

      /* Here we ignore found connections which are not stored in the `cmesh` object. */
      if (cmesh_dual_itree > -1) {
        EXPECT_EQ (conn_dual_itree, cmesh_dual_itree) << "Neighboring trees do not match.";
        EXPECT_EQ (conn_dual_iface, cmesh_dual_iface) << "Dual faces do not match.";
        EXPECT_EQ (conn_orientation, cmesh_orientation) << "Face orientations do not match.";
      }
    }
  }

  t8_cmesh_destroy (&cmesh2);
}

TEST (t8_cmesh_set_join_by_vertices, test_cmesh_set_join_by_vertices)
{
  /* Some defaults. */
  sc_MPI_Comm comm = sc_MPI_COMM_WORLD;
  const int do_partition = 0;

  /*
   * Tests from `t8code` that are not included in AllCmeshsParam.
   */

  {
    const int periodic = 0;
    t8_cmesh_t cmesh;
    t8_cmesh_init (&cmesh);
    t8_cmesh_new_hypercube_hybrid (cmesh, sc_MPI_COMM_WORLD, periodic);
    test_with_cmesh (cmesh);
    t8_cmesh_destroy (&cmesh);
  }

  {  // Test additionally to AllCmeshsParams with interesting boundary_coords.
    const double boundary_coords[24] = { 1, 0, 0, 4, 0, 0, 0, 6, 0, 5, 5, 0, -1, -2, 8, 9, 0, 10, 0, 8, 9, 10, 10, 10 };

    t8_eclass_t eclass = T8_ECLASS_HEX;
    const int use_axis_aligned = 0;
    t8_cmesh_t cmesh;
    t8_cmesh_init (&cmesh);
    t8_cmesh_new_hypercube_pad (cmesh, eclass, comm, boundary_coords, 2, 2, 2, use_axis_aligned);
    test_with_cmesh (cmesh);
    t8_cmesh_destroy (&cmesh);
  }

  {
    t8_cmesh_t cmesh;
    t8_cmesh_init (&cmesh);
    t8_cmesh_new_brick_2d (cmesh, 3, 4, 1, 1, comm);
    test_with_cmesh (cmesh);
    t8_cmesh_destroy (&cmesh);
  }

  {
    t8_cmesh_t cmesh;
    t8_cmesh_init (&cmesh);
    t8_cmesh_new_brick_3d (cmesh, 3, 4, 5, 1, 1, 1, comm);
    test_with_cmesh (cmesh);
    t8_cmesh_destroy (&cmesh);
  }

  /* {
   *   // Side note: This test does not work.
   *   // Problem: Crashes with memory error.
   *   // Reason: There are no tree vertices set. Thus it cannot work.
   *   int num_trees = 1;
   *   t8_eclass_t eclass = T8_ECLASS_HEX;
   *   t8_cmesh_t  cmesh = t8_cmesh_new_bigmesh (eclass, num_trees,
   *                                         sc_MPI_COMM_WORLD);
   *   test_with_cmesh(cmesh);
   *   t8_cmesh_destroy (&cmesh);
   * }
   */

  /*
   * Tests with 2D and 3D example meshes from `p4est`.
   */

  /* Make sure that p4est is properly initialized. If not, do it here. */
  if (!sc_package_is_registered (p4est_package_id)) {
    p4est_init (NULL, SC_LP_ESSENTIAL);
  }

  {
    p4est_connectivity_t *p4_conn = p4est_connectivity_new_disk_nonperiodic ();
    t8_cmesh_t cmesh;
    t8_cmesh_init (&cmesh);
    t8_cmesh_new_from_p4est (cmesh, p4_conn, comm, do_partition);
    p4est_connectivity_destroy (p4_conn);
    test_with_cmesh (cmesh);
    t8_cmesh_destroy (&cmesh);
  }

  {
    p4est_connectivity_t *p4_conn = p4est_connectivity_new_brick (3, 3, 0, 1);
    t8_cmesh_t cmesh;
    t8_cmesh_init (&cmesh);
    t8_cmesh_new_from_p4est (cmesh, p4_conn, comm, do_partition);
    p4est_connectivity_destroy (p4_conn);
    test_with_cmesh (cmesh);
    t8_cmesh_destroy (&cmesh);
  }

  {
    p4est_connectivity_t *p4_conn = p4est_connectivity_new_icosahedron ();
    t8_cmesh_t cmesh;
    t8_cmesh_init (&cmesh);
    t8_cmesh_new_from_p4est (cmesh, p4_conn, comm, do_partition);
    p4est_connectivity_destroy (p4_conn);
    test_with_cmesh (cmesh);
    t8_cmesh_destroy (&cmesh);
  }

  {
    p4est_connectivity_t *p4_conn = p4est_connectivity_new_star ();
    t8_cmesh_t cmesh;
    t8_cmesh_init (&cmesh);
    t8_cmesh_new_from_p4est (cmesh, p4_conn, comm, do_partition);
    p4est_connectivity_destroy (p4_conn);
    test_with_cmesh (cmesh);
    t8_cmesh_destroy (&cmesh);
  }

  {
    p4est_connectivity_t *p4_conn = p4est_connectivity_new_moebius ();
    t8_cmesh_t cmesh;
    t8_cmesh_init (&cmesh);
    t8_cmesh_new_from_p4est (cmesh, p4_conn, comm, do_partition);
    p4est_connectivity_destroy (p4_conn);
    test_with_cmesh (cmesh);
    t8_cmesh_destroy (&cmesh);
  }

  {
    p4est_connectivity_t *p4_conn = p4est_connectivity_new_pillow ();
    t8_cmesh_t cmesh;
    t8_cmesh_init (&cmesh);
    t8_cmesh_new_from_p4est (cmesh, p4_conn, comm, do_partition);
    p4est_connectivity_destroy (p4_conn);
    test_with_cmesh (cmesh);
    t8_cmesh_destroy (&cmesh);
  }

  {
    p4est_connectivity_t *p4_conn = p4est_connectivity_new_corner ();
    t8_cmesh_t cmesh;
    t8_cmesh_init (&cmesh);
    t8_cmesh_new_from_p4est (cmesh, p4_conn, comm, do_partition);
    p4est_connectivity_destroy (p4_conn);
    test_with_cmesh (cmesh);
    t8_cmesh_destroy (&cmesh);
  }

  {
    p8est_connectivity_t *p8_conn = p8est_connectivity_new_brick (3, 3, 3, 0, 0, 0);
    t8_cmesh_t cmesh;
    t8_cmesh_init (&cmesh);
    t8_cmesh_new_from_p8est (cmesh, p8_conn, comm, do_partition);
    test_with_cmesh (cmesh);
    p8est_connectivity_destroy (p8_conn);
    t8_cmesh_destroy (&cmesh);

    /* Note: Other p8est connectivity examples, like
     *   p8est_connectivity_t *p8_conn = p8est_connectivity_new_shell ();
     *   p8est_connectivity_t *p8_conn = p8est_connectivity_new_sphere ();
     * cannot work since the vertices are net set properly. This has to be fixed in `p4est`. */
  }

  {
    auto testfile = t8_test_data_dir / "test_cube_unstructured_1.inp";
    p8est_connectivity_t *p8_conn = p8est_connectivity_read_inp (testfile.c_str ());
    t8_cmesh_t cmesh;
    t8_cmesh_init (&cmesh);
    t8_cmesh_new_from_p8est (cmesh, p8_conn, comm, do_partition);
    test_with_cmesh (cmesh);
    p8est_connectivity_destroy (p8_conn);
    t8_cmesh_destroy (&cmesh);
  }

  {
    auto testfile = t8_test_data_dir / "test_cube_unstructured_2.inp";
    p8est_connectivity_t *p8_conn = p8est_connectivity_read_inp (testfile.c_str ());
    t8_cmesh_t cmesh;
    t8_cmesh_init (&cmesh);
    t8_cmesh_new_from_p8est (cmesh, p8_conn, comm, do_partition);
    test_with_cmesh (cmesh);
    p8est_connectivity_destroy (p8_conn);
    t8_cmesh_destroy (&cmesh);
  }
}

struct t8_cmesh_set_join_by_vertices_class: public testing::TestWithParam<cmesh_example_base *>
{
 protected:
  void
  SetUp () override
  {
    size_t found = GetParam ()->name.find (std::string ("bigmesh"));
    if (found != std::string::npos) {
      /* skip bigmeshes as they do not have vertices from which to build the connectivity */
      GTEST_SKIP ();
    }

    cmesh = GetParam ()->cmesh_create ();
    if (cmesh->set_partition) {
      /* skip partitioned cmeshes, as they do not have all necessary vertex information for the neighbors */
      GTEST_SKIP ();
    }
  }

  void
  TearDown () override
  {
    if (cmesh != NULL) {
      t8_cmesh_destroy (&cmesh);
    }
  }

  int cmesh_id;
  t8_cmesh_t cmesh = NULL;
};

TEST_P (t8_cmesh_set_join_by_vertices_class, test_cmesh_set_join_by_vertices_parameterized)
{
  test_with_cmesh (cmesh);
}

/* Test all cmeshes over all different inputs we get through their id */
INSTANTIATE_TEST_SUITE_P (t8_cmesh_set_join_by_vertices, t8_cmesh_set_join_by_vertices_class, AllCmeshsParam,
                          pretty_print_base_example);

/* This test exercises the id-based fast path of `t8_cmesh_join_by_vertices`, which the tests
 * above never reach (they always build a cmesh with no global vertex ids set, so the
 * coordinate-fallback path is taken instead). Requesting vertex connectivity makes
 * `t8_cmesh_commit` derive global vertex ids first; the subsequent automatic-face-joining step
 * then finds those ids in the stash and matches faces by id instead of by coordinate. */
TEST (t8_cmesh_set_join_by_vertices, test_automatic_face_joining_with_ids)
{
  /* Two quads sharing one edge: tree 0 spans [0,1]x[0,1], tree 1 spans [1,2]x[0,1]. */
  t8_cmesh_t cmesh;
  t8_cmesh_init (&cmesh);

  t8_cmesh_set_tree_class (cmesh, 0, T8_ECLASS_QUAD);
  t8_cmesh_set_tree_class (cmesh, 1, T8_ECLASS_QUAD);

  const double vertices0[12] = { 0, 0, 0, 1, 0, 0, 0, 1, 0, 1, 1, 0 };
  const double vertices1[12] = { 1, 0, 0, 2, 0, 0, 1, 1, 0, 2, 1, 0 };
  t8_cmesh_set_tree_vertices (cmesh, 0, vertices0, 4);
  t8_cmesh_set_tree_vertices (cmesh, 1, vertices1, 4);

  t8_cmesh_enable_vertex_conn (cmesh);
  t8_cmesh_activate_automatic_face_joining (cmesh);
  t8_cmesh_commit (cmesh, sc_MPI_COMM_WORLD);

  /* Tree 0's face 1 (x=1 side, vertices 1 and 3) connects to tree 1's face 0 (x=1 side, vertices 0 and 2). */
  int dual_face, orientation;
  const t8_locidx_t dual_tree = t8_cmesh_get_face_neighbor (cmesh, 0, 1, &dual_face, &orientation);
  EXPECT_EQ (dual_tree, 1);
  EXPECT_EQ (dual_face, 0);

  /* Global vertex ids should have been derived automatically: 6 unique corners total. */
  EXPECT_EQ (t8_cmesh_get_num_global_vertices (cmesh), 6);
  EXPECT_EQ (t8_cmesh_get_global_vertex_of_tree (cmesh, 0, 1), t8_cmesh_get_global_vertex_of_tree (cmesh, 1, 0));
  EXPECT_EQ (t8_cmesh_get_global_vertex_of_tree (cmesh, 0, 3), t8_cmesh_get_global_vertex_of_tree (cmesh, 1, 2));

  t8_cmesh_destroy (&cmesh);
}

/* Activating only automatic face joining (without requesting vertex connectivity) must still
 * join faces correctly via the coordinate fallback, but must not expose any global vertex ids on
 * the committed cmesh: the two features are independent, and the id-based speedup taken
 * internally by `t8_cmesh_join_by_vertices` is not a user-visible side effect. */
TEST (t8_cmesh_set_join_by_vertices, test_automatic_face_joining_without_vertex_conn)
{
  /* Two quads sharing one edge: tree 0 spans [0,1]x[0,1], tree 1 spans [1,2]x[0,1]. */
  t8_cmesh_t cmesh;
  t8_cmesh_init (&cmesh);

  t8_cmesh_set_tree_class (cmesh, 0, T8_ECLASS_QUAD);
  t8_cmesh_set_tree_class (cmesh, 1, T8_ECLASS_QUAD);

  const double vertices0[12] = { 0, 0, 0, 1, 0, 0, 0, 1, 0, 1, 1, 0 };
  const double vertices1[12] = { 1, 0, 0, 2, 0, 0, 1, 1, 0, 2, 1, 0 };
  t8_cmesh_set_tree_vertices (cmesh, 0, vertices0, 4);
  t8_cmesh_set_tree_vertices (cmesh, 1, vertices1, 4);

  t8_cmesh_activate_automatic_face_joining (cmesh);
  t8_cmesh_commit (cmesh, sc_MPI_COMM_WORLD);

  int dual_face, orientation;
  const t8_locidx_t dual_tree = t8_cmesh_get_face_neighbor (cmesh, 0, 1, &dual_face, &orientation);
  EXPECT_EQ (dual_tree, 1);
  EXPECT_EQ (dual_face, 0);

  EXPECT_EQ (t8_cmesh_get_vertex_conn_status (cmesh), 0);

  t8_cmesh_destroy (&cmesh);
}
