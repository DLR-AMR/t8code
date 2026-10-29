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

/** \file t8_gtest_boundaries.cxx
 * Tests for the boundaries of all dimensions of the standalone scheme: boundary counts, children,
 * first and last descendants and neighbors at boundaries, and their agreement with the face functions.
 */

#include <gtest/gtest.h>
#include <t8_eclass.h>
#include <test/t8_gtest_schemes.hxx>
#include <test/t8_gtest_custom_assertion.hxx>
#include <test/t8_gtest_macros.hxx>
#include "t8_gtest_dfs_base.hxx"

/* Compute the first or last descendant of an element at a boundary by descending through the first or last
 * child at the boundary, see also t8_test_manual_first_last_face_descendant. */
static void
t8_test_manual_boundary_descendant (const t8_scheme *scheme, const t8_element_t *element, const t8_eclass_t tree_class,
                                    const int bdy_dim, const int bdy_id, const int desc_level, const int last,
                                    t8_element_t *desc)
{
  const int max_children = scheme->get_max_num_boundary_children (tree_class);
  t8_element_t **children = T8_TESTSUITE_ALLOC (t8_element_t *, max_children);
  int *child_bdy_ids = T8_TESTSUITE_ALLOC (int, max_children);
  scheme->element_new (tree_class, max_children, children);
  scheme->element_copy (tree_class, element, desc);
  int current_bdy_id = bdy_id;
  for (int ilevel = scheme->element_get_level (tree_class, element); ilevel < desc_level; ilevel++) {
    const int num_children = scheme->element_get_num_boundary_children (tree_class, desc, bdy_dim, current_bdy_id);
    scheme->element_get_children_at_boundary (tree_class, desc, bdy_dim, current_bdy_id, children, num_children,
                                              child_bdy_ids);
    const int ichild = last ? num_children - 1 : 0;
    scheme->element_copy (tree_class, children[ichild], desc);
    current_bdy_id = child_bdy_ids[ichild];
  }
  scheme->element_destroy (tree_class, max_children, children);
  T8_TESTSUITE_FREE (children);
  T8_TESTSUITE_FREE (child_bdy_ids);
}

struct class_boundaries: public TestDFS
{
 private:
  void
  check_element () override
  {
    const int dim = t8_eclass_to_dimension[eclass];
    const int level = scheme->element_get_level (eclass, element);

    /* The Euler characteristic of the boundary of a dim-dimensional polytope is 1 - (-1)^dim */
    int euler = 0;
    for (int bdy_dim = 0; bdy_dim < dim; bdy_dim++) {
      euler += (bdy_dim % 2 ? -1 : 1) * scheme->element_get_num_boundaries (eclass, element, bdy_dim);
    }
    EXPECT_EQ (euler, 1 - (dim % 2 ? -1 : 1));

    /* The boundaries of dimension dim - 1 are the faces */
    const int num_faces = scheme->element_get_num_faces (eclass, element);
    EXPECT_EQ (scheme->element_get_num_boundaries (eclass, element, dim - 1), num_faces);
    for (int iface = 0; iface < num_faces; iface++) {
      const int num_face_children = scheme->element_get_num_face_children (eclass, element, iface);
      ASSERT_EQ (scheme->element_get_num_boundary_children (eclass, element, dim - 1, iface), num_face_children);
      int face_child_ids[T8_TEST_MAX_CHILDREN], child_bdy_ids[T8_TEST_MAX_CHILDREN];
      scheme->element_get_children_at_face (eclass, element, iface, face_children, num_face_children, face_child_ids);
      scheme->element_get_children_at_boundary (eclass, element, dim - 1, iface, bdy_children, num_face_children,
                                                child_bdy_ids);
      for (int ichild = 0; ichild < num_face_children; ichild++) {
        EXPECT_ELEM_EQ (scheme, eclass, face_children[ichild], bdy_children[ichild]);
        EXPECT_EQ (child_bdy_ids[ichild], scheme->element_face_get_child_face (eclass, element, iface, ichild));
      }
      if (level < max_desc_level) {
        scheme->element_get_first_descendant_face (eclass, element, iface, desc_a, max_desc_level);
        scheme->element_get_first_descendant_boundary (eclass, element, dim - 1, iface, desc_b, max_desc_level);
        EXPECT_ELEM_EQ (scheme, eclass, desc_a, desc_b);
        scheme->element_get_last_descendant_face (eclass, element, iface, desc_a, max_desc_level);
        scheme->element_get_last_descendant_boundary (eclass, element, dim - 1, iface, desc_b, max_desc_level);
        EXPECT_ELEM_EQ (scheme, eclass, desc_a, desc_b);
      }
    }

    for (int bdy_dim = 0; bdy_dim < dim; bdy_dim++) {
      const int num_boundaries = scheme->element_get_num_boundaries (eclass, element, bdy_dim);
      for (int bdy_id = 0; bdy_id < num_boundaries; bdy_id++) {
        check_children (bdy_dim, bdy_id);
        check_descendants (bdy_dim, bdy_id, level);
        check_neighbors (bdy_dim, bdy_id);
      }
    }
  }

  /* The children at a boundary are children of the element in SFC order, and a child boundary is not more than
   * the number of boundaries of the child. */
  void
  check_children (const int bdy_dim, const int bdy_id)
  {
    if (scheme->element_get_level (eclass, element) >= scheme->get_maxlevel (eclass)) {
      return;
    }
    const int num_children = scheme->element_get_num_boundary_children (eclass, element, bdy_dim, bdy_id);
    ASSERT_GT (num_children, 0);
    ASSERT_LE (num_children, scheme->get_max_num_boundary_children (eclass));
    int child_bdy_ids[T8_TEST_MAX_CHILDREN];
    scheme->element_get_children_at_boundary (eclass, element, bdy_dim, bdy_id, bdy_children, num_children,
                                              child_bdy_ids);
    int last_child_id = -1;
    for (int ichild = 0; ichild < num_children; ichild++) {
      scheme->element_get_parent (eclass, bdy_children[ichild], desc_a);
      EXPECT_ELEM_EQ (scheme, eclass, desc_a, element);
      const int child_id = scheme->element_get_child_id (eclass, bdy_children[ichild]);
      EXPECT_GT (child_id, last_child_id);
      last_child_id = child_id;
      EXPECT_GE (child_bdy_ids[ichild], 0);
      EXPECT_LT (child_bdy_ids[ichild], scheme->element_get_num_boundaries (eclass, bdy_children[ichild], bdy_dim));
    }

    /* The children touching the boundary in any dimension are children of the element in SFC order, and the ones
     * touching it with full dimension are the children at the boundary. */
    const int num_touching = scheme->element_get_num_boundary_touching_children (eclass, element, bdy_dim, bdy_id);
    ASSERT_GE (num_touching, num_children);
    ASSERT_LE (num_touching, scheme->get_max_num_boundary_touching_children (eclass));
    int touching_dims[T8_TEST_MAX_CHILDREN], touching_ids[T8_TEST_MAX_CHILDREN];
    scheme->element_get_children_touching_boundary (eclass, element, bdy_dim, bdy_id, touching_children, num_touching,
                                                    touching_dims, touching_ids);
    int ifull = 0;
    last_child_id = -1;
    for (int ichild = 0; ichild < num_touching; ichild++) {
      scheme->element_get_parent (eclass, touching_children[ichild], desc_a);
      EXPECT_ELEM_EQ (scheme, eclass, desc_a, element);
      const int child_id = scheme->element_get_child_id (eclass, touching_children[ichild]);
      EXPECT_GT (child_id, last_child_id);
      last_child_id = child_id;
      ASSERT_GE (touching_dims[ichild], 0);
      ASSERT_LE (touching_dims[ichild], bdy_dim);
      EXPECT_LT (touching_ids[ichild],
                 scheme->element_get_num_boundaries (eclass, touching_children[ichild], touching_dims[ichild]));
      if (touching_dims[ichild] == bdy_dim) {
        ASSERT_LT (ifull, num_children);
        EXPECT_ELEM_EQ (scheme, eclass, touching_children[ichild], bdy_children[ifull]);
        EXPECT_EQ (touching_ids[ichild], child_bdy_ids[ifull]);
        ifull++;
      }
    }
    EXPECT_EQ (ifull, num_children);
  }

  /* The first and last descendant at a boundary are reached by descending through the first and last child at
   * the boundary, which also checks that they keep the type of the element. */
  void
  check_descendants (const int bdy_dim, const int bdy_id, const int level)
  {
    for (int desc_level = level + 1; desc_level <= max_desc_level; desc_level++) {
      for (int last = 0; last < 2; last++) {
        t8_test_manual_boundary_descendant (scheme, element, eclass, bdy_dim, bdy_id, desc_level, last, desc_a);
        if (last) {
          scheme->element_get_last_descendant_boundary (eclass, element, bdy_dim, bdy_id, desc_b, desc_level);
        }
        else {
          scheme->element_get_first_descendant_boundary (eclass, element, bdy_dim, bdy_id, desc_b, desc_level);
        }
        EXPECT_ELEM_EQ (scheme, eclass, desc_a, desc_b);
      }
      /* The first descendant does not come after the last one in SFC order */
      scheme->element_get_first_descendant_boundary (eclass, element, bdy_dim, bdy_id, desc_a, desc_level);
      scheme->element_get_last_descendant_boundary (eclass, element, bdy_dim, bdy_id, desc_b, desc_level);
      EXPECT_LE (scheme->element_get_linear_id (eclass, desc_a, desc_level),
                 scheme->element_get_linear_id (eclass, desc_b, desc_level));
    }
  }

  /* The neighbors at a boundary have the element as neighbor at their boundary, and at vertices they share
   * the vertex. */
  void
  check_neighbors (const int bdy_dim, const int bdy_id)
  {
    const int max_neighbors = scheme->get_max_num_boundary_neighbors (eclass);
    int neigh_bdy_ids[T8_TEST_MAX_NEIGHBORS], back_bdy_ids[T8_TEST_MAX_NEIGHBORS];
    ASSERT_LE (max_neighbors, T8_TEST_MAX_NEIGHBORS);
    const int num_neighbors
      = scheme->element_get_boundary_neighbors_inside (eclass, element, bdy_dim, bdy_id, neighbors, neigh_bdy_ids);
    ASSERT_LT (num_neighbors, max_neighbors);
    for (int ineigh = 0; ineigh < num_neighbors; ineigh++) {
      EXPECT_FALSE (scheme->element_is_equal (eclass, neighbors[ineigh], element));
      EXPECT_EQ (scheme->element_get_level (eclass, neighbors[ineigh]), scheme->element_get_level (eclass, element));
      const int num_back = scheme->element_get_boundary_neighbors_inside (
        eclass, neighbors[ineigh], bdy_dim, neigh_bdy_ids[ineigh], back_neighbors, back_bdy_ids);
      int found = 0;
      for (int iback = 0; iback < num_back; iback++) {
        found += scheme->element_is_equal (eclass, back_neighbors[iback], element) && back_bdy_ids[iback] == bdy_id;
      }
      EXPECT_EQ (found, 1) << "The element is not a neighbor of its neighbor at boundary " << bdy_id
                           << " of dimension " << bdy_dim;
      if (bdy_dim == 0) {
        double coords[3], neigh_coords[3];
        scheme->element_get_vertex_reference_coords (eclass, element, bdy_id, coords);
        scheme->element_get_vertex_reference_coords (eclass, neighbors[ineigh], neigh_bdy_ids[ineigh], neigh_coords);
        for (int idim = 0; idim < t8_eclass_to_dimension[eclass]; idim++) {
          EXPECT_NEAR (coords[idim], neigh_coords[idim], 1e-12);
        }
      }
    }
  }

 protected:
  static constexpr int T8_TEST_MAX_CHILDREN = 10;
  static constexpr int T8_TEST_MAX_NEIGHBORS = 32;

  void
  SetUp () override
  {
    dfs_test_setup ();
    scheme->element_new (eclass, 1, &desc_a);
    scheme->element_new (eclass, 1, &desc_b);
    scheme->element_new (eclass, T8_TEST_MAX_CHILDREN, face_children);
    scheme->element_new (eclass, T8_TEST_MAX_CHILDREN, bdy_children);
    scheme->element_new (eclass, T8_TEST_MAX_CHILDREN, touching_children);
    scheme->element_new (eclass, T8_TEST_MAX_NEIGHBORS, neighbors);
    scheme->element_new (eclass, T8_TEST_MAX_NEIGHBORS, back_neighbors);
  }
  void
  TearDown () override
  {
    scheme->element_destroy (eclass, 1, &desc_a);
    scheme->element_destroy (eclass, 1, &desc_b);
    scheme->element_destroy (eclass, T8_TEST_MAX_CHILDREN, face_children);
    scheme->element_destroy (eclass, T8_TEST_MAX_CHILDREN, bdy_children);
    scheme->element_destroy (eclass, T8_TEST_MAX_CHILDREN, touching_children);
    scheme->element_destroy (eclass, T8_TEST_MAX_NEIGHBORS, neighbors);
    scheme->element_destroy (eclass, T8_TEST_MAX_NEIGHBORS, back_neighbors);
    dfs_test_teardown ();
  }

#if T8_TEST_LEVEL_INT >= 1
  const int max_dfs_level = 2;
#else
  const int max_dfs_level = 3;
#endif
  const int max_desc_level = 5;
  t8_element_t *desc_a;
  t8_element_t *desc_b;
  t8_element_t *face_children[T8_TEST_MAX_CHILDREN];
  t8_element_t *bdy_children[T8_TEST_MAX_CHILDREN];
  t8_element_t *touching_children[T8_TEST_MAX_CHILDREN];
  t8_element_t *neighbors[T8_TEST_MAX_NEIGHBORS];
  t8_element_t *back_neighbors[T8_TEST_MAX_NEIGHBORS];

 public:
  void
  run ()
  {
    check_recursive_dfs_to_max_lvl (max_dfs_level);
  }
};

TEST_P (class_boundaries, t8_check_boundaries)
{
  run ();
}

/* Boundaries are only implemented for the standalone scheme (scheme id 1), and vertices have no boundaries */
INSTANTIATE_TEST_SUITE_P (t8_gtest_boundaries, class_boundaries,
                          ::testing::Combine (::testing::Values (1),
                                              ::testing::Range (T8_ECLASS_LINE, T8_ECLASS_COUNT)),
                          print_all_schemes);
