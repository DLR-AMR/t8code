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

/** \file t8_gtest_cmesh_boundary_connectivity.cxx
 * Tests for \ref t8_cmesh_get_boundary_neighbors: on meshes with linear geometry, the neighbors of a tree boundary
 * are exactly the other tree boundaries with the same vertex coordinates, and the corner permutation maps corners
 * onto corners with the same coordinates. Periodic meshes are checked by counting. We also check that the root
 * boundaries of the standalone scheme run along the t8code vertex numbering of edges and faces, so that the corner
 * permutation can be applied to boundary elements directly.
 */

#include <gtest/gtest.h>
#include <t8_eclass.h>
#include <t8_element.h>
#include <t8_cmesh/t8_cmesh.h>
#include <t8_cmesh/t8_cmesh.hxx>
#include <t8_cmesh/t8_cmesh_examples.h>
#include <t8_schemes/t8_standalone/t8_standalone.hxx>
#include <array>
#include <cmath>
#include <functional>
#include <set>
#include <tuple>

/* The tree vertex of corner i of a boundary of the given dimension, in t8code numbering. */
static int
t8_test_boundary_vertex (const t8_eclass_t eclass, const int bdy_dim, const int bdy_id, const int corner)
{
  const int dim = t8_eclass_to_dimension[eclass];
  if (bdy_dim == 0) {
    return bdy_id;
  }
  if (bdy_dim == dim - 1) {
    return t8_face_vertex_to_tree_vertex[eclass][bdy_id][corner];
  }
  return t8_edge_vertex_to_tree_vertex[eclass][bdy_id][corner];
}

static int
t8_test_num_boundaries (const t8_eclass_t eclass, const int bdy_dim)
{
  const int dim = t8_eclass_to_dimension[eclass];
  if (bdy_dim == 0) {
    return t8_eclass_num_vertices[eclass];
  }
  if (bdy_dim == dim - 1) {
    return t8_eclass_num_faces[eclass];
  }
  return t8_eclass_num_edges[eclass];
}

static int
t8_test_num_boundary_corners (const t8_eclass_t eclass, const int bdy_dim, const int bdy_id)
{
  if (bdy_dim == 0) {
    return 1;
  }
  if (bdy_dim == t8_eclass_to_dimension[eclass] - 1) {
    return t8_eclass_num_vertices[t8_eclass_face_types[eclass][bdy_id]];
  }
  return 2;
}

using t8_test_point = std::array<double, 3>;

/* The coordinates of the corners of a tree boundary, rounded to avoid floating point differences. */
static std::vector<t8_test_point>
t8_test_boundary_coords (const t8_cmesh_t cmesh, const t8_locidx_t ltreeid, const int bdy_dim, const int bdy_id)
{
  const t8_eclass_t eclass = t8_cmesh_get_tree_class (cmesh, ltreeid);
  const double *vertices = t8_cmesh_get_tree_vertices (cmesh, ltreeid);
  std::vector<t8_test_point> coords;
  for (int icorner = 0; icorner < t8_test_num_boundary_corners (eclass, bdy_dim, bdy_id); icorner++) {
    const int vertex = t8_test_boundary_vertex (eclass, bdy_dim, bdy_id, icorner);
    t8_test_point point;
    for (int idim = 0; idim < 3; idim++) {
      point[idim] = std::round (vertices[3 * vertex + idim] * 1e8) / 1e8;
    }
    coords.push_back (point);
  }
  return coords;
}

struct cmesh_boundary_connectivity: public testing::TestWithParam<std::pair<const char *, std::function<t8_cmesh_t ()>>>
{
 protected:
  void
  SetUp () override
  {
    cmesh = GetParam ().second ();
  }
  void
  TearDown () override
  {
    t8_cmesh_destroy (&cmesh);
  }
  t8_cmesh_t cmesh;
};

/* The neighbors are exactly the other tree boundaries with the same corner coordinates,
 * and the permutation maps each corner onto the corner with the same coordinates. */
TEST_P (cmesh_boundary_connectivity, matches_coordinates)
{
  const t8_locidx_t num_trees = t8_cmesh_get_num_local_trees (cmesh);
  for (t8_locidx_t itree = 0; itree < num_trees; itree++) {
    const t8_eclass_t eclass = t8_cmesh_get_tree_class (cmesh, itree);
    const int dim = t8_eclass_to_dimension[eclass];
    for (int bdy_dim = 0; bdy_dim < dim; bdy_dim++) {
      for (int bdy_id = 0; bdy_id < t8_test_num_boundaries (eclass, bdy_dim); bdy_id++) {
        const auto coords = t8_test_boundary_coords (cmesh, itree, bdy_dim, bdy_id);
        const std::set<t8_test_point> coord_set (coords.begin (), coords.end ());

        /* All other tree boundaries with the same corners */
        std::set<std::pair<t8_locidx_t, int>> expected;
        for (t8_locidx_t jtree = 0; jtree < num_trees; jtree++) {
          const t8_eclass_t jeclass = t8_cmesh_get_tree_class (cmesh, jtree);
          for (int jbdy = 0; jbdy < t8_test_num_boundaries (jeclass, bdy_dim); jbdy++) {
            const auto jcoords = t8_test_boundary_coords (cmesh, jtree, bdy_dim, jbdy);
            if ((jtree != itree || jbdy != bdy_id) && std::set<t8_test_point> (jcoords.begin (), jcoords.end ()) == coord_set) {
              expected.insert ({ jtree, jbdy });
            }
          }
        }

        std::set<std::pair<t8_locidx_t, int>> found;
        for (const auto &neighbor : t8_cmesh_get_boundary_neighbors (cmesh, itree, bdy_dim, bdy_id)) {
          found.insert ({ neighbor.gtreeid, neighbor.bdy_id });
          const auto neigh_coords = t8_test_boundary_coords (cmesh, neighbor.gtreeid, bdy_dim, neighbor.bdy_id);
          for (size_t icorner = 0; icorner < coords.size (); icorner++) {
            EXPECT_EQ (coords[icorner], neigh_coords[neighbor.perm[icorner]])
              << "tree " << itree << " bdy_dim " << bdy_dim << " bdy_id " << bdy_id << " corner " << icorner
              << " neighbor tree " << neighbor.gtreeid << " bdy_id " << neighbor.bdy_id;
          }
        }
        EXPECT_EQ (found, expected) << "tree " << itree << " bdy_dim " << bdy_dim << " bdy_id " << bdy_id;
      }
    }
  }
}

static const std::pair<const char *, std::function<t8_cmesh_t ()>> t8_test_cmeshes[] = {
  { "line", [] () { return t8_cmesh_new_hypercube_pad (T8_ECLASS_LINE, sc_MPI_COMM_WORLD, std::array<double, 6> { 0, 0, 0, 1, 0, 0 }.data (), 3, 0, 0, 0); } },
  { "quad", [] () { return t8_cmesh_new_hypercube_pad (T8_ECLASS_QUAD, sc_MPI_COMM_WORLD, std::array<double, 12> { 0, 0, 0, 1, 0, 0, 0, 1, 0, 1, 1, 0 }.data (), 3, 2, 0, 0); } },
  { "triangle", [] () { return t8_cmesh_new_hypercube_pad (T8_ECLASS_TRIANGLE, sc_MPI_COMM_WORLD, std::array<double, 12> { 0, 0, 0, 1, 0, 0, 0, 1, 0, 1, 1, 0 }.data (), 3, 2, 0, 0); } },
  { "hex", [] () { return t8_cmesh_new_hypercube_pad (T8_ECLASS_HEX, sc_MPI_COMM_WORLD, std::array<double, 24> { 0, 0, 0, 1, 0, 0, 0, 1, 0, 1, 1, 0, 0, 0, 1, 1, 0, 1, 0, 1, 1, 1, 1, 1 }.data (), 2, 2, 2, 0); } },
  { "tet", [] () { return t8_cmesh_new_hypercube_pad (T8_ECLASS_TET, sc_MPI_COMM_WORLD, std::array<double, 24> { 0, 0, 0, 1, 0, 0, 0, 1, 0, 1, 1, 0, 0, 0, 1, 1, 0, 1, 0, 1, 1, 1, 1, 1 }.data (), 2, 2, 2, 0); } },
  { "prism", [] () { return t8_cmesh_new_hypercube_pad (T8_ECLASS_PRISM, sc_MPI_COMM_WORLD, std::array<double, 24> { 0, 0, 0, 1, 0, 0, 0, 1, 0, 1, 1, 0, 0, 0, 1, 1, 0, 1, 0, 1, 1, 1, 1, 1 }.data (), 2, 2, 2, 0); } },
  { "pyramid", [] () { return t8_cmesh_new_hypercube (T8_ECLASS_PYRAMID, sc_MPI_COMM_WORLD, 0, 0, 0); } },
  { "hybrid", [] () { return t8_cmesh_new_hypercube_hybrid (sc_MPI_COMM_WORLD, 0, 0); } },
  { "full_hybrid", [] () { return t8_cmesh_new_full_hybrid (sc_MPI_COMM_WORLD); } },
  { "prism_cake", [] () { return t8_cmesh_new_prism_cake (sc_MPI_COMM_WORLD, 5); } },
  { "pyramid_cake", [] () { return t8_cmesh_new_pyramid_cake (sc_MPI_COMM_WORLD, 5); } },
};

INSTANTIATE_TEST_SUITE_P (t8_gtest_cmesh_boundary_connectivity, cmesh_boundary_connectivity,
                          testing::ValuesIn (t8_test_cmeshes),
                          [] (const auto &info) { return std::string (info.param.first); });

/* Each boundary is a neighbor of its neighbors, with the inverse permutation. */
struct cmesh_boundary_connectivity_symmetry: public cmesh_boundary_connectivity
{
};

TEST_P (cmesh_boundary_connectivity_symmetry, is_symmetric)
{
  for (t8_locidx_t itree = 0; itree < t8_cmesh_get_num_local_trees (cmesh); itree++) {
    const t8_eclass_t eclass = t8_cmesh_get_tree_class (cmesh, itree);
    const int dim = t8_eclass_to_dimension[eclass];
    for (int bdy_dim = 0; bdy_dim < dim; bdy_dim++) {
      for (int bdy_id = 0; bdy_id < t8_test_num_boundaries (eclass, bdy_dim); bdy_id++) {
        const int num_corners = t8_test_num_boundary_corners (eclass, bdy_dim, bdy_id);
        for (const auto &neighbor : t8_cmesh_get_boundary_neighbors (cmesh, itree, bdy_dim, bdy_id)) {
          int num_back = 0;
          for (const auto &back : t8_cmesh_get_boundary_neighbors (cmesh, neighbor.gtreeid, bdy_dim, neighbor.bdy_id)) {
            if (back.gtreeid == itree && back.bdy_id == bdy_id) {
              num_back++;
              for (int icorner = 0; icorner < num_corners; icorner++) {
                EXPECT_EQ (back.perm[neighbor.perm[icorner]], icorner)
                  << "tree " << itree << " bdy_dim " << bdy_dim << " bdy_id " << bdy_id;
              }
            }
          }
          EXPECT_EQ (num_back, 1) << "tree " << itree << " bdy_dim " << bdy_dim << " bdy_id " << bdy_id;
        }
      }
    }
  }
}

/* Meshes whose face connections do not match their coordinates are only checked for symmetry. */
static const std::pair<const char *, std::function<t8_cmesh_t ()>> t8_test_orientation_cmeshes[] = {
  { "prism_cake_funny_oriented", [] () { return t8_cmesh_new_prism_cake_funny_oriented (sc_MPI_COMM_WORLD); } },
  { "tet_orientation_test", [] () { return t8_cmesh_new_tet_orientation_test (sc_MPI_COMM_WORLD); } },
  { "periodic_hybrid", [] () { return t8_cmesh_new_periodic_hybrid (sc_MPI_COMM_WORLD); } },
  { "periodic_tri", [] () { return t8_cmesh_new_periodic_tri (sc_MPI_COMM_WORLD); } },
  { "periodic_hex_brick", [] () { return t8_cmesh_new_brick_3d (2, 3, 2, 1, 1, 1, sc_MPI_COMM_WORLD); } },
};

/* In a partitioned cmesh, each local tree has the same neighbors as in the replicated cmesh. */
struct cmesh_boundary_connectivity_partitioned: public cmesh_boundary_connectivity
{
};

TEST_P (cmesh_boundary_connectivity_partitioned, matches_replicated)
{
  /* Partition the trees uniformly, so that with several processes trees are shared */
  t8_cmesh_t partitioned;
  t8_cmesh_init (&partitioned);
  t8_cmesh_ref (cmesh);
  t8_cmesh_set_derive (partitioned, cmesh);
  t8_cmesh_set_partition_uniform (partitioned, 1, t8_scheme_new_standalone ());
  t8_cmesh_commit (partitioned, sc_MPI_COMM_WORLD);
  t8_cmesh_build_boundary_connectivity (partitioned, sc_MPI_COMM_WORLD);

  using t8_test_neighbor = std::tuple<t8_gloidx_t, int, int, std::vector<int>>;
  for (t8_locidx_t itree = 0; itree < t8_cmesh_get_num_local_trees (partitioned); itree++) {
    const t8_gloidx_t gtreeid = t8_cmesh_get_global_id (partitioned, itree);
    const t8_eclass_t eclass = t8_cmesh_get_tree_class (partitioned, itree);
    const int dim = t8_eclass_to_dimension[eclass];
    for (int bdy_dim = 0; bdy_dim < dim; bdy_dim++) {
      for (int bdy_id = 0; bdy_id < t8_test_num_boundaries (eclass, bdy_dim); bdy_id++) {
        const int num_corners = t8_test_num_boundary_corners (eclass, bdy_dim, bdy_id);
        auto neighbor_set = [&] (const std::vector<t8_cmesh_boundary_neighbor> &neighbors) {
          std::set<t8_test_neighbor> result;
          for (const auto &neighbor : neighbors) {
            result.insert ({ neighbor.gtreeid, neighbor.eclass, neighbor.bdy_id,
                             std::vector<int> (neighbor.perm, neighbor.perm + num_corners) });
          }
          return result;
        };
        EXPECT_EQ (neighbor_set (t8_cmesh_get_boundary_neighbors (partitioned, itree, bdy_dim, bdy_id)),
                   neighbor_set (t8_cmesh_get_boundary_neighbors (cmesh, gtreeid, bdy_dim, bdy_id)))
          << "tree " << gtreeid << " bdy_dim " << bdy_dim << " bdy_id " << bdy_id;
      }
    }
  }
  t8_cmesh_destroy (&partitioned);
}

INSTANTIATE_TEST_SUITE_P (t8_gtest_cmesh_boundary_connectivity_geometric, cmesh_boundary_connectivity_partitioned,
                          testing::ValuesIn (t8_test_cmeshes),
                          [] (const auto &info) { return std::string (info.param.first); });
INSTANTIATE_TEST_SUITE_P (t8_gtest_cmesh_boundary_connectivity_orientation, cmesh_boundary_connectivity_partitioned,
                          testing::ValuesIn (t8_test_orientation_cmeshes),
                          [] (const auto &info) { return std::string (info.param.first); });

INSTANTIATE_TEST_SUITE_P (t8_gtest_cmesh_boundary_connectivity_geometric, cmesh_boundary_connectivity_symmetry,
                          testing::ValuesIn (t8_test_cmeshes),
                          [] (const auto &info) { return std::string (info.param.first); });
INSTANTIATE_TEST_SUITE_P (t8_gtest_cmesh_boundary_connectivity_orientation, cmesh_boundary_connectivity_symmetry,
                          testing::ValuesIn (t8_test_orientation_cmeshes),
                          [] (const auto &info) { return std::string (info.param.first); });


/* In the periodic unit cube of one tree, all vertices are identified, edges with the parallel edges in the same
 * direction and faces with the opposite face. */
TEST (cmesh_boundary_connectivity_periodic, unit_cube)
{
  for (int dim = 1; dim <= 3; dim++) {
    t8_cmesh_t cmesh = t8_cmesh_new_periodic (sc_MPI_COMM_WORLD, dim);
    const t8_eclass_t eclass = t8_cmesh_get_tree_class (cmesh, 0);
    for (int bdy_dim = 0; bdy_dim < dim; bdy_dim++) {
      for (int bdy_id = 0; bdy_id < t8_test_num_boundaries (eclass, bdy_dim); bdy_id++) {
        const auto neighbors = t8_cmesh_get_boundary_neighbors (cmesh, 0, bdy_dim, bdy_id);
        /* The number of copies of a k-dimensional boundary in the cube is 2^(dim - k) */
        EXPECT_EQ ((int) neighbors.size (), (1 << (dim - bdy_dim)) - 1);
        for (const auto &neighbor : neighbors) {
          EXPECT_EQ (neighbor.gtreeid, 0);
          if (bdy_dim > 0) {
            /* Parallel boundaries have the same direction and are mapped by a translation. */
            for (int icorner = 0; icorner < t8_test_num_boundary_corners (eclass, bdy_dim, bdy_id); icorner++) {
              EXPECT_EQ (neighbor.perm[icorner], icorner);
            }
            if (bdy_dim == 1 && dim == 3) {
              EXPECT_EQ (neighbor.bdy_id / 4, bdy_id / 4);
            }
            else {
              EXPECT_EQ (neighbor.bdy_id / 2, bdy_id / 2);
            }
          }
        }
      }
    }
    t8_cmesh_destroy (&cmesh);
  }
}

/* The first and last descendant of the standalone root element at a tree boundary contain the first and last
 * corner of the boundary in t8code numbering, so the root boundaries run along the t8code corner order. */
TEST (cmesh_boundary_connectivity_standalone, root_boundary_corners)
{
  const t8_scheme *scheme = t8_scheme_new_standalone ();
  for (const t8_eclass_t eclass : { T8_ECLASS_LINE, T8_ECLASS_QUAD, T8_ECLASS_TRIANGLE, T8_ECLASS_HEX, T8_ECLASS_TET,
                                    T8_ECLASS_PRISM, T8_ECLASS_PYRAMID }) {
    const int dim = t8_eclass_to_dimension[eclass];
    t8_element_t *root, *desc;
    scheme->element_new (eclass, 1, &root);
    scheme->element_new (eclass, 1, &desc);
    scheme->set_to_root (eclass, root);
    const int level = 3;
    for (int bdy_dim = 0; bdy_dim < dim; bdy_dim++) {
      for (int bdy_id = 0; bdy_id < t8_test_num_boundaries (eclass, bdy_dim); bdy_id++) {
        const int num_corners = t8_test_num_boundary_corners (eclass, bdy_dim, bdy_id);
        for (const int last : { 0, 1 }) {
          if (last) {
            scheme->element_get_last_descendant_boundary (eclass, root, bdy_dim, bdy_id, desc, level);
          }
          else {
            scheme->element_get_first_descendant_boundary (eclass, root, bdy_dim, bdy_id, desc, level);
          }
          const int vertex = t8_test_boundary_vertex (eclass, bdy_dim, bdy_id, last ? num_corners - 1 : 0);
          bool contains = false;
          for (int icorner = 0; icorner < scheme->element_get_num_corners (eclass, desc); icorner++) {
            double coords[3] = { 0, 0, 0 };
            scheme->element_get_vertex_reference_coords (eclass, desc, icorner, coords);
            bool equal = true;
            for (int idim = 0; idim < dim; idim++) {
              equal &= coords[idim] == t8_element_corner_ref_coords[eclass][vertex][idim];
            }
            contains |= equal;
          }
          EXPECT_TRUE (contains) << t8_eclass_to_string[eclass] << " bdy_dim " << bdy_dim << " bdy_id " << bdy_id
                                 << (last ? " last" : " first");
        }
      }
    }
    scheme->element_destroy (eclass, 1, &root);
    scheme->element_destroy (eclass, 1, &desc);
  }
  scheme->unref ();
}
