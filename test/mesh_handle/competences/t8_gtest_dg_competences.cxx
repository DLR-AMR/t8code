/*
This file is part of t8code.
t8code is a C library to manage a collection (a forest) of multiple
connected adaptive space-trees of general element classes in parallel.

Copyright (C) 2026 the developers

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
 * \file t8_gtest_dg_competences.cxx
 * Checks that the competences for discontinuous Galerkin methods work as expected.
 */
#include <gtest/gtest.h>
#include <iostream>
#include <t8.h>

#include <mesh_handle/mesh.hxx>
#include <mesh_handle/competences/dg_competences.hxx>
#include <mesh_handle/competence_pack.hxx>
#include <mesh_handle/constructor_wrappers.hxx>

/** Store the rank on each element. */
struct data_per_element
{
  int rank;  ///< Rank of the element.
};
/** Callback function for the mesh handle to decide for refining or coarsening of (a family of) elements.
 * The adaptation criterion is to refine every element with even id.
 * The function header fits the definition of \ref TMesh::adapt_callback_type_with_userdata.
 * \tparam TMeshClass    The mesh handle class.
 * \param [in] mesh      The mesh that should be adapted.
 * \param [in] elements  One element or a family of elements to consider for adaptation.
 * \return 1 if the first entry in \a elements should be refined,
 *        -1 if the family \a elements shall be coarsened,
 *         0 else.
 */
template <typename TMeshClass>
int
mesh_adapt_callback_test_refine_second ([[maybe_unused]] const TMeshClass& mesh,
                                        std::span<const typename TMeshClass::element_class> elements)
{
  if ((elements[0].get_element_handle_id ()) % 2 == 0) {
    return 1;
  }
  return 0;
}

/** TODO
 */
TEST (t8_gtest_dg_competences, remote_ranks)
{
  const int level = 2;
  using namespace t8_mesh_handle;
  using mesh_class
    = mesh<data_element_competences, union_competence_packs_type<mesh_competence_pack<remote_ranks_mesh_competence>,
                                                                 data_mesh_competences<data_per_element>>>;
  auto mesh = handle_hypercube_hybrid_uniform_default<mesh_class> (level, sc_MPI_COMM_WORLD, true, true, false);
  mesh->set_adapt (mesh_adapt_callback_test_refine_second<mesh_class>);
  mesh->set_partition ();
  mesh->set_ghost ();
  mesh->commit ();

  if ((mesh->get_dimension () > 1) && (mesh->get_num_local_elements () > 1)) {
    // Ensure that we actually test with ghost elements.
    ASSERT_GT (mesh->get_num_ghosts (), 0);
  }

  int mpirank;
  int mpiret = sc_MPI_Comm_rank (sc_MPI_COMM_WORLD, &mpirank);
  SC_CHECK_MPI (mpiret);

  // Set local rank for all local mesh elements.
  std::vector<data_per_element> element_data (mesh->get_num_local_elements (), { mpirank });
  mesh->set_element_data (std::move (element_data));
  // Get element data and check that the remote ranks functionality works as expected.
  mesh->exchange_ghost_data ();
  mesh->set_rank_vector ();
  for (const auto& elem : *mesh) {
    EXPECT_EQ (elem.get_element_data ().rank, mpirank);
    EXPECT_EQ (LOCAL_RANK, mesh->get_rank (elem.get_element_handle_id ()));
  }
  for (t8_locidx_t ighost = mesh->get_num_local_elements ();
       ighost < mesh->get_num_local_elements () + mesh->get_num_ghosts (); ighost++) {
    EXPECT_EQ ((*mesh)[ighost].get_element_data ().rank, mesh->get_rank ((*mesh)[ighost].get_element_handle_id ()));
  }
}

TEST (t8_gtest_dg_competences, face_vector_mesh_competence)
{
  const int level = 1;
  using namespace t8_mesh_handle;
  using mesh_class = mesh<element_competence_pack<cache_neighbors>, dg_mesh_competences>;
  auto mesh = handle_hypercube_hybrid_uniform_default<mesh_class> (level, sc_MPI_COMM_WORLD, true, true, false);
  mesh->set_adapt (mesh_adapt_callback_test_refine_second<mesh_class>);
  mesh->set_partition ();
  mesh->set_ghost ();
  mesh->commit ();

  if ((mesh->get_dimension () > 1) && (mesh->get_num_local_elements () > 1)) {
    ASSERT_GT (mesh->get_num_ghosts (), 0);
  }
  // Set the face vector using the competence function.
  mesh->set_unique_face_vector ();

  const auto& faces = mesh->get_unique_face_vector ();
  const auto& element_face_vector = mesh->get_element_face_vector ();

  const t8_locidx_t num_local = mesh->get_num_local_elements ();
  const t8_locidx_t num_ghosts = mesh->get_num_ghosts ();
  ASSERT_EQ (static_cast<t8_locidx_t> (element_face_vector.size ()), num_local + num_ghosts);

  // Check element_face_vector first.
  for (const auto& elem : *mesh) {
    const int handle_id = elem.get_element_handle_id ();
    EXPECT_EQ (static_cast<int> (element_face_vector[handle_id].size ()), elem.get_num_faces ())
      << "Element " << handle_id << " has wrong number of face indices.";
    // Check that the element_face_vector points to valid face indices in the faces vector.
    for (int iface = 0; iface < elem.get_num_faces (); ++iface) {
      EXPECT_GE (element_face_vector[handle_id][iface], 0);
      EXPECT_LE (element_face_vector[handle_id][iface], static_cast<int> (faces.size ()));
    }
  }
  for (t8_locidx_t ighost = num_local; ighost < num_local + num_ghosts; ighost++) {
    const int element_face_vec_size = element_face_vector[ighost].size ();
    EXPECT_GE (element_face_vec_size, 1) << "Ghost element " << ighost << " must have at least one face index.";
    for (int iface = 0; iface < element_face_vec_size; ++iface) {
      EXPECT_GE (element_face_vector[ighost][iface], 0);
      EXPECT_LE (element_face_vector[ighost][iface], static_cast<int> (faces.size ()));
    }
  }

  // Check faces vector.
  for (int iface = 0; iface < static_cast<int> (faces.size ()); ++iface) {
    const auto& face = faces[iface];
    // Check that the face has at least one side and get the element and the neighbors of the first side.
    EXPECT_GE (face.sides.size (), 1) << "Face " << iface << " must have at least one side.";
    const auto elem_first = (*mesh)[face.sides[0].element_id];
    if (elem_first.is_ghost_element ()) {
      continue;
    }
    std::vector<int> dual_faces;
    auto neighs = elem_first.get_face_neighbors (face.sides[0].local_face_id, dual_faces);
    // Check for the case that this face is a boundary.
    if (face.type == face_type::BOUNDARY) {
      EXPECT_EQ (face.sides.size (), 1) << "BOUNDARY face must have exactly 1 side.";
      EXPECT_EQ (face.sides[0].rank, LOCAL_RANK) << "BOUNDARY side must be local.";
      EXPECT_FALSE (elem_first.is_ghost_element ()) << "BOUNDARY side must be local.";
      EXPECT_EQ (neighs.size (), 0) << "BOUNDARY side must not have neighbors.";
    }
    else if (face.type == face_type::CONFORMAL || face.type == face_type::MPI_CONFORMAL) {
      EXPECT_EQ (face.sides.size (), 2) << "CONFORMAL face must have exactly 2 sides.";
      EXPECT_EQ (face.sides[0].rank, LOCAL_RANK) << "CONFORMAL primary side must be local.";
      if (face.type == face_type::CONFORMAL) {
        EXPECT_EQ (face.sides[1].rank, LOCAL_RANK) << "CONFORMAL secondary side must be local.";
      }
      else {
        EXPECT_EQ (face.sides[1].rank, mesh->get_rank (face.sides[1].element_id))
          << "MPI_CONFORMAL secondary side must be remote.";
      }
      EXPECT_EQ (neighs.size (), 1) << "CONFORMAL side must have exactly one neighbor.";
      EXPECT_EQ (neighs[0]->get_element_handle_id (), face.sides[1].element_id)
        << "CONFORMAL side neighbor must match the second side of the face.";
      EXPECT_EQ (dual_faces[0], face.sides[1].local_face_id) << "CONFORMAL side neighbor face index must match.";
      EXPECT_EQ (face.sides[0].orientation, face.sides[1].orientation)
        << "ATTENTION if this is correct we can pull orientation information to the face instead of face sides.";
    }
    else if (face.type == face_type::MORTAR || face.type == face_type::MPI_MORTAR) {
      EXPECT_GT (face.sides.size (), 2) << "MORTAR face must have more than 2 sides.";
      // Check that the first side is the large side (num_neighs>1).
      EXPECT_EQ (neighs.size (), face.sides.size () - 1)
        << "MORTAR side must have as many neighbors as small sides of the face.";
      bool has_remote = std::any_of (face.sides.begin (), face.sides.end (),
                                     [] (const face_side& s) { return s.rank != LOCAL_RANK; });
      EXPECT_EQ (has_remote, face.type == face_type::MPI_MORTAR);
      for (int ineigh = 0; ineigh < static_cast<int> (neighs.size ()); ++ineigh) {
        // The neighbors must not have the same order as the face sides.
        auto neigh_face_side
          = std::find_if (face.sides.begin (), face.sides.end (), [&neighs, ineigh] (const face_side& s) {
              return (neighs[ineigh]->get_element_handle_id () == s.element_id);
            });
        EXPECT_FALSE (neigh_face_side == face.sides.end ()) << "MORTAR side neighbor must be found in the face sides.";
        if (!neighs[ineigh]->is_ghost_element ()) {
          EXPECT_EQ (neighs[ineigh]->get_face_neighbors (dual_faces[ineigh]).size (), 1)
            << "MORTAR side neighbor must have exactly one neighbor on this face.";
        }
        EXPECT_EQ (dual_faces[ineigh], (*neigh_face_side).local_face_id)
          << "MORTAR side neighbor face index must match.";
        EXPECT_EQ (face.sides[0].orientation, (*neigh_face_side).orientation)
          << "ATTENTION if this is correct we can pull orientation information to the face instead of face sides.";
      }
    }

    // This ensures that the element_face_vector is filled correctly and together with the check that
    for (const auto& side : face.sides) {
      bool has_face_in_element_face_vec
        = std::any_of (element_face_vector[side.element_id].begin (), element_face_vector[side.element_id].end (),
                       [iface] (const int faceidx) { return iface == faceidx; });
      EXPECT_TRUE (has_face_in_element_face_vec)
        << "Face index in element_face_vector does not match face index in faces vector.";
    }
  }
}
