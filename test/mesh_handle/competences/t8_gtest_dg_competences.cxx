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
 * Checks that the competences for discontinuous Galerkin methods defined in \ref dg_competences.hxx work as expected.
 */
#include <gtest/gtest.h>
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
 * The function header fits the definition of \ref TMesh::adapt_callback_type.
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

/** Check the competence remote_ranks_mesh_competence for correctness. The ranks are set as data first, exchanged for 
 * ghost elements and then checked against the competence functionality.
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

  const t8_locidx_t num_local = mesh->get_num_local_elements ();
  const t8_locidx_t num_ghosts = mesh->get_num_ghosts ();
  if ((mesh->get_dimension () > 1) && (num_local > 1)) {
    // Ensure that we actually test with ghost elements.
    ASSERT_GT (num_ghosts, 0);
  }

  int mpirank;
  int mpiret = sc_MPI_Comm_rank (sc_MPI_COMM_WORLD, &mpirank);
  SC_CHECK_MPI (mpiret);

  // Set local rank for all local mesh elements.
  std::vector<data_per_element> element_data (num_local, { mpirank });
  mesh->set_element_data (std::move (element_data));
  // Get element data and check that the remote ranks competence works as expected.
  mesh->exchange_ghost_data ();
  mesh->fill_rank_vector ();
  for (const auto& elem : *mesh) {
    EXPECT_EQ (elem.get_element_data ().rank, mpirank);
    EXPECT_EQ (LOCAL_RANK, mesh->get_rank (elem.get_element_handle_id ()));
  }
  for (t8_locidx_t ighost = num_local; ighost < num_local + num_ghosts; ighost++) {
    EXPECT_EQ ((*mesh)[ighost].get_element_data ().rank, mesh->get_rank ((*mesh)[ighost].get_element_handle_id ()));
  }
}

/** Check the competence face_vector_mesh_competence for correctness. 
 * The test checks that the face vector and the element-face vector are consistent with each other and with the 
 * neighbors of the elements.
 */
TEST (t8_gtest_dg_competences, face_vector_mesh_competence)
{
  const int level = 2;
  using namespace t8_mesh_handle;
  using mesh_class = mesh<element_competence_pack<cache_neighbors>, dg_mesh_competences>;
  auto mesh = handle_hypercube_hybrid_uniform_default<mesh_class> (level, sc_MPI_COMM_WORLD, true, true, false);
  mesh->set_adapt (mesh_adapt_callback_test_refine_second<mesh_class>);
  mesh->set_partition ();
  mesh->set_ghost ();
  mesh->commit ();

  // Set the face vector using the competence function.
  mesh->fill_unique_face_vector ();

  const auto& faces = mesh->get_unique_face_vector ();
  const auto& element_face_vector = mesh->get_element_face_vector ();

  const t8_locidx_t num_local = mesh->get_num_local_elements ();
  const t8_locidx_t num_ghosts = mesh->get_num_ghosts ();
  ASSERT_EQ (static_cast<t8_locidx_t> (element_face_vector.size ()), num_local + num_ghosts);

  // Check element_face_vector first.
  for (t8_locidx_t ielem = num_local; ielem < num_local + num_ghosts; ielem++) {
    EXPECT_EQ (static_cast<int> (element_face_vector[ielem].size ()), (*mesh)[ielem].get_num_faces ());
    // Check that the element_face_vector points to valid face indices in the faces vector.
    for (int iface = 0; iface < (*mesh)[ielem].get_num_faces (); ++iface) {
      if (element_face_vector[ielem][iface] == -1) {
        EXPECT_TRUE ((*mesh)[ielem].is_ghost_element ());
      }
      else {
        EXPECT_GE (element_face_vector[ielem][iface], 0);
        EXPECT_LE (element_face_vector[ielem][iface], static_cast<int> (faces.size ()));
      }
    }
  }

  // Check faces vector.
  for (int iface = 0; iface < static_cast<int> (faces.size ()); ++iface) {
    const auto& face = faces[iface];
    // Check that the face has at least one side and get the element and the neighbors of the first side.
    EXPECT_GE (face.sides.size (), 1);
    const auto elem_first = (*mesh)[face.sides[0].element_id];

    std::vector<int> dual_faces;
    auto neighs = elem_first.get_face_neighbors (face.sides[0].local_face_id, dual_faces);
    /* --- BOUNDARY --- */
    switch (face.type) {
    case face_type::BOUNDARY:
      EXPECT_EQ (face.sides.size (), 1) << "BOUNDARY face must have exactly 1 side.";
      EXPECT_EQ (face.sides[0].rank, LOCAL_RANK) << "BOUNDARY side must be local.";
      EXPECT_FALSE (elem_first.is_ghost_element ()) << "BOUNDARY side element must be local.";
      EXPECT_EQ (neighs.size (), 0) << "BOUNDARY side must not have neighbors.";
      EXPECT_EQ (face.orientation,
                 t8_forest_leaf_face_orientation (mesh->get_forest (), elem_first.get_local_tree_id (),
                                                  t8_forest_get_scheme (mesh->get_forest ()),
                                                  elem_first.get_forest_element (), face.sides[0].local_face_id))
        << "BOUNDARY side face orientation must match the one computed from the forest.";
      break;
    /* --- CONFORMAL --- */
    case face_type::CONFORMAL:
      EXPECT_EQ (face.sides[1].rank, LOCAL_RANK) << "CONFORMAL secondary side must be local.";
      [[fallthrough]];
    case face_type::MPI_CONFORMAL:
      if (face.type == face_type::MPI_CONFORMAL) {
        EXPECT_EQ (face.sides[1].rank, mesh->get_rank (face.sides[1].element_id))
          << "MPI_CONFORMAL secondary side must be remote.";
      }
      EXPECT_EQ (face.sides.size (), 2) << "CONFORMAL face must have exactly 2 sides.";
      EXPECT_EQ (face.sides[0].rank, LOCAL_RANK) << "CONFORMAL primary side must be local.";

      EXPECT_EQ (neighs.size (), 1) << "CONFORMAL side must have exactly one neighbor.";
      EXPECT_EQ (neighs[0]->get_element_handle_id (), face.sides[1].element_id)
        << "CONFORMAL side neighbor must match the second side of the face.";
      EXPECT_EQ (dual_faces[0], face.sides[1].local_face_id) << "CONFORMAL side neighbor face index must match.";
      EXPECT_EQ (face.orientation,
                 t8_forest_leaf_face_orientation (mesh->get_forest (), neighs[0]->get_local_tree_id (),
                                                  t8_forest_get_scheme (mesh->get_forest ()),
                                                  neighs[0]->get_forest_element (), dual_faces[0]))
        << "CONFORMAL side face orientation must match.";
      break;
    /* --- MORTAR --- */
    case face_type::MORTAR:
      // Check that first element is the large mortar.
      // For MPI_MORTARs, the large mortar could have only one neighbor.
      EXPECT_GT (neighs.size (), 1) << "MORTAR face must have more than 2 sides.";
      [[fallthrough]];
    case face_type::MPI_MORTAR: {
      bool has_remote = std::any_of (face.sides.begin (), face.sides.end (),
                                     [] (const face_side& s) { return s.rank != LOCAL_RANK; });
      EXPECT_EQ (has_remote, face.type == face_type::MPI_MORTAR);
      EXPECT_EQ (neighs[0]->get_level (), elem_first.get_level () + 1)
        << "MORTAR first element must be the large mortar.";
      EXPECT_EQ (neighs.size (), face.sides.size () - 1)
        << "MORTAR side must have as many neighbors as small sides of the face.";

      for (int ineigh = 0; ineigh < static_cast<int> (neighs.size ()); ++ineigh) {
        // The neighbors do not necessarily have the same order as the face sides.
        auto neigh_face_side
          = std::find_if (face.sides.begin (), face.sides.end (), [&neighs, ineigh] (const face_side& s) {
              return (neighs[ineigh]->get_element_handle_id () == s.element_id);
            });
        EXPECT_FALSE (neigh_face_side == face.sides.end ()) << "MORTAR side neighbor must be found in the face sides.";
        EXPECT_EQ (neighs[ineigh]->get_face_neighbors (dual_faces[ineigh]).size (), 1)
          << "MORTAR side neighbor must have exactly one neighbor on this face.";
        EXPECT_EQ (dual_faces[ineigh], (*neigh_face_side).local_face_id)
          << "MORTAR side neighbor face index must match.";
        EXPECT_EQ (face.orientation,
                   t8_forest_leaf_face_orientation (mesh->get_forest (), neighs[ineigh]->get_local_tree_id (),
                                                    t8_forest_get_scheme (mesh->get_forest ()),
                                                    neighs[ineigh]->get_forest_element (), dual_faces[ineigh]))
          << "MORTAR side face orientation must match.";
      }
      break;
    }
    default:
      EXPECT_TRUE (false) << "The face type should be valid such that this code is unreachable.";
    }

    // This ensures that the element_face_vector is filled correctly.
    for (const auto& side : face.sides) {
      EXPECT_EQ (element_face_vector[side.element_id][side.local_face_id], iface);
    }
  }
}
