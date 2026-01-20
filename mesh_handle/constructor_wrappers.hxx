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

/** \file constructor_wrappers.hxx
 * Wrapper to construct a mesh handle instance from a cmesh. 
 * Additionally, a very small fraction of example coarse meshes is wrapped directly for easy access and 
 * for the use in tests.
 * See \ref t8_cmesh_examples.h for more exemplary cmeshes to be put into the wrapper constructors.
 */

#pragma once

#include <t8.h>
#include <t8_forest/t8_forest_general.h>
#include <t8_cmesh/t8_cmesh.h>
#include <t8_schemes/t8_default/t8_default.hxx>
#include <t8_cmesh/t8_cmesh_examples.h>
#include <memory>

namespace t8_mesh_handle
{

/** Build a uniformly refined mesh handle on a coarse mesh using a scheme.
 * \param [in] cmesh         A coarse mesh.
 * \param [in] scheme        Refinement scheme to use.
 * \param [in] level         An initial uniform refinement level.
 * \param [in] comm          MPI communicator to use.
 * \param [in] do_face_ghost If true, a layer of ghost elements is created.
 * \tparam TMesh             The mesh handle class.
 * \return Unique pointer to a uniformly refined mesh handle with coarse mesh \a cmesh and refinement level \a level.
 */
template <typename TMesh>
std::unique_ptr<TMesh>
handle_new_uniform (const t8_cmesh_t cmesh, const t8_scheme *scheme, const int level, const sc_MPI_Comm comm,
                    const bool do_face_ghost = false)
{
  t8_forest_t forest = t8_forest_new_uniform (cmesh, scheme, level, do_face_ghost, sc_MPI_COMM_WORLD);
  return std::make_unique<TMesh> (forest);
}

/** Build a uniformly refined mesh handle on a coarse mesh using a default scheme.
 * \param [in] cmesh         A coarse mesh.
 * \param [in] level         An initial uniform refinement level.
 * \param [in] comm          MPI communicator to use.
 * \param [in] do_face_ghost If true, a layer of ghost elements is created.
 * \tparam TMesh             The mesh handle class.
 * \return Unique pointer to a uniformly refined mesh handle with coarse mesh \a cmesh and refinement level \a level.
 */
template <typename TMesh>
std::unique_ptr<TMesh>
handle_new_uniform_default (const t8_cmesh_t cmesh, const int level, const sc_MPI_Comm comm,
                            const bool do_face_ghost = false)
{
  return handle_new_uniform<TMesh> (cmesh, t8_scheme_new_default (), level, comm, do_face_ghost);
}

// --- A very small fraction of example coarse meshes is wrapped here for easy access and for the use in tests. ---
/** Hybercube with 6 Tets, 6 Prism, 4 Hex. Refined uniformly to given level using the default scheme.
 * \param [in] level         An initial uniform refinement level.
 * \param [in] comm          MPI communicator to use.
 * \param [in] do_partition  If non-zero create a partitioned cmesh.
 * \param [in] do_face_ghost If true, a layer of ghost elements is created.
 * \param [in] periodic      If non-zero create a periodic cmesh in each direction.
 * \tparam TMesh             The mesh handle class.
 * \return Unique pointer to a uniformly refined mesh handle initially consisting of 6 Tets, 6 prism and 4 hex.
 *         Together, they form a cube.
*/
template <typename TMesh>
std::unique_ptr<TMesh>
handle_hypercube_hybrid_uniform_default (const int level, const sc_MPI_Comm comm, const bool do_partition = false,
                                         const bool do_face_ghost = false, const bool periodic = false)
{
  t8_cmesh_t cmesh = t8_cmesh_new_hypercube_hybrid (comm, do_partition, periodic);
  return handle_new_uniform_default<TMesh> (cmesh, level, comm, do_face_ghost);
}

/** Construct hybercube from one primitive tree class. Refined uniformly to given level using the default scheme.
 * \param [in] eclass        This element class determines the dimension and the number of trees needed to construct a cube.
 * \param [in] level         An initial uniform refinement level.
 * \param [in] comm          MPI communicator to use.
 * \param [in] do_partition  If non-zero create a partitioned cmesh.
 * \param [in] do_face_ghost If true, a layer of ghost elements is created.
 * \param [in] periodic      If non-zero create a periodic cmesh in each direction. Not possible with \a eclass pyramid.
 * \tparam TMesh             The mesh handle class.
 * \return Unique pointer to a uniformly refined mesh handle hypercube.
*/
template <typename TMesh>
std::unique_ptr<TMesh>
handle_hypercube_uniform_default (t8_eclass_t eclass, const int level, const sc_MPI_Comm comm,
                                  const bool do_partition = false, const bool do_face_ghost = false,
                                  const bool periodic = false)
{
  // Broadcast option is hidden from the user.
  t8_cmesh_t cmesh = t8_cmesh_new_hypercube (eclass, comm, 0, do_partition, periodic);
  return handle_new_uniform_default<TMesh> (cmesh, level, comm, do_face_ghost);
}

}  // namespace t8_mesh_handle
