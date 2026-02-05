/*
  This file is part of t8code.
  t8code is a C library to manage a collection (a forest) of multiple
  connected adaptive space-trees of general element classes in parallel.

  Copyright (C) 2024 the developers

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

/** \file t8_cmesh_params.hxx
 *  File to store the parameters that are used for our parameterized cmesh-tests.
 */
#ifndef T8_CMESH_PARAMS_HXX
#define T8_CMESH_PARAMS_HXX

#include <t8_eclass.h>
#include <vector>
#include <numeric>
#include <algorithm>
#include <t8_geometry/t8_geometry_base.hxx>
#include <t8_geometry/t8_geometry_implementations/t8_geometry_linear.hxx>
#include <t8_geometry/t8_geometry_implementations/t8_geometry_linear_axis_aligned.hxx>

#if T8_TEST_LEVEL_INT >= 2
#define T8_CMESH_MAX_NUM_OF_TREES 2 /**< Maximum number of trees in the cmesh. */
// Note: This is not the number of prisms used in the test but the number of different settings tested.
#define T8_CMESH_MAX_NUM_OF_PRISMS 2 /**< Maximum number of different settings tested for prisms. */
#define T8_CMESH_MAX_NUM_XYZ_TREES 1 /**< Maximum number of trees in xyz direction in the cmesh. */
#elif T8_TEST_LEVEL_INT >= 1
#define T8_CMESH_MAX_NUM_OF_TREES 5  /**< Maximum number of trees in the cmesh. */
#define T8_CMESH_MAX_NUM_OF_PRISMS 5 /**< Maximum number of different settings tested for prisms. */
#define T8_CMESH_MAX_NUM_XYZ_TREES 2 /**< Maximum number of trees in xyz direction in the cmesh. */
#else
#define T8_CMESH_MAX_NUM_OF_TREES 10  /**< Maximum number of trees in the cmesh. */
#define T8_CMESH_MAX_NUM_OF_PRISMS 10 /**< Maximum number of different settings tested for prisms. */
#define T8_CMESH_MAX_NUM_XYZ_TREES 3  /**< Maximum number of trees in xyz direction in the cmesh. */
#endif

namespace cmesh_params
{
/** Function to convert an mpi communicator to a string.
 * \param [in] comm The mpi communicator to convert.
 */
std::string
comm_to_string (const sc_MPI_Comm &comm)
{
  int mpi_ret;
  sc_MPI_Comm_compare (comm, sc_MPI_COMM_WORLD, &mpi_ret);
  if (mpi_ret == sc_MPI_SUCCESS) {
    return std::string ("sc_MPI_COMM_WORLD");
  }
  return std::string ("No_String_for_this_communicator");
}

/** Creates a vector filled with sequentially increasing values.
 * \tparam T The type of the vector elements.
 * \param [in] size The number of elements in the resulting vector.
 * \param [in] start The starting value for the sequence.
 * \return A vector of type T containing sequential values starting from \p start. 
 */
template <typename T>
std::vector<T>
filled_vector (const size_t size, T start)
{
  std::vector<T> tmp (size);
  std::iota (tmp.begin (), tmp.end (), start);
  return tmp;
}

/** Define bounds of a unit cube in 3D. */
const double cube_bounds[24] = { 0, 0, 0, 1, 0, 0, 0, 1, 0, 1, 1, 0, 0, 0, 1, 1, 0, 1, 0, 1, 1, 1, 1, 1 };

/** Vector containing parameters for cmesh boundaries. */
std::vector<const double *> boundaries = { cube_bounds };

/** Vector with all possibilities for usage of axis aligned. */
std::vector<int> use_axis_aligned = { 0, 1 };

/** Vector with parameters for a large mesh. */
std::vector<int> large_mesh = filled_vector (T8_CMESH_MAX_NUM_OF_TREES, 1);

/** Vector with parameters for the elements per dimension. */
std::vector<int> elems_per_dim = filled_vector (T8_CMESH_MAX_NUM_XYZ_TREES, 1);

/** Vector containing all mpi communicators used for the tests. */
std::vector<sc_MPI_Comm> my_comms = { sc_MPI_COMM_WORLD };

/** Vector containing all eclasses of dimension 1, 2 and 3. */
std::vector<t8_eclass_t> eclasses = { T8_ECLASS_VERTEX, T8_ECLASS_LINE, T8_ECLASS_QUAD,  T8_ECLASS_TRIANGLE,
                                      T8_ECLASS_HEX,    T8_ECLASS_TET,  T8_ECLASS_PRISM, T8_ECLASS_PYRAMID };

/** Vector containing all eclasses. */
std::vector<t8_eclass_t> all_eclasses
  = { T8_ECLASS_ZERO, T8_ECLASS_VERTEX, T8_ECLASS_LINE,    T8_ECLASS_QUAD,  T8_ECLASS_TRIANGLE, T8_ECLASS_HEX,
      T8_ECLASS_TET,  T8_ECLASS_PRISM,  T8_ECLASS_PYRAMID, T8_ECLASS_COUNT, T8_ECLASS_INVALID };

/** Vector with parameters for broadcast. */
std::vector<int> do_bcast = { 0, 1 };
/** Vector with parameters for partition. */
std::vector<int> partition = { 0, 1 };
/** Vector with parameters for periodic. */
std::vector<int> periodic = { 0, 1 };
/** Currently a dummy vector for examples that have periodic argument but not fully support it yet. */
std::vector<int> no_periodic = { 0 };
/** Vector with parameters for the dimensions. */
std::vector<int> dims = { 1, 2, 3 };
/** Vector with parameters for the number of prisms. */
std::vector<int> num_prisms = filled_vector (T8_CMESH_MAX_NUM_OF_PRISMS, 3);
}  // namespace cmesh_params

#endif /* T8_CMESH_PARAMS_HXX */
