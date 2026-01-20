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

#ifndef T8_CMESH_NEW_BIGMESH_PARAM_HXX
#define T8_CMESH_NEW_BIGMESH_PARAM_HXX

#include "t8_cmesh_generator/t8_gtest_cmesh_cartestian_product.hxx"
#include "t8_cmesh_generator/t8_cmesh_parameterized_examples/t8_cmesh_params.hxx"
#include "t8_cmesh/t8_cmesh_examples.h"
#include <t8_eclass.h>

namespace new_bigmesh
{
std::function<t8_cmesh_t (t8_eclass_t, int, sc_MPI_Comm)> bigmesh = t8_cmesh_new_bigmesh;

std::string
make_param_string (const t8_eclass_t eclass, const int num_trees, const sc_MPI_Comm comm)
{
  std::string delimiter = std::string ("_");
  std::string params = delimiter + t8_eclass_to_string[eclass] + delimiter + std::to_string (num_trees) + delimiter
                       + cmesh_params::comm_to_string (comm);
  return params;
}

std::function<std::string (const t8_eclass_t, const int, const sc_MPI_Comm)> make_param_string_wrapper
  = make_param_string;

example_set *cmesh_example
  = (example_set *) new cmesh_cartesian_product_params<decltype (cmesh_params::eclasses.begin ()),
                                                       decltype (cmesh_params::large_mesh.begin ()),
                                                       decltype (cmesh_params::my_comms.begin ())> (
    std::make_pair (cmesh_params::eclasses.begin (), cmesh_params::eclasses.end ()),
    std::make_pair (cmesh_params::large_mesh.begin (), cmesh_params::large_mesh.end ()),
    std::make_pair (cmesh_params::my_comms.begin (), cmesh_params::my_comms.end ()), bigmesh, make_param_string_wrapper,
    "t8_cmesh_new_bigmesh");
}  // namespace new_bigmesh

#endif /* T8_CMESH_NEW_BIGMESH_PARAM_HXX */
