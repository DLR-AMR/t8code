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

#include "test/t8_cmesh_generator/t8_gtest_cmesh_cartestian_product.hxx"
#include "test/t8_cmesh_generator/t8_cmesh_parameterized_examples/t8_cmesh_params.hxx"
#include <t8_cmesh/t8_cmesh_examples.h>

#ifndef T8_CMESH_NEW_PRISM_CAKE_PARAM_HXX
#define T8_CMESH_NEW_PRISM_CAKE_PARAM_HXX

namespace new_prism_cake
{
std::function<t8_cmesh_t (sc_MPI_Comm, int)> prism_cake = t8_cmesh_new_prism_cake;

std::string
make_param_string (const sc_MPI_Comm &comm, const int &num_prisms)
{
  std::string delimiter = std::string ("_");
  std::string params = delimiter + cmesh_params::comm_to_string (comm) + delimiter + std::to_string (num_prisms);
  return params;
}

std::function<std::string (const sc_MPI_Comm &, const int &)> make_param_string_wrapper = make_param_string;

example_set *cmesh_example
  = (example_set *) new cmesh_cartesian_product_params<decltype (cmesh_params::my_comms.begin ()),
                                                       decltype (cmesh_params::num_prisms.begin ())> (
    std::make_pair (cmesh_params::my_comms.begin (), cmesh_params::my_comms.end ()),
    std::make_pair (cmesh_params::num_prisms.begin (), cmesh_params::num_prisms.end ()), prism_cake,
    make_param_string_wrapper, "t8_cmesh_new_prism_cake");
}  // namespace new_prism_cake

#endif /* T8_CMESH_NEW_PRISM_CAKE_PARAM_HXX */
