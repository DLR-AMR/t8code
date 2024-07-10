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

#ifndef T8_CMESH_NEW_COMM
#define T8_CMESH_NEW_COMM

#include "test/t8_cmesh_generator/t8_gtest_cmesh_cartestian_product.hxx"
#include "test/t8_cmesh_generator/t8_cmesh_parametrized_examples/t8_cmesh_params.hxx"
#include <t8_cmesh/t8_cmesh_examples.h>

namespace new_cmesh_comm
{
std::string
make_param_string (const sc_MPI_Comm &comm)
{
  std::string delimiter = std::string ("_");
  std::string params = delimiter + cmesh_params::comm_to_string (comm);
  return params;
}

std::function<std::string (const sc_MPI_Comm)> print_function = make_param_string;

std::vector<std::function<t8_cmesh_t (sc_MPI_Comm)>> cmesh_functions = { t8_cmesh_new_periodic_tri,
                                                                         t8_cmesh_new_periodic_hybrid,
                                                                         t8_cmesh_new_periodic_line_more_trees,
                                                                         t8_cmesh_new_line_zigzag,
                                                                         t8_cmesh_new_prism_deformed,
                                                                         t8_cmesh_new_pyramid_deformed,
                                                                         t8_cmesh_new_prism_cake_funny_oriented,
                                                                         t8_cmesh_new_prism_geometry,
                                                                         t8_cmesh_new_tet_orientation_test,
                                                                         t8_cmesh_new_hybrid_gate,
                                                                         t8_cmesh_new_hybrid_gate_deformed,
                                                                         t8_cmesh_new_full_hybrid };

std::vector<std::string> names = { "t8_cmesh_new_periodic_tri_",
                                   "t8_cmesh_new_periodic_hybrid_",
                                   "t8_cmesh_new_periodic_line_more_trees_",
                                   "t8_cmesh_new_line_zigzag_",
                                   "t8_cmesh_new_prism_deformed_",
                                   "t8_cmesh_new_pyramid_deformed_",
                                   "t8_cmesh_new_prism_cake_funny_oriented_",
                                   "t8_cmesh_new_prism_geometry_",
                                   "t8_cmesh_new_tet_orientation_test_",
                                   "t8_cmesh_new_hybrid_gate_",
                                   "t8_cmesh_new_hybrid_gate_deformed_",
                                   "t8_cmesh_new_full_hybrid_" };

example_set *cmesh_example
  = (example_set *) new cmesh_cartesian_product_params<decltype (cmesh_params::my_comms.begin ())> (
    std::make_pair (cmesh_params::my_comms.begin (), cmesh_params::my_comms.end ()), cmesh_functions, print_function,
    names);
}  // namespace new_cmesh_comm

#endif /* T8_CMESH_NEW_COMM */