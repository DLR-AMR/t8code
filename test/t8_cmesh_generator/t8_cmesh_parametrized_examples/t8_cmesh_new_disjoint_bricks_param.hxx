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
#include "test/t8_cmesh_generator/t8_cmesh_parametrized_examples/t8_cmesh_params.hxx"
#include "t8_cmesh/t8_cmesh_examples.h"

#ifndef T8_CMESH_NEW_DISJOINT_BRICKS_PARAM_HXX
#define T8_CMESH_NEW_DISJOINT_BRICKS_PARAM_HXX

namespace new_disjoint_bricks
{
std::function<t8_cmesh_t (t8_gloidx_t, t8_gloidx_t, t8_gloidx_t, int, int, int, sc_MPI_Comm)> disjoint_bricks
  = t8_cmesh_new_disjoint_bricks;

std::string
make_param_string (const t8_gloidx_t num_x, const t8_gloidx_t num_y, const t8_gloidx_t num_z, const int num_x_periodic,
                   const int num_y_periodic, const int num_z_periodic, sc_MPI_Comm comm)
{
  std::string delimiter = std::string ("_");
  std::string params = delimiter + std::to_string (num_x) + delimiter + std::to_string (num_y) + delimiter
                       + std::to_string (num_z) + delimiter + std::to_string (num_x_periodic) + delimiter
                       + std::to_string (num_y_periodic) + delimiter + std::to_string (num_z_periodic) + delimiter
                       + cmesh_params::comm_to_string (comm);
  return params;
}

std::function<std::string (const t8_gloidx_t, const t8_gloidx_t, const t8_gloidx_t, const int, const int, const int,
                           const sc_MPI_Comm)>
  make_param_string_wrapper = make_param_string;

example_set *cmesh_example = (example_set *) new cmesh_cartesian_product_params<
  decltype (cmesh_params::elems_per_dim.begin ()), decltype (cmesh_params::elems_per_dim.begin ()),
  decltype (cmesh_params::elems_per_dim.begin ()), decltype (cmesh_params::periodic.begin ()),
  decltype (cmesh_params::periodic.begin ()), decltype (cmesh_params::periodic.begin ()),
  decltype (cmesh_params::my_comms.begin ())> (
  std::make_pair (cmesh_params::elems_per_dim.begin (), cmesh_params::elems_per_dim.end ()),
  std::make_pair (cmesh_params::elems_per_dim.begin (), cmesh_params::elems_per_dim.end ()),
  std::make_pair (cmesh_params::elems_per_dim.begin (), cmesh_params::elems_per_dim.end ()),
  std::make_pair (cmesh_params::periodic.begin (), cmesh_params::periodic.end ()),
  std::make_pair (cmesh_params::periodic.begin (), cmesh_params::periodic.end ()),
  std::make_pair (cmesh_params::periodic.begin (), cmesh_params::periodic.end ()),
  std::make_pair (cmesh_params::my_comms.begin (), cmesh_params::my_comms.end ()), disjoint_bricks,
  make_param_string_wrapper, "t8_cmesh_new_disjoint_brick_");
}  // namespace new_disjoint_bricks

#endif /* T8_CMESH_NEW_DISJOINT_BRICKS_PARAM_HXX */
