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

#ifndef T8_CMESH_NEW_HYPERCUBE_PARAM_HXX
#define T8_CMESH_NEW_HYPERCUBE_PARAM_HXX

#include <t8_eclass.h>
#include "test/t8_cmesh_generator/t8_gtest_cmesh_cartestian_product.hxx"
#include "test/t8_cmesh_generator/t8_cmesh_parametrized_examples/t8_cmesh_params.hxx"
#include <t8_cmesh/t8_cmesh_examples.h>

namespace new_hypercube_cmesh
{
std::string
make_param_string (const t8_eclass_t eclass, const sc_MPI_Comm comm, const int do_bcast, const int do_partition,
                   const int periodic)
{
  std::string delimiter = std::string ("_");
  std::string bcast = do_bcast ? std::string ("bcast") : std::string ("noBcast");
  std::string partition = do_partition ? std::string ("partition") : std::string ("noPartition");
  std::string periodic_string = periodic ? std::string ("periodic") : std::string ("noPeriodic");
  std::string params
    = delimiter + t8_eclass_to_string[eclass] + delimiter + bcast + delimiter + partition + delimiter + periodic_string;

  return params;
}

std::function<std::string (const t8_eclass_t, const sc_MPI_Comm, const int, const int, const int)> param_to_string
  = make_param_string;

std::vector<t8_eclass_t> periodic_eclasses = { T8_ECLASS_VERTEX, T8_ECLASS_LINE, T8_ECLASS_QUAD, T8_ECLASS_TRIANGLE,
                                               T8_ECLASS_HEX,    T8_ECLASS_TET,  T8_ECLASS_PRISM };

std::function<t8_cmesh_t (t8_eclass_t, sc_MPI_Comm, int, int, int)> cmesh_wrapper = t8_cmesh_new_hypercube;

std::vector<t8_eclass_t> nonperiodic_eclasses = { T8_ECLASS_PYRAMID };

example_set *cmesh_example = (example_set *) new cmesh_cartesian_product_params<
  decltype (periodic_eclasses.begin ()), decltype (cmesh_params::my_comms.begin ()),
  decltype (cmesh_params::do_bcast.begin ()), decltype (cmesh_params::partition.begin ()),
  decltype (cmesh_params::periodic.begin ())> (
  std::make_pair (periodic_eclasses.begin (), periodic_eclasses.end ()),
  std::make_pair (cmesh_params::my_comms.begin (), cmesh_params::my_comms.end ()),
  std::make_pair (cmesh_params::do_bcast.begin (), cmesh_params::do_bcast.end ()),
  std::make_pair (cmesh_params::no_partition.begin (), cmesh_params::no_partition.end ()),
  std::make_pair (cmesh_params::periodic.begin (), cmesh_params::periodic.end ()), cmesh_wrapper, param_to_string,
  "t8_cmesh_new_hypercube_");

example_set *cmesh_example_pyra = (example_set *) new cmesh_cartesian_product_params<
  decltype (periodic_eclasses.begin ()), decltype (cmesh_params::my_comms.begin ()),
  decltype (cmesh_params::do_bcast.begin ()), decltype (cmesh_params::partition.begin ()),
  decltype (cmesh_params::no_periodic.begin ())> (
  std::make_pair (nonperiodic_eclasses.begin (), nonperiodic_eclasses.end ()),
  std::make_pair (cmesh_params::my_comms.begin (), cmesh_params::my_comms.end ()),
  std::make_pair (cmesh_params::do_bcast.begin (), cmesh_params::do_bcast.end ()),
  std::make_pair (cmesh_params::no_partition.begin (), cmesh_params::no_partition.end ()),
  std::make_pair (cmesh_params::no_periodic.begin (), cmesh_params::no_periodic.end ()), cmesh_wrapper, param_to_string,
  "t8_cmesh_new_hypercube_");
}  // namespace new_hypercube_cmesh

#endif /* T8_CMESH_NEW_HYPERCUBE_PARAM_HXX */
