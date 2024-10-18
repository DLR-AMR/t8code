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
#include <t8_cmesh/t8_cmesh_examples.h>
#include <t8_eclass.h>

#ifndef T8_CMESH_NEW_FROM_CLASS_PARAM
#define T8_CMESH_NEW_FROM_CLASS_PARAM

namespace new_from_class
{

std::string
make_param_string (const t8_eclass_t &eclass, const sc_MPI_Comm &comm)
{
  std::string delimiter = std::string ("_");
  std::string params
    = delimiter + std::string (t8_eclass_to_string[eclass]) + delimiter + cmesh_params::comm_to_string (comm);
  return params;
}

std::function<std::string (const t8_eclass_t &, const sc_MPI_Comm &)> print_function = make_param_string;

std::function<t8_cmesh_t (t8_eclass_t, sc_MPI_Comm)> new_from_class_wrapper = t8_cmesh_new_from_class;

example_set *cmesh_example
  = (example_set *) new cmesh_cartesian_product_params<decltype (cmesh_params::eclasses.begin ()),
                                                       decltype (cmesh_params::my_comms.begin ())> (
    std::make_pair (cmesh_params::eclasses.begin (), cmesh_params::eclasses.end ()),
    std::make_pair (cmesh_params::my_comms.begin (), cmesh_params::my_comms.end ()), new_from_class_wrapper,
    print_function, "t8_cmesh_new_from_class_");

}  // namespace new_from_class

#endif /* T8_CMESH_NEW_FROM_CLASS_PARAM */
