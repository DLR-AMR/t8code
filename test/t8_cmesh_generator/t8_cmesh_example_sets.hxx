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

#ifndef T8_GTEST_CMESH_COMM_CREATOR_HXX
#define T8_GTEST_CMESH_COMM_CREATOR_HXX

#include <vector>

#include "test/t8_cmesh_generator/t8_cmesh_parameterized_examples/t8_cmesh_new_prism_cake_param.hxx"
#include "test/t8_cmesh_generator/t8_cmesh_parameterized_examples/t8_cmesh_new_from_class_param.hxx"
#include "test/t8_cmesh_generator/t8_cmesh_parameterized_examples/t8_cmesh_new_bigmesh_param.hxx"
#include "test/t8_cmesh_generator/t8_cmesh_parameterized_examples/t8_cmesh_new_disjoint_bricks_param.hxx"
#include "test/t8_cmesh_generator/t8_cmesh_parameterized_examples/t8_cmesh_new_comm.hxx"
#include "test/t8_cmesh_generator/t8_cmesh_parameterized_examples/t8_cmesh_new_hypercube_pad.hxx"
#include "test/t8_cmesh_generator/t8_cmesh_parameterized_examples/t8_cmesh_new_hypercube_param.hxx"
#include "test/t8_cmesh_generator/t8_cmesh_parameterized_examples/t8_cmesh_new_empty.hxx"
#include "test/t8_cmesh_generator/t8_cmesh_parameterized_examples/t8_cmesh_new_periodic.hxx"
#include "test/t8_cmesh_generator/t8_gtest_cmesh_cartestian_product.hxx"
#include "test/t8_cmesh_generator/t8_gtest_cmesh_sum_of_sets.hxx"
#include "test/t8_gtest_schemes.hxx"

T8_EXTERN_C_BEGIN ();

/**
 * lambda to pass to an INSTANTIATE_TEST_SUITE_P to print the current cmesh_example_base
 * 
 */
auto pretty_print_base_example = [] (const testing::TestParamInfo<cmesh_example_base *> &info) {
  std::string name;
  info.param->param_to_string (name);
  return name;
};

auto pretty_print_base_example_scheme = [] (const testing::TestParamInfo<std::tuple<int, cmesh_example_base *>> &info) {
  std::string name;
  std::get<1> (info.param)->param_to_string (name);
  name += std::string ("_") + t8_scheme_to_string[std::get<0> (info.param)];
  return name;
};

namespace cmesh_list
{
std::vector<example_set *> cart_prod_vec = { new_from_class::cmesh_example,
                                             new_prism_cake::cmesh_example,
                                             new_bigmesh::cmesh_example,
                                             new_cmesh_comm::cmesh_example,
                                             new_hypercube_pad::cmesh_example_non_periodic_boundaries,
                                             new_hypercube_pad::cmesh_example_periodic_boundaries,
                                             new_hypercube_cmesh::cmesh_example,
                                             new_hypercube_cmesh::cmesh_example_pyra,
                                             new_disjoint_bricks::cmesh_example,
                                             new_empty::cmesh_example,
                                             new_periodic::cmesh_example };

cmesh_sum_of_sets cmesh_sums (cart_prod_vec);

}  // namespace cmesh_list

#define AllCmeshsParam \
  ::testing::ValuesIn (cmesh_list::cmesh_sums.cmesh_examples.begin (), cmesh_list::cmesh_sums.cmesh_examples.end ())

T8_EXTERN_C_END ();

#endif /* T8_GTEST_CMESH_COMM_CREATOR_HXX */
