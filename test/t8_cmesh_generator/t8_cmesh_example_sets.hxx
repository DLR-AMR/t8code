/*
This file is part of t8code.
t8code is a C library to manage a collection (a forest) of multiple
connected adaptive space-trees of general element classes in parallel.

Copyright (C) 2023 the developers

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

#include "test/t8_cmesh_generator/t8_cmesh_parametrized_examples/t8_cmesh_new_prism_cake_param.hxx"
#include "test/t8_cmesh_generator/t8_cmesh_parametrized_examples/t8_cmesh_new_from_class_param.hxx"
#include "test/t8_cmesh_generator/t8_gtest_cmesh_cartestian_product.hxx"
#include "test/t8_cmesh_generator/t8_gtest_cmesh_sum_of_sets.hxx"

T8_EXTERN_C_BEGIN ();

namespace cmesh_list
{
std::vector<parameter_cartesian_product *> cart_prod_vec
  = { new_from_class::cmesh_example, new_prism_cake::cmesh_example };

cmesh_sum_of_sets cmesh_sums (cart_prod_vec);
cmesh_sum_of_sets cbegin = cmesh_sums.begin ();
cmesh_sum_of_sets cend = cmesh_sums.end ();
cmesh_sum_of_sets cstep = cmesh_sums.step ();
}  // namespace cmesh_list

#define AllCmeshsParam ::testing::Range (cmesh_list::cbegin, cmesh_list::cend, cmesh_list::cstep)

T8_EXTERN_C_END ();

#endif /* T8_GTEST_CMESH_COMM_CREATOR_HXX */