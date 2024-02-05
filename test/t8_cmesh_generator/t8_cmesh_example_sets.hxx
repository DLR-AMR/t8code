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
#include <t8_cmesh/t8_cmesh_examples.h>
#include "test/t8_cmesh_generator/t8_gtest_cmesh_cartestian_product.hxx"
#include "test/t8_cmesh_generator/t8_gtest_cmesh_sum_of_sets.hxx"
#include <t8_cmesh/t8_cmesh_geometry.h>
#include <t8_geometry/t8_geometry_implementations/t8_geometry_linear.h>
#include <t8_geometry/t8_geometry_implementations/t8_geometry_linear_axis_aligned.h>
#include <t8_geometry/t8_geometry_base.h>

T8_EXTERN_C_BEGIN ();

namespace cmesh_list
{
std::vector<sc_MPI_Comm> my_comms = { sc_MPI_COMM_WORLD };
std::vector<t8_eclass_t> eclasses = { T8_ECLASS_ZERO, T8_ECLASS_LINE, T8_ECLASS_TRIANGLE };

std::function<t8_cmesh_t (t8_eclass_t, sc_MPI_Comm)> new_from_class_wrapper = t8_cmesh_new_from_class;
parameter_cartesian_product *new_form_class_prod
  = (parameter_cartesian_product *) new cmesh_parameter_combinations<decltype (eclasses.begin ()),
                                                                     decltype (my_comms.begin ())> (
    std::make_pair (eclasses.begin (), eclasses.end ()), std::make_pair (my_comms.begin (), my_comms.end ()),
    new_from_class_wrapper, "t8_cmesh_new_from_class");

std::vector<int> num_prisms = { 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13 };
std::function<t8_cmesh_t (sc_MPI_Comm, int)> prism_cake = t8_cmesh_new_prism_cake;
parameter_cartesian_product *new_prism_cake
  = (parameter_cartesian_product *) new cmesh_parameter_combinations<decltype (my_comms.begin ()),
                                                                     decltype (num_prisms.begin ())> (
    std::make_pair (my_comms.begin (), my_comms.end ()), std::make_pair (num_prisms.begin (), num_prisms.end ()),
    prism_cake, "t8_cmesh_new_prism_cake");

std::vector<parameter_cartesian_product *> cart_prod_vec = { new_form_class_prod, new_prism_cake };

cmesh_sum_of_sets cmesh_sums (cart_prod_vec);
cmesh_sum_of_sets cbegin = cmesh_sums.begin ();
cmesh_sum_of_sets cend = cmesh_sums.end ();
cmesh_sum_of_sets cstep = cmesh_sums.step ();
}  // namespace cmesh_list

#define AllCmeshsParam ::testing::Range (cmesh_list::cbegin, cmesh_list::cend, cmesh_list::cstep)

T8_EXTERN_C_END ();

#endif /* T8_GTEST_CMESH_COMM_CREATOR_HXX */