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

/**
 * File to store the parameters that are used for our parameterized cmesh-tests
 */

#include <vector>
#include <t8_eclass.h>

#ifndef T8_CMESH_PARAMS_HXX
#define T8_CMESH_PARAMS_HXX

namespace cmesh_params
{

std::vector<sc_MPI_Comm> my_comms = { sc_MPI_COMM_WORLD };
std::vector<t8_eclass_t> eclasses = { T8_ECLASS_VERTEX, T8_ECLASS_LINE, T8_ECLASS_QUAD,  T8_ECLASS_TRIANGLE,
                                      T8_ECLASS_HEX,    T8_ECLASS_TET,  T8_ECLASS_PRISM, T8_ECLASS_PYRAMID };

std::vector<t8_eclass_t> all_eclasses
  = { T8_ECLASS_ZERO, T8_ECLASS_VERTEX, T8_ECLASS_LINE,    T8_ECLASS_QUAD,  T8_ECLASS_TRIANGLE, T8_ECLASS_HEX,
      T8_ECLASS_TET,  T8_ECLASS_PRISM,  T8_ECLASS_PYRAMID, T8_ECLASS_COUNT, T8_ECLASS_INVALID };

std::vector<int> num_prisms = { 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13 };
}  // namespace cmesh_params

#endif /* T8_CMESH_PARAMS_HXX */