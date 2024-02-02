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

#ifndef T8_GTEST_CMESH_NUM_ELEM_CREATOR_HXX
#define T8_GTEST_CMESH_NUM_ELEM_CREATOR_HXX

#include <t8_cmesh.h>
#include <t8_cmesh/t8_cmesh_examples.h>
#include <vector>
#include <algorithm>
#include "test/t8_cmesh_generator/t8_gtest_cmesh_creator_base.hxx"

T8_EXTERN_C_BEGIN ();

#define tree_vector std::vector<int> (10)

#define comms \
  std::vector<sc_MPI_Comm> \
  { \
    sc_MPI_COMM_SELF, SC3_MPI_COMM_WORLD \
  }

std::iota (tree_vector.begin (), tree_vector.end (), 2);

cmesh_creator<int, sc_MPI_Comm> num_elem_example (t8_cmesh_new_prism_cake,
                                                  std::make_pair (tree_vector.begin (), tree_vector.end (),
                                                                  std::make_pair (comms.begin (), comms.end ())));

T8_EXTERN_C_END ();

#endif /* T8_GTEST_CMESH_NUM_ELEM_CREATOR_HXX */