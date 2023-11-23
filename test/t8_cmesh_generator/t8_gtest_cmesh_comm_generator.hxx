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

#ifndef T8_GTEST_CMESH_COMM_GENERATOR_HXX
#define T8_GTEST_CMESH_COMM_GENERATOR_HXX

#include <vector>
#include <t8_cmesh/t8_cmesh_examples.h>

/* A function creating a cmesh getting a communicator */
typedef t8_cmesh_t (*t8_cmesh_w_comm) (sc_MPI_Comm comm);

/* List of all functions that have a num_elems-parameter*/
const std::vector<t8_cmesh_w_comm> cmeshes_with_comm = { t8_cmesh_new_periodic_tri,
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
                                                         t8_cmesh_new_full_hybrid,
                                                         t8_cmesh_new_pyramid_cake };

#endif /* T8_GTEST_CMESH_COMM_GENERATOR_HXX */