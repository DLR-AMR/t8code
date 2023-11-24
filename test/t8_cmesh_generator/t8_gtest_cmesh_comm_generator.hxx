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
#include "test/t8_cmesh_generator/t8_gtest_cmesh_generator_base.hxx"

T8_EXTERN_C_BEGIN ();
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
                                                         t8_cmesh_new_full_hybrid };

struct all_cmeshes_with_comm: public cmesh_generator
{
 public:
  /* Constructor */
  all_cmeshes_with_comm (int creator, sc_MPI_Comm comm): create_func (cmeshes_with_comm[creator])
  {
  }

  all_cmeshes_with_comm (): create_func (cmeshes_with_comm[0])
  {
  }

  /* Copy-Constructor */
  all_cmeshes_with_comm (const all_cmeshes_with_comm &other): create_func (other.create_func)
  {
  }

  /* overloading the < operator for gtest-ranges */
  bool
  operator< (const all_cmeshes_with_comm &other)
  {
    return current_creator < other.current_creator;
  }

  /* The next cmesh is the cmesh with on more element until max_num_trees is reached,
     * then we go to the next constructor. */
  all_cmeshes_with_comm
  operator+ (const all_cmeshes_with_comm &step)
  {
    current_creator += step.current_creator;
    return all_cmeshes_with_comm (current_creator, comm);
  }

  /* To save memory, the cmesh is not created by default */
  void
  create_cmesh ()
  {
    cmesh = create_func (comm);
  }

  void
  get_step (all_cmeshes_with_comm *step)
  {
    step->current_creator = 1;
    step->comm = sc_MPI_COMM_WORLD;
  }

  void
  get_first (all_cmeshes_with_comm *first)
  {
    first->current_creator = 0;
    first->comm = sc_MPI_COMM_WORLD;
    first->cmesh = NULL;
    first->create_func = cmeshes_with_comm[0];
  }

  void
  set_last ()
  {
    current_creator = cmeshes_with_comm.size () - 1;
  }

  int
  is_at_last ()
  {
    return (long unsigned int) current_creator == cmeshes_with_comm.size () - 1;
  }

  /* Destruktor */
  ~all_cmeshes_with_comm ()
  {
    /* unref the cmesh only if it has been created */
    if (cmesh != NULL) {
      t8_cmesh_unref (&cmesh);
    }
  }

  t8_cmesh_w_comm create_func = cmeshes_with_comm[0];
};

T8_EXTERN_C_END ();

#endif /* T8_GTEST_CMESH_COMM_GENERATOR_HXX */