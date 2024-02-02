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
#include "test/t8_cmesh_generator/t8_gtest_cmesh_creator_base.hxx"

T8_EXTERN_C_BEGIN ();
#if 0
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

/**
 * A class to generate all cmeshes using only a communicator
 * 
 */
class all_cmeshes_with_comm: public cmesh_creator {
 public:
  /* overloading the < operator for gtest-ranges */
  bool
  operator< (const cmesh_creator &other)
  {
    return current_creator < other.current_creator;
  }

  /* The next cmesh is the cmesh with on more element until max_num_trees is reached,
     * then we go to the next constructor. */
  void
  addition (const std::shared_ptr<cmesh_creator> step)
  {
    current_creator += 1;
    T8_ASSERT ((unsigned long int) current_creator < cmeshes_with_comm.size ());
  }

  /* To save memory, the cmesh is not created by default */
  void
  create_cmesh ()
  {
    cmesh = cmeshes_with_comm[current_creator](comm);
  }

  void
  set_first ()
  {
    current_creator = 0;
  }

  void
  set_last ()
  {
    current_creator = cmeshes_with_comm.size () - 1;
  }

  bool
  is_at_last ()
  {
    return (long unsigned int) current_creator == cmeshes_with_comm.size () - 1;
  }

  /* Destructor */
  ~all_cmeshes_with_comm ()
  {
    /* unref the cmesh only if it has been created */
    if (cmesh != NULL) {
      t8_cmesh_unref (&cmesh);
    }
  }
};
#endif
T8_EXTERN_C_END ();

#endif /* T8_GTEST_CMESH_COMM_CREATOR_HXX */