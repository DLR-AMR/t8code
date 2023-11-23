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

/**
 * \file This file provides a generator for cmeshes for googletest
 * 
 */

#ifndef T8_GTEST_CMESH_GENERATOR_HXX
#define T8_GTEST_CMESH_GENERATOR_HXX

#include <t8_cmesh.h>
#include <t8_cmesh/t8_cmesh_examples.h>
#include <vector>

#ifdef T8_ENABLE_LESS_TESTS
#define MAX_NUM_TREES 5
#else
#define MAX_NUM_TREES 10

#endif

/* A function creating a cmesh getting a communicator and a number of elements to create */
typedef t8_cmesh_t (*t8_cmesh_w_num_elem_create) (sc_MPI_Comm comm, int num_elems);

/* List of all functions that have a num_elems-parameter*/
const std::vector<t8_cmesh_w_num_elem_create> cmeshes_with_num_elems
  = { t8_cmesh_new_long_brick_pyramid, t8_cmesh_new_prism_cake, t8_cmesh_new_pyramid_cake };

class all_cmeshes_with_num_elem {
 public:
  /* Constructor */
  all_cmeshes_with_num_elem (int creator, sc_MPI_Comm comm, int num_elems)
    : num_elems (num_elems), current_creator (creator), create_func (cmeshes_with_num_elems[creator]), comm (comm)
  {
  }

  //all_cmeshes_with_num_elem() : num_elems(1), current_creator(0), create_func(cmeshes_with_num_elems[0]),
  //comm(sc_MPI_COMM_WORLD){}

  /* Copy-Constructor */
  all_cmeshes_with_num_elem (const all_cmeshes_with_num_elem &other)
    : num_elems (other.num_elems), current_creator (other.current_creator), create_func (other.create_func),
      comm (other.comm), cmesh (NULL)
  {
  }

  /* overloading the < operator for gtest-ranges */
  bool
  operator< (const all_cmeshes_with_num_elem &other)
  {
    if (current_creator == other.current_creator) {
      return num_elems < other.num_elems;
    }
    else {
      return current_creator < other.current_creator;
    }
  }

  /* The next cmesh is the cmesh with on more element until max_num_trees is reached,
     * then we go to the next constructor. */
  all_cmeshes_with_num_elem
  operator+ (const all_cmeshes_with_num_elem &step)
  {
    if (num_elems + step.num_elems > MAX_NUM_TREES) {
      num_elems = (num_elems + step.num_elems) % (MAX_NUM_TREES + 1);
      current_creator++;
    }
    else {
      num_elems += step.num_elems;
    }
    if (current_creator > 0 && num_elems <= 2)
      num_elems = 3;
    return all_cmeshes_with_num_elem (current_creator, comm, num_elems);
  }

  /* To save memory, the cmesh is not created by default */
  void
  create_cmesh ()
  {
    cmesh = create_func (comm, num_elems);
  }

  /* The cmesh is referenced and returned */
  t8_cmesh_t
  get_cmesh ()
  {
    t8_cmesh_ref (cmesh);
    return cmesh;
  }

  /* Destruktor */
  ~all_cmeshes_with_num_elem ()
  {
    /* unref the cmesh only if it has been created */
    if (cmesh != NULL) {
      t8_cmesh_unref (&cmesh);
    }
  }
  int num_elems;
  int current_creator;
  t8_cmesh_w_num_elem_create create_func;
  sc_MPI_Comm comm;
  t8_cmesh_t cmesh = NULL;
};

#endif /* T8_GTEST_CMESH_GENERATOR_HXX */