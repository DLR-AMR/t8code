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
#include "test/t8_cmesh_generator/t8_gtest_cmesh_creator_base.hxx"

T8_EXTERN_C_BEGIN ();

#define MAX_NUM_TREES 5

/* A function creating a cmesh getting a communicator and a number of elements to create */
typedef t8_cmesh_t (*t8_cmesh_with_num_trees) (sc_MPI_Comm comm, int num_trees);

/* List of all functions that have a num_trees-parameter*/
const std::vector<t8_cmesh_with_num_trees> cmeshes_with_num_trees
  = { t8_cmesh_new_long_brick_pyramid, t8_cmesh_new_prism_cake, t8_cmesh_new_pyramid_cake };

/**
 * A class that creates all cmeshes using a number of trees and a communicator
 */
class all_cmeshes_with_num_trees: public cmesh_creator {
 public:
  /* overloading the < operator for gtest-ranges */
  bool
  operator< (const cmesh_creator &other)
  {
    if (current_creator == other.current_creator) {
      return num_trees < other.num_trees;
    }
    else {
      return current_creator < other.current_creator;
    }
  }

  /* The next cmesh is the cmesh with on more element until max_num_trees is reached,
     * then we go to the next constructor. */
  void
  addition (const std::shared_ptr<cmesh_creator> step)
  {
    if (num_trees + step->num_trees > MAX_NUM_TREES) {
      num_trees = (num_trees + step->num_trees) % (MAX_NUM_TREES + 1);
      current_creator++;
      T8_ASSERT ((unsigned long int) current_creator < cmeshes_with_num_trees.size ());
    }
    else {
      num_trees += step->num_trees;
    }
    if (current_creator > 0 && num_trees <= 2) {
      num_trees = 3;
    }
  }

  /* To save memory, the cmesh is not created by default */
  void
  create_cmesh ()
  {
    cmesh = cmeshes_with_num_trees[current_creator](comm, num_trees);
  }

  void
  set_first ()
  {
    current_creator = 0;
    num_trees = 1;
  }

  void
  set_last ()
  {
    current_creator = cmeshes_with_num_trees.size () - 1;
    num_trees = MAX_NUM_TREES;
  }

  bool
  is_at_last ()
  {
    return (long unsigned int) current_creator == cmeshes_with_num_trees.size () - 1 && num_trees == MAX_NUM_TREES;
  }
  /* Destructor */
  ~all_cmeshes_with_num_trees ()
  {
    /* unref the cmesh only if it has been created */
    if (cmesh != NULL) {
      t8_cmesh_unref (&cmesh);
    }
  }
};

T8_EXTERN_C_END ();

#endif /* T8_GTEST_CMESH_NUM_ELEM_CREATOR_HXX */