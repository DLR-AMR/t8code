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

#ifndef T8_GTEST_cmesh_sum_of_sets_HXX
#define T8_GTEST_cmesh_sum_of_sets_HXX

#include "test/t8_cmesh_generator/t8_gtest_cmesh_cartestian_product.hxx"

T8_EXTERN_C_BEGIN ();

/**
 * A class that holds multiple ways to create a cmesh.
 * 
 */
class cmesh_sum_of_sets {
 public:
  cmesh_sum_of_sets () {};
  /**
   * Construct a new cmesh sum of sets object, that will generate cmeshes given by \ref cmesh_set
   * 
   * \param[in] cmesh_set A vector of \ref parameter_cartesian_product 
   */
  cmesh_sum_of_sets (std::vector<example_set *> cmesh_set)
  {
    for (size_t icreator = 0; icreator < cmesh_set.size (); icreator++) {
      cmesh_examples.insert (cmesh_examples.end (), cmesh_set[icreator]->example_all_combination.begin (),
                             cmesh_set[icreator]->example_all_combination.end ());
    }
  }

  /**
   * Destroy the cmesh generator cxx object
   * 
   */
  ~cmesh_sum_of_sets ()
  {
  }

 public:
  std::vector<cmesh_example_base *> cmesh_examples;
};

T8_EXTERN_C_END ();
#endif /* T8_GTEST_cmesh_sum_of_sets_HXX */