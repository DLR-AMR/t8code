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

#ifndef T8_GTEST_CMESH_GENERATOR_HXX
#define T8_GTEST_CMESH_GENERATOR_HXX

#include "test/t8_cmesh_generator/t8_gtest_cmesh_creator_base.hxx"
#include "test/t8_cmesh_generator/t8_gtest_cmesh_comm_creator.hxx"
#include "test/t8_cmesh_generator/t8_gtest_cmesh_num_elem_creator.hxx"

T8_EXTERN_C_BEGIN ();

typedef enum cmesh_types {
  CMESH_ZERO = 0,

  CMESH_WITH_COMM = CMESH_ZERO,
  CMESH_WITH_NUM_TREES = 1,

  CMESH_NUM_TYPES = 2
} cmesh_types_t;

/**
 * A class that can create all cmeshes from all cmesh_generators
 * 
 */
class cmesh_generator {
 public:
  /**
   * Construct a new cmesh generator cxx object
   * 
   */
  cmesh_generator (): current_generator (CMESH_ZERO)
  {
    generators[CMESH_WITH_COMM] = std::shared_ptr<all_cmeshes_with_comm> (new all_cmeshes_with_comm ());
    generators[CMESH_WITH_NUM_TREES] = std::shared_ptr<all_cmeshes_with_num_trees> (new all_cmeshes_with_num_trees ());
  }

  /**
   * Construct a new cmesh generator cxx object given a specific generator
   * 
   * \param[in] generator which cmesh_generator to use, has to be smaller than CMESH_NUM_TYPES
   */
  cmesh_generator (int generator): cmesh_generator ()
  {
    T8_ASSERT (generator < (int) CMESH_NUM_TYPES);
    current_generator = generator;
  }

  /**
   * Copy-constructor
   * 
   * \param other[in] cmesh_generator to copy from
   */
  cmesh_generator (const cmesh_generator &other): current_generator (other.current_generator)
  {
    T8_ASSERT (other.generators != NULL);
    for (int igen = CMESH_ZERO; igen < CMESH_NUM_TYPES; igen++) {
      T8_ASSERT (other.generators[igen] != NULL);
    }

    for (int igen = CMESH_ZERO; igen < CMESH_NUM_TYPES; igen++) {
      generators[igen] = other.generators[igen];
    }
  }

  /**
   * Call the creator of the current generator
   */
  void
  create_cmesh ()
  {
    generators[current_generator]->create_cmesh ();
  }

  /**
   * Unref the cmesh of the current generator
   * 
   */
  void
  unref_cmesh ()
  {
    generators[current_generator]->unref_cmesh ();
  }

  /**
   * Get the cmesh of the current generator
   * 
   * \return If created, the cmesh of the current generator, NULL otherwise. 
   */
  t8_cmesh_t
  get_cmesh ()
  {
    return generators[current_generator]->get_cmesh ();
  }

  /**
   * Compare two cmesh_generator. If the current_generator is equal, the generators are compared. 
   * If not, the object with the smaller current_generator is considered smaller. 
   * 
   * \param[in] other the cmesh_generator to compare with.
   * \return true if both are equal
   * \return false ow
   */
  bool
  operator< (const cmesh_generator &other)
  {
    if (current_generator == other.current_generator) {
      return generators[current_generator] < other.generators[current_generator];
    }
    else {
      return current_generator < other.current_generator;
    }
  }

  /**
   * + operator to be able to use ::testing:Range from the GoogleTestSuite. 
   * 
   * \param[in] step cmesh_generator describing by how far to step forward
   * \return cmesh_generator 
   */
  cmesh_generator
  operator+ (const cmesh_generator &step)
  {
    cmesh_generator tmp (*this);
    if (generators[current_generator]->is_at_last ()) {
      tmp.current_generator++;
      tmp.generators[current_generator]->set_first ();
    }
    else {
      tmp.generators[current_generator]->addition (step.generators[current_generator]);
    }
    return cmesh_generator (tmp);
  }

  /**
   * Set the object to the last creator of the last generator
   * 
   */
  cmesh_generator
  get_last ()
  {
    cmesh_generator tmp;
    tmp.current_generator = CMESH_NUM_TYPES - 1;
    tmp.generators[current_generator]->set_last ();
    return tmp;
  }

  /**
   * Destroy the cmesh generator cxx object
   * 
   */
  ~cmesh_generator ()
  {
  }

  /* The generator that is active. */
  int current_generator;
  /* Holds all generators. 
   * TODO: only use a single shared_ptr to a cmesh_generator and update */
  std::shared_ptr<cmesh_creator> generators[CMESH_NUM_TYPES];
};

T8_EXTERN_C_END ();
#endif /* T8_GTEST_CMESH_GENERATOR_HXX */