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

#ifndef T8_GTEST_CMESH_GENERATOR_CXX_HXX
#define T8_GTEST_CMESH_GENERATOR__CXX_HXX

#include "test/t8_cmesh_generator/t8_gtest_cmesh_generator_base.hxx"
#include "test/t8_cmesh_generator/t8_gtest_cmesh_comm_generator.hxx"
#include "test/t8_cmesh_generator/t8_gtest_cmesh_num_elem_generator.hxx"

T8_EXTERN_C_BEGIN ();

typedef enum cmesh_types {
  CMESH_ZERO = 0,

  CMESH_WITH_COMM = CMESH_ZERO,
  CMESH_WITH_NUM_TREES = 1,

  CMESH_NUM_TYPES = 2
} cmesh_types_t;

class cmesh_generator_cxx {
 public:
  cmesh_generator_cxx (): current_generator (CMESH_ZERO)
  {
    generators = T8_ALLOC_ZERO (cmesh_generator *, CMESH_NUM_TYPES);
    t8_refcount_init (&rc);

    generators[CMESH_WITH_COMM] = new all_cmeshes_with_comm ();
    generators[CMESH_WITH_NUM_TREES] = new all_cmeshes_with_num_elem ();
  }

  /* Avoid code duplication here*/
  cmesh_generator_cxx (int generator)
  {
    T8_ASSERT (generator < (int) CMESH_NUM_TYPES);
    generators = T8_ALLOC_ZERO (cmesh_generator *, CMESH_NUM_TYPES);
    t8_refcount_init (&rc);

    generators[CMESH_WITH_COMM] = new all_cmeshes_with_comm ();
    generators[CMESH_WITH_NUM_TREES] = new all_cmeshes_with_num_elem ();
    current_generator = generator;
  }

  cmesh_generator_cxx (cmesh_generator_cxx &other): current_generator (other.current_generator)
  {
    T8_ASSERT (other.generators != NULL);
    sc_refcount_ref (&other.rc);

    generators = other.generators;
  }

  void
  cmesh_generator_cxx_destroy ()
  {
    for (int igen = CMESH_ZERO; igen < CMESH_NUM_TYPES; igen++) {
      if (generators != NULL) {
        delete generators[igen];
      }
    }
    T8_FREE (generators);
  }

  void
  create_cmesh ()
  {
    generators[current_generator]->create_cmesh ();
  }

  t8_cmesh_t
  get_cmesh ()
  {
    return generators[current_generator]->get_cmesh ();
  }

  bool
  operator< (const cmesh_generator_cxx &other)
  {
    if (current_generator == other.current_generator) {
      return generators[current_generator] < other.generators[current_generator];
    }
    else {
      return current_generator < other.current_generator;
    }
  }

  cmesh_generator_cxx
  operator+ (const cmesh_generator_cxx &step)
  {
    cmesh_generator_cxx tmp (*this);
    if (generators[current_generator]->is_at_last ()) {
      tmp.current_generator++;
      tmp.generators[current_generator]->set_first ();
    }
    else {
      tmp.generators[current_generator]->addition (step.generators[current_generator]);
    }
    return tmp;
  }

  void
  set_last ()
  {
    current_generator = CMESH_NUM_TYPES - 1;
    generators[current_generator]->set_last ();
  }

  ~cmesh_generator_cxx ()
  {
    T8_ASSERT (generators != NULL);
    if (sc_refcount_unref (&rc)) {
      cmesh_generator_cxx_destroy ();
    }
  }

  int current_generator;
  sc_refcount_t rc;
  cmesh_generator **generators;
};

T8_EXTERN_C_END ();
#endif /* T8_GTEST_CMESH_GENERATOR_CXX_HXX */