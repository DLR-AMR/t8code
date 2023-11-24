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
#include <vector>
#include "test/t8_cmesh_generator/t8_gtest_cmesh_generator_cxx.hxx"

T8_EXTERN_C_BEGIN ();

struct generate_all_cmeshes
{
 public:
  generate_all_cmeshes (): current_cmesh_type (CMESH_ZERO)
  {
    generator = cmesh_generator_init ();
  };

  generate_all_cmeshes (const int cmesh_type): current_cmesh_type (cmesh_type)
  {
    generator = cmesh_generator_init ();
  }

  generate_all_cmeshes (const int cmesh_type, cmesh_generator_cxx *gen)
  {
    cmesh_generator_cxx_ref (gen);
    generator = gen;
    current_cmesh_type = cmesh_type;
  }

  generate_all_cmeshes (const generate_all_cmeshes &other): current_cmesh_type (other.current_cmesh_type)
  {
    cmesh_generator_cxx_ref (other.generator);
    generator = other.generator;
  };

  bool
  operator< (const generate_all_cmeshes &other)
  {
    if (current_cmesh_type == other.current_cmesh_type) {
      return generator->generators[current_cmesh_type] < other.generator->generators[current_cmesh_type];
    }
    else {
      return current_cmesh_type < other.current_cmesh_type;
    }
  }

  generate_all_cmeshes
  get_last ()
  {
    const cmesh_types_t last = CMESH_NUM_TYPES;
    generate_all_cmeshes tmp = generate_all_cmeshes (last);
    tmp.generator->generators[last]->set_last ();
    return tmp;
  }

  generate_all_cmeshes
  operator+ (const generate_all_cmeshes &step)
  {
    if (generator->generators[current_cmesh_type]->is_at_last ()) {
      current_cmesh_type++;
    }
    else {
      generate_all_cmeshes tmp = generate_all_cmeshes (current_cmesh_type);
      generator->generators[current_cmesh_type]->get_step (tmp.generator->generators[current_cmesh_type]);
      generator->generators[current_cmesh_type]
        = *generator->generators[current_cmesh_type] + *tmp.generator->generators[current_cmesh_type];
    }
    return generate_all_cmeshes (current_cmesh_type, generator);
  }

  void
  create_cmesh ()
  {
    generator->generators[current_cmesh_type]->create_cmesh ();
  }

  t8_cmesh_t
  get_cmesh ()
  {
    return generator->generators[current_cmesh_type]->get_cmesh ();
  }

  ~generate_all_cmeshes ()
  {
  }

  int current_cmesh_type = (int) CMESH_ZERO;
  cmesh_generator_cxx *generator;
};

T8_EXTERN_C_END ();

#endif /* T8_GTEST_CMESH_GENERATOR_HXX */