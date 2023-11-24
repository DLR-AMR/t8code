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

typedef struct cmesh_generator_t
{
  sc_refcount_t rc;

  cmesh_generator *generators[CMESH_NUM_TYPES];
} cmesh_generator_cxx;

cmesh_generator_cxx *
cmesh_generator_init (void)
{
  cmesh_generator_cxx *gen = T8_ALLOC_ZERO (cmesh_generator_cxx, 1);
  t8_refcount_init (&gen->rc);

  gen->generators[CMESH_WITH_COMM] = new all_cmeshes_with_comm ();
  gen->generators[CMESH_WITH_NUM_TREES] = new all_cmeshes_with_num_elem ();
  return gen;
};

void
cmesh_generator_cxx_destroy (cmesh_generator_cxx *gen)
{
  for (int igen = CMESH_ZERO; igen < CMESH_NUM_TYPES; igen++) {
    if (gen->generators[igen] != NULL) {
      delete gen->generators[igen];
    }
  }
  T8_FREE (gen);
};

void
cmesh_generator_cxx_ref (cmesh_generator_cxx *gen)
{
  T8_ASSERT (gen != NULL);

  sc_refcount_ref (&gen->rc);
};

void
cmesh_generator_cxx_unref (cmesh_generator_cxx **pgen)
{
  cmesh_generator_cxx *gen;

  T8_ASSERT (pgen != NULL);
  gen = *pgen;
  T8_ASSERT (gen != NULL);

  if (sc_refcount_unref (&gen->rc)) {
    cmesh_generator_cxx_destroy (gen);
    *pgen = NULL;
  }
};

T8_EXTERN_C_END ();
#endif /* T8_GTEST_CMESH_GENERATOR_CXX_HXX */