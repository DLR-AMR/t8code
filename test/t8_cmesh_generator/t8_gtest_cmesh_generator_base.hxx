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

#ifndef T8_GTEST_CMESH_GENERATOR_BASE_HXX
#define T8_GTEST_CMESH_GENERATOR_BASE_HXX

#include <t8_cmesh.h>

struct cmesh_generator
{
 public:
  cmesh_generator (int creator, sc_MPI_Comm comm): current_creator (creator), comm (comm), cmesh (NULL) {};

  cmesh_generator (): current_creator (0), comm (sc_MPI_COMM_WORLD), cmesh (NULL) {};

  cmesh_generator (const cmesh_generator &other)
    : current_creator (other.current_creator), comm (other.comm), cmesh (NULL) {};

  virtual bool
  operator< (const cmesh_generator &other)
    = 0;

  virtual cmesh_generator
  operator+ (const cmesh_generator &step)
    = 0;

  virtual void
  create_cmesh ()
    = 0;

  virtual void
  get_first (cmesh_generator *first)
    = 0;

  virtual void
  set_last ()
    = 0;

  virtual int
  is_at_last ()
    = 0;

  virtual void
  get_step (cmesh_generator *step)
    = 0;

  /* The cmesh is referenced and returned */
  t8_cmesh_t
  get_cmesh ()
  {
    t8_cmesh_ref (cmesh);
    return cmesh;
  }

  virtual ~cmesh_generator () {};

  int current_creator = 0;
  sc_MPI_Comm comm = sc_MPI_COMM_WORLD;
  t8_cmesh_t cmesh = NULL;
};

#endif /* T8_GTEST_CMESH_GENERATOR_BASE_HXX */