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

#ifndef T8_GTEST_CMESH_CREATOR_BASE_HXX
#define T8_GTEST_CMESH_CREATOR_BASE_HXX

#include <t8_cmesh.h>

/**
 * A base class for cmesh_creators. 
 * 
 */
class cmesh_creator {
 public:
  /**
   *  Construct a new cmesh generator object
   *  Is set to the first creator, using sc_MPI_COMM_WORLD
   */
  cmesh_creator (): current_creator (0), comm (sc_MPI_COMM_WORLD), cmesh (NULL) {};

  /**
   * Construct a new cmesh generator object, with a creator and comm specified
   * 
   * \param[in] creator The creator to use
   * \param[in] comm The communicator to use
   */
  cmesh_creator (int creator, sc_MPI_Comm comm): current_creator (creator), comm (comm), cmesh (NULL) {};

  /**
   *  Copy-constructor of a new cmesh generator object
   * 
   * \param other[in] The cmesh_creator to copy
   */
  cmesh_creator (const cmesh_creator &other)
    : current_creator (other.current_creator), num_trees (other.num_trees), comm (other.comm)
  {
    if (cmesh != NULL) {
      t8_cmesh_ref (other.cmesh);
      cmesh = other.cmesh;
    }
    else {
      cmesh = NULL;
    }
  };

  /**
   * Compare two cmesh_creator
   * 
   * \param[in] other the cmesh_creator to compare with
   * \return true if both cmesh_creators are equal
   * \return false if the cmesh_creators differ
   */
  virtual bool
  operator< (const cmesh_creator &other)
    = 0;

  /**
   * Function to increase to cmesh_creator to the next creator function
   * 
   * \param step A cmesh_creator describing by how much to increase 
   */
  virtual void
  addition (const std::shared_ptr<cmesh_creator> step)
    = 0;

  /**
   * Create the cmesh. Initially no cmesh is created. Calling this function will create the cmesh
   */
  virtual void
  create_cmesh ()
    = 0;

  /**
   * Set the creator to the first creator function (the function that creates a cmesh) with
   * the first set of parameters. 
   */
  virtual void
  set_first ()
    = 0;

  /**
   * Set the creator to the last creator function (the function that creates a cmesh) with
   * the last set of parameters. 
   */
  virtual void
  set_last ()
    = 0;

  /**
   * Check if the current set of parameters is the last set of parameters to use. 
   * 
   * \return true if the current setup is the last
   * \return false ow
   */
  virtual bool
  is_at_last ()
    = 0;

  /**
   * Get the cmesh object, if it has been created
   * 
   * \return the cmesh, if it has been created, NULL if not. 
   */
  t8_cmesh_t
  get_cmesh ()
  {
    if (cmesh != NULL) {
      t8_cmesh_ref (cmesh);
      return cmesh;
    }
    return NULL;
  }

  /**
   * Dereference the current cmesh. Will be deleted if the last reference is unset. 
   * 
   */
  void
  unref_cmesh ()
  {
    if (cmesh != NULL) {
      t8_cmesh_unref (&cmesh);
    }
  }

  /**
   * Destroy the cmesh generator object
   * 
   */
  virtual ~cmesh_creator () {};

  /* Defines the current creator function */
  int current_creator = 0;
  /* If needed, the number of trees to create by the creator function*/
  int num_trees = 1;
  /* The communicator to use. */
  sc_MPI_Comm comm = sc_MPI_COMM_WORLD;
  /* The cmesh. */
  t8_cmesh_t cmesh = NULL;
};

#endif /* T8_GTEST_CMESH_CREATOR_BASE_HXX */