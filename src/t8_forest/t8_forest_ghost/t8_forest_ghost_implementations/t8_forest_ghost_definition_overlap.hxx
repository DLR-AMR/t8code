/*
  This file is part of t8code.
  t8code is a C library to manage a collection (a forest) of multiple
  connected adaptive space-trees of general element classes in parallel.

  Copyright (C) 2025 the developers

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

/** \file t8_forest_ghost_definition_overlap.hxx
 *  Implements a class of define ghost for PUMA.
 */

#ifndef T8_FOREST_GHOST_DEFINITION_OVERLAP_HXX
#define T8_FOREST_GHOST_DEFINITION_OVERLAP_HXX
 
#include <t8_forest/t8_forest_ghost/t8_forest_ghost_definition_base.hxx>
#include <t8_forest/t8_forest_types.h>
#include <vector>
#include <array>
 
/**
 * Base class for a ghost definition with elements with stretch factors.
 */
struct t8_forest_ghost_definition_overlap : public t8_forest_ghost_definition
{
  public:
  /** Base constructor with no arguments. We need this since it
   * is called from derived class constructors. */
  t8_forest_ghost_definition_overlap () : t8_forest_ghost_definition (T8_GHOST_USER_DEFINED)
  {
  }

  /** Constructor for ghost definition overlap with uniform stretch factor. */
  t8_forest_ghost_definition_overlap (std::array<double, 3> stretch_factors) : t8_forest_ghost_definition (T8_GHOST_USER_DEFINED) , _uniform_stretch_factor(stretch_factors)
  {
    _has_uniform_stretch_factor = true;
  }

  /** Get the uniform stretch factor of the ghost definition. 
   * \return uniform stretch factor in every dimension, or NULL if ther is no uniform stretch factor.
  */
  std::array<double, 3>
  get_uniform_stretch_factor() const;
  /** Set a uniform stretch factor for each dimension. */
  void
  set_uniform_stretch_factor(std::array<double, 3> stretch_factors);
  /** Get if the ghost definition uses a uniform stretch factor 
   * or computes an communicate max factors for every process. */
  bool
  has_uniform_stretch_factor() const;
  /** Change the strategy of the ghost definition from a uniform stretch factor
   * to compute and communicate max stretch factors for each process.
   */
  void
  unable_uniform_stretch_factor();

  /** Create one layer of ghost elements for a forest.
   * \param [in,out]    forest     The forest.
   * \return T8_SUBROUTINE_SUCCESS if successful, T8_SUBROUTINE_FAILURE if not.
   * \a forest must be committed before calling this function.
   */
  virtual bool
  do_ghost (t8_forest_t forest) override;

  // void
  // communicate_ownerships (t8_forest_t forest) override;

  std::array<double, 3> 
  get_uniform_stretch_factors() const;

  void
  set_uniform_stretch_factors(std::array<double, 3> stretch_factors);

  protected:
  /** Each process uses a the minimal covers from the communication 
   * to check which of its elements overlap with those of the others 
   * and creates lists of remote ghost elements for them. 
   */
  virtual void
  search_for_ghost_elements (t8_forest_t forest);

  /** Build a cover (coarsest possible grid of the local elements of a single process)
   * for each process.
   * \param [in]        forest     Committed forest with elements with uniform stretch factor
   */
  void
  build_all_cover (t8_forest_t forest);

  /** Build a cover of the local elements.
   * \param [in]        forest     committed forest.
   * \note the element scheme must support a stretch factor.
   */
  // void
  // build_owen_cover (t8_forest_t forest);

  /** Compute with build_owen_cover the cover of local process.
   * Exchange list of covers between the processes.
   * \param [in]        forest     The forest.
   */
  // void
  // communicate_covers (t8_forest_t forest);

  void
  communicate_ownerships (t8_forest_t forest) override;

  void 
  communicate_max_stretch_factor (t8_forest_t forest);

  /**
   * If memory was allocated for the offset array in communicate_ownerships it is released here.
   * Use memory_flag for this.
   */
  virtual void
  clean_up (t8_forest_t forest) override;

  bool _has_uniform_stretch_factor = false;

  std::array<double, 3> _uniform_stretch_factor;

  std::vector<std::vector<t8_element_t *>> _list_of_covers{};

  t8_shmem_array_t _max_stretch_factors = NULL;
};

#endif /* T8_FOREST_GHOST_DEFINITION_OVERLAP_HXX */