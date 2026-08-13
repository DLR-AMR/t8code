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

#ifndef T8_FOREST_GHOST_DEFINITION_OVERLAP_FOR_TESTING_HXX
#define T8_FOREST_GHOST_DEFINITION_OVERLAP_FOR_TESTING_HXX

#include <t8_forest/t8_forest_ghost/t8_forest_ghost_implementations/t8_forest_ghost_definition_overlap.hxx>
#include <t8_forest/t8_forest_types.h>
#include <array>

/* For writing the cover in the vtk use a struct with data. */
struct data_struct_if_element_has_cover
{
  /** every element get a cover number. 
   *    -1 = element has no cover
   *     n = element is covered by the n-th element of the cover.
  */
  int covernumber;
};

/**
 * Test class for a ghost definition with elements with stretch factors.
 * The class exists to run tests on the ghost definition overlap.
 */
struct t8_forest_ghost_definition_overlap_for_testing: public t8_forest_ghost_definition_overlap
{
 public:
  /** Base constructor with no arguments. We need this since it
   * is called from derived class constructors. */
  t8_forest_ghost_definition_overlap_for_testing (): t8_forest_ghost_definition_overlap ()
  {
  }

  /** Constructor for ghost definition overlap with uniform stretch factor. */
  t8_forest_ghost_definition_overlap_for_testing (std::array<double, 3> stretch_factors)
    : t8_forest_ghost_definition_overlap (stretch_factors)
  {
  }

  void
  set_write_forest_as_vtk (bool write_as_vtk)
  {
    _write_forest_as_vtk = write_as_vtk;
  }

  /**
   *  Check if every local leaf element have on cover element as ancestor.
   * \param [in]      forest    The forest on which the check should be computed.
   * \return a bool if the cover has the desired property.
   */
  bool
  check_cover (t8_forest_t forest);

  /**
   * \param [in]    forest    The forest to print.
   * \param [in]    cover     A sorted list of cover elements.
   * \note It is assumed that the elements in the cover are sorted in ascending order by ther linear id.
   */
  void
  write_forest_with_cover_as_vtk (t8_forest_t forest, std::vector<t8_element_t *> cover);

 protected:
  struct data_struct_if_element_has_cover *
  init_data_if_element_has_cover (t8_forest_t forest, const t8_scheme_c *eclass_scheme, const t8_eclass_t tree_class,
                                  std::vector<t8_element_t *> &cover);

  bool _write_forest_as_vtk { false };
};

#endif /* T8_FOREST_GHOST_DEFINITION_OVERLAP_FOR_TESTING_HXX */
