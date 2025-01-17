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

#ifndef T8_FOREST_GHOST_STENCIL_HXX
#define T8_FOREST_GHOST_STENCIL_HXX

#include <t8_forest/t8_forest_ghost_interface.h>
#include <t8_forest/t8_forest_ghost_interface.hxx>

T8_EXTERN_C_BEGIN ();

struct t8_forest_ghost_stencil: t8_forest_ghost_definition
{
 public:
  /**
   * Only constructor of t8_forest_ghost_stencil.
   * Derived class of t8_forest_ghost_definition.
   * Hase userdefined as type.
   */
  t8_forest_ghost_stencil (): t8_forest_ghost_definition (T8_GHOST_USERDEFINED) {};

  /** Create one layer of ghost elements for a forest.
   * The neighborhood for this is defined by an stencil.
   * \see add_stencil_to_ghost
   * \param [in,out]    forest     The forest.
   * \a forest must be committed before calling this function.
   */
  void
  do_ghost (t8_forest_t forest) override;

 protected:
  /**
   * Add this stencil (elements N and F) for element E to ghost.
   * 
   *           N
   *           |
   *       N - F - N
   *       |   |   |
   *   N - F - E - F - N
   *       |   |   |
   *       N - F - N
   *           |
   *           N
   * \param[in] forest          a commit uniform forest
   * \param[in] element         element E
   * \param[in] eclass_scheme   
   * \param[in] level           level of the uniform forest
   * \param[in] tree_class      expect T8_ECLASS_QUAD
   * \param[in] ltreeid         local tree id, expect 0
   * \param[in] ielement        local index of element E
   * \note some parameters currently expect specific values
   */
  void
  add_stencil_to_ghost (t8_forest_t forest, const t8_element_t *element, const t8_eclass_scheme_c *eclass_scheme,
                        const int level, const t8_eclass_t tree_class, const t8_locidx_t ltreeid,
                        const t8_locidx_t ielement);
};

T8_EXTERN_C_END ();

#endif /* !T8_FOREST_GHOST_STENCIL_HXX */
