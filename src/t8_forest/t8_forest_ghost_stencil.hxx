/*
  This file is part of t8code.
  t8code is a C library to manage a collection (a forest) of multiple
  connected adaptive space-trees of general element classes in parallel.

  Copyright (C) 2015 the developers

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

struct t8_forest_ghost_stencil : t8_forest_ghost_interface
{
  public:
  t8_forest_ghost_stencil() : t8_forest_ghost_interface(T8_GHOST_USERDEFINED){};
  
  void 
  do_ghost(t8_forest_t forest) override;

  protected:
  void add_stencil_to_ghost(t8_forest_t forest, const t8_element_t * element, t8_eclass_scheme_c *eclass_scheme, int level,
                            t8_eclass_t tree_class, t8_locidx_t ltreeid, t8_locidx_t ielement);
};

T8_EXTERN_C_END ();

#endif /* !T8_FOREST_GHOST_STENCIL_HXX */
