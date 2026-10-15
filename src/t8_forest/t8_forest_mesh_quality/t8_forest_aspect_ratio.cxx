/*
  This file is part of t8code.
  t8code is a C library to manage a collection (a forest) of multiple
  connected adaptive space-trees of general element classes in parallel.

  Copyright (C) 2026 the developers

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

/** \file TODO
 */

#include <t8.h>
#include <t8_forest/t8_forest_general.h>
#include <t8_eclass/t8_eclass.h>
#include "t8_forest_aspect_ratio.hxx"
#include <t8_schemes/t8_scheme.hxx>
#include <t8_forest/t8_forest_geometrical.h> /* This is only needed for the dummy call of the volume function!! */

/* Aspect ratio helper that works only for quads. */
double
t8_forest_element_aspect_ratio_quad (t8_forest_t forest, t8_locidx_t ltreeid, const t8_element_t *element)
{  // Dummy call of the element volume just to do "something" with the inputs.
  return t8_forest_element_volume (forest, ltreeid, element);
}

double
t8_forest_element_aspect_ratio (t8_forest_t forest, t8_locidx_t ltreeid, const t8_element_t *element)
{
  // I think its a good idea to have one function for each element type. Maybe you dont need this but i implemented it this way for now.
  // You can have a look at t8_eclass.h for all the shapes that we support. Its totally fine if you start only with quads and extend this step by step.
  if (t8_forest_get_scheme (forest)->element_get_shape (t8_forest_get_tree_class (forest, ltreeid), element)
      == T8_ECLASS_QUAD) {
    return t8_forest_element_aspect_ratio_quad (forest, ltreeid, element);
  }
  return 0.;
}

double
t8_forest_average_aspect_ratio ([[maybe_unused]] t8_forest_t forest)
{
  return 0.;
}
