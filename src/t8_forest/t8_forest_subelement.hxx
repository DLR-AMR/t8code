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

/** \file t8_forest_subelement.hxx
 *  Functionality to handle subelements in a forest.
 */
#pragma once

#include <t8.h>
#include "t8_forest_general.h"
#include <t8_eclass/t8_eclass.h>
#include <t8_schemes/t8_scheme.hxx>
#include <t8_element/t8_element.h>

/** Remove hanging nodes from the forest by transitioning elements with hanging nodes into subelements. 
* \param [in] forest The input forest, which may contain hanging nodes. 
*        \a forest must be committed before calling this function. Please note that the scheme provided with the 
*         forest has to be a fitting subelement scheme.
* \return A new forest with the same number of trees and the same connectivity, but conformal without hanging nodes. 
*/
t8_forest_t
t8_forest_remove_hanging_nodes (t8_forest_t forest);

/** Remove all subelements from a forest. This is required to restore the original mesh using only recursive refinement
*    and to be able to adapt again.
* \param [in] forest The input forest which may contain subelements.
* \return A new forest with the same number of trees and the same connectivity, but without subelements. 
*/
t8_forest_t
t8_forest_discard_subelements (t8_forest_t forest);

/** Check if a forest contains subelements.
* \param [in] forest The forest to be checked.
* \return true if there are subelements in the forest, false otherwise.
*/
bool
t8_forest_has_subelements (const t8_forest_t forest);
