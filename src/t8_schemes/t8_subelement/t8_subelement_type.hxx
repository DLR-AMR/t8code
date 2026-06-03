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

/** \file t8_subelement_type.hxx
 * Definition of the element class of a subelement. A subelement always contains of an underlying element and
 * subelement type and id defining how the underlying element is transitioned into a subelement.
 */

#pragma once

/** Definition of the subelement class. A subelement always has an underlying element.
 * With the type, it is defined if the element is further defined into subelements (e.g. for hanging node resolution).
 * Type 0 means no subelement and the subelement is just the  underlying element. 
 * For hanging node resolution, the type encodes which faces are hanging and therefore the number of subelements in
 * which the element is transitioned. 
 * Accordingly, the subelement id is between 0 and num_subelement - 1.
 * \tparam TUnderlyingElement The type of the underlying element. For example a standalone element.
 */
template <typename TUnderlyingElement>
struct t8_subelement_element
{
  TUnderlyingElement element; /**< Standalone element of the subelement. */
  int subelement_type
    = 0; /**< Type of the transition cell a subelement is associated to (default is 0, meaning no subelement). */
  int subelement_id = 0; /**< Id of the children subelement the given element is (default is 0). */
};
