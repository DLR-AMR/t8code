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

/** \file t8_subelement.hxx
 * 
 */

#pragma once
#include <t8_schemes/t8_standalone/t8_standalone_elements.hxx>
#include <t8_eclass/t8_eclass.h>

#define T8_SUBELEMENT_FACES 3
#define T8_SUB_QUAD_MAX_SUBELEMENT_TYPE 14
#define T8_SUB_QUAD_MIN_SUBELEMENT_TYPE 1
#define T8_SUB_QUAD_MAX_SUBELEMENT_ID 6
#define T8_SUB_QUAD_MIN_SUBELEMENT_ID 0

struct t8_subelement_element
{
  t8_standalone_element<T8_ECLASS_QUAD> element;
  int subelement_type
    = 0; /* saves the information, which type of transition cell a subelement is associated to (default is 0, meaning no subelement). */
  int subelement_id = 0; /* saves the information, what children subelement the given element is (default is 0) */
};
