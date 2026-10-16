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

/** \file t8_extruded_element.hxx
 * Element type of the extruded (quasi-3D) schemes.
 * An extruded element is the Cartesian product of a 2D base element and a 1D line in the extrusion direction.
 */

#pragma once

#include <p4est.h>
#include <t8_schemes/t8_default/t8_default_line/t8_dline.h>
#include <type_traits>

/** An extruded element: a 2D base element times a line in the extrusion (z) direction.
 * Refinement only acts on the base element. The line always spans the whole tree height, i.e. it is the root line
 * TODO extruded: ??
 * (level 0, x = 0). The only exception are face neighbors across the bottom/top face constructed outside of the
 * tree, whose line is shifted by +-T8_DLINE_ROOT_LEN.
 * \note The base element must be the first member, so that a pointer to an extruded element can be passed to the
 *       functions of the base scheme.
 * \tparam TBaseElem The element type of the 2D base scheme.
 */
template <class TBaseElem>
struct t8_extruded_element
{
  TBaseElem base;  /**< The 2D base element. */
  t8_dline_t line; /**< The element in the extrusion direction. */
};

/** An extruded hex element, i.e. a quad extruded in z-direction. */
using t8_extruded_hex_element = t8_extruded_element<p4est_quadrant_t>;

// TODO extruded: ??
static_assert (std::is_standard_layout_v<t8_extruded_hex_element>,
               "An extruded element must be standard layout so that its base element is at offset zero.");
