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

/** \file mesh_handle.hxx
 * Convenience header that includes the complete set of public headers of the mesh handle.
 * Include this file if you want to use the mesh handle and do not want to care about which header files to include.
 */

#pragma once

#include "concepts.hxx"
#include "mesh.hxx"
#include "element.hxx"
#include "mesh_io.hxx"
#include "constructor_wrappers.hxx"

#include "competence_pack.hxx"
#include "competences/element_data_competences.hxx"
#include "competences/dg_competences.hxx"
#include "competences/cache_element_competences.hxx"
