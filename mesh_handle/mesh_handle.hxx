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

/** \file mesh_handle.hxx
 * Convenience header that includes the complete set of public headers of the mesh handle.
 * Include this file if you want to use the mesh handle and do not want to care about which header files to include.
 */

#pragma once

#include "concepts.hxx"             /* Concepts to constraint template parameters related to the mesh handle. */
#include "mesh.hxx"                 /* The mesh class. */
#include "element.hxx"              /* Class of the elements of the mesh. */
#include "mesh_io.hxx"              /* In- and output of meshes. */
#include "constructor_wrappers.hxx" /* Wrapper to construct a mesh handle instance from a cmesh instead of forest. */

#include "competence_pack.hxx" /* Mesh and element competence packs to extend the functionality of the mesh. */
#include "competences/element_data_competences.hxx" /* Competences to use element data with the mesh. */
#include "competences/dg_competences.hxx"           /* Competences useful for discontinuous Galerkin methods. */
#include "competences/cache_element_competences.hxx" /* Competences to cache element properties instead of recalculation. */
