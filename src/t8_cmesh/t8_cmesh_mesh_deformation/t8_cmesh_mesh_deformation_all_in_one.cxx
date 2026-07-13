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

/** \file t8_cmesh_mesh_deformation_all_in_one.cxx
 *  This file implements an hardcoded example for CAD-based mesh deformation.
 */

#include "t8_cmesh_mesh_deformation_all_in_one.h"
#if T8_ENABLE_OCC
#include <t8_cad/t8_cad_handle.hxx>
#include <t8_cmesh/t8_cmesh_mesh_deformation/t8_cmesh_mesh_deformation.hxx>
#include <t8_cmesh/t8_cmesh_mesh_deformation/t8_rbf.hxx>
#endif

#include <filesystem>
#include <iostream>
#include <vector>

namespace fs = std::filesystem;

static std::vector<fs::path>
findBrepFiles (const char* folder)
{
  std::vector<fs::path> files;

  for (const auto& entry : fs::directory_iterator (folder)) {
    if (entry.is_regular_file () && entry.path ().extension () == ".brep") {
      files.push_back (entry.path ());
    }
  }

  return files;
}

T8_EXTERN_C_BEGIN ();

int
t8_forest_hacky_deformation ([[maybe_unused]] t8_forest_t forest, [[maybe_unused]] const char* brep_folder)
{
  static size_t ifile = 0;
  t8_productionf ("ifile = %li\n", ifile);

  auto files = findBrepFiles (brep_folder);

  t8_productionf ("num_files = %li\n", files.size ());
  if (ifile == files.size ()) {
    return 0;
  }
  t8_cmesh_t cmesh = t8_forest_get_cmesh (forest);

  /* Load CAD geometry from .brep file. */
  auto cad = std::make_shared<t8_cad_handle> (std::filesystem::path (files[ifile]).replace_extension ().c_str ());

  /* Initialize the deformation object for the given mesh. */
  t8_cmesh_mesh_deformation deformation (cmesh);

  /** Save the input RBF type. T8_RBF_CP_C2 or T8_RBF_TPS. */
  t8_rbf_function_type rbf_type = T8_RBF_TPS;

  auto displacements = deformation.calculate_displacement_surface_vertices (cad.get (), rbf_type, 4);

  deformation.apply_vertex_displacements (displacements, cad, rbf_type);

  ++ifile;
  return 1;
}

T8_EXTERN_C_END ();
