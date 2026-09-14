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

/** \file t8_openfoam_reader.hxx
 * Class for reading OpenFOAM meshes into t8code forests.
 */

#pragma once

#include <t8_cmesh/t8_cmesh.hxx>
#include <t8_forest/t8_forest_general.h>
#include <t8_types/t8_vec.hxx>
#include <t8_geometry/t8_geometry_implementations/t8_geometry_linear.hxx>
#include <t8_eclass/t8_eclass.h>

#include <string_view>
#include <filesystem>
#include <vector>
#include <optional>
#include <fstream>
#include <string>
#include <sstream>
#include <algorithm>
#include <span>
#include <unordered_set>

/**
 * Reads an OpenFOAM case into t8code data structures.
 */
struct t8_openfoam_reader
{
 public:
  /** Abbreviation for std::filesystem::path */
  using t8_path = std::filesystem::path;

  /**
   * Constructor of OpenFOAM reader.
   * \param [in] foamfile  Path to the *.foam file inside the OpenFOAM case directory.
   * \param [in] comm      The communicator to use for the obtained forest.
   */
  t8_openfoam_reader (t8_path foamfile, [[maybe_unused]] sc_MPI_Comm comm)
    : m_case_dir (get_case_dir (foamfile)), m_cmesh (nullptr)
  {
    SC_CHECK_ABORTF (std::filesystem::exists (foamfile), "ERROR: Foam file does not exist: %s", foamfile.c_str ());
  };

  /**
   * Reads the OpenFOAM mesh.
   * \note: When this class is finished it will return a forest instead of a cmesh.
   * This is necessary to link the mesh data to the forest and adaptation and interpolation of the data can be performed.
   * \return  A cmesh containing the OpenFOAM mesh.
   */
  t8_cmesh_t
  read ();

 private:
  /**
   * Compute the path to the OpenFOAM case (parent folder of *.foam file)
   * \param [in] foamfile   Path to the *.foam file
   * \return                Path to the OpenFOAM case
   */
  inline t8_path
  get_case_dir (t8_path foamfile)
  {
    t8_path path { foamfile };
    return path.parent_path ();
  }

  /**
   * Skips the OpenFOAM file header and positions the stream at the first numeric value.
   * Also checks if the file is in the ascii format.
   * \param[in,out] input_stream    The file stream to read from.
   * \return                        True on success, false if header not found or format not ascii.
   */
  bool
  skip_openfoam_header (std::istream& input_stream);

  /**
   * Reads an OpenFOAM label list into a vector.
   * These lists can have one of three encodings (to my knowledge):
   * expanded list:
   * \code
   * <num_values>
   * (
   * <val_0>
   * ...
   * <val_n>
   * )
   *
   * compact list:
   * <num_values>(val_0 ... val_n)
   *
   * uniform list:
   * <num_values>{val} //all values are the same
   * \endcode
   *
   * \param [in,out] input_stream   Stream positioned at the list start.
   * \param [out] values            Parsed list values.
   * \return                        True on success.
   */
  bool
  read_openfoam_label_list (std::istream& input_stream, std::vector<size_t>& values);

  /**
   * Reads in an OpenFOAM points file and stores the points in the member variable m_points.
   * \param [in] points_file The path to the points file.
   * \return                True on success.
   */
  bool
  read_points (const t8_path& points_file);

  /**
   * Reads in an OpenFOAM faces file and stores the faces in the member variable m_face_points.
   * \param [in] faces_file  The path to the faces file.
   * \return                True on success.
   */
  bool
  read_faces (const t8_path& faces_file);

  /**
   * Reads in an OpenFOAM owner file, assigns the faces to their cells and
   * also fills one half of the neighborhoods.
   * \param [in] owner_file  The path to the owner file.
   * \return                True on success.
   */
  bool
  read_owner (const t8_path& owner_file);

  /**
   * Reads in an OpenFOAM neighbor file, assigns the faces to their cells and
   * also fills the other half of the neighborhoods.
   * \param [in] neighbor_file  The path to the points file.
   * \return                    True on success.
   */
  bool
  read_neighbor (const t8_path& neighbor_file);

  /**
   * Build the cmesh from the raw mesh data.
   * \return  True on success
   */
  bool
  build_cmesh ();

  /** Path to the OpenFOAM case. */
  t8_path m_case_dir;
  /** The cmesh to build. */
  t8_cmesh_t m_cmesh;

  /** Holds all points of the mesh. point_id -> x, y, z */
  std::vector<t8_3D_vec> m_points;
  /** Holds all faces of the mesh. face_id -> point_id 0, ..., point_id n */
  std::vector<std::vector<size_t>> m_face_points;
  /** Holds all cells and their faces. Second value shows of normal points outwards (1) or inwards (0).
   * cell_id -> face_id 0, ..., face_id n */
  std::vector<std::vector<std::pair<size_t, bool>>> m_cell_faces;
};
