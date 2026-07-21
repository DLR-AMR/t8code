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
  t8_openfoam_reader (t8_path foamfile, sc_MPI_Comm comm)
    : m_case_dir (get_case_dir (foamfile)), m_cmesh (nullptr), m_comm (comm)
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
   * Compute the cell shape of the OpenFOAM cell
   * \param [in] cell_id  The id of the cell
   * \return              The eclass of the cell
   */
  t8_eclass_t
  get_cell_eclass (const size_t cell_id)
  {
    /* Count points of the cell. The points are counted for every face,
    so a hex for example has 6 (faces) * 4 (vertices per face) = 24 points. */
    int num_points = 0;
    for (const auto& face : m_cell_faces[cell_id]) {
      num_points += m_face_points[face.first].size ();
    }

    /* OpenFOAM only has 3D cells, so we do not check 0-2 dimensional cells. */
    switch (m_cell_faces[cell_id].size ()) {
    case 4:
      /* This cell is hopefully a tet with 4*3=12 points. */
      if (num_points == 12)
        return T8_ECLASS_TET;
      break;
    case 5:
      /* This cell could either be a prism (3 * 4 + 2 * 3 = 18) or pyramid (4 * 3 + 1 * 4 = 16). */
      if (num_points == 18)
        return T8_ECLASS_PRISM;
      if (num_points == 16)
        return T8_ECLASS_PYRAMID;
      break;
    case 6:
      /* This cell is hopefully a hex with 6 * 4 = 24 points. */
      if (num_points == 24)
        return T8_ECLASS_HEX;
      break;
    default:
      break;
    }
    return T8_ECLASS_INVALID;
  }

  /**
   * Given a face eclass and orientation converts the face id of a t8code face vertex id to the corresponding
   * OpenFOAM face id.
   * \param [in] face_class         The eclass of the face (triangle or quad)
   * \param [in] orientation        The orientation of the face (int and bool, the first is the OF vertex,
   *                                which corresponds to the t8code vertex 0 and the second is 1, if both face
   *                                normals point in the same direction and 0 otherwise)
   * \param [in] t8_face_vertex_id  The t8 face vertex id to convert
   * \return                        The converted OF face vertex id
   */
  constexpr static u_int8_t
  t8_face_vertex_to_of_point (t8_eclass_t face_class, std::pair<u_int8_t, u_int8_t> orientation,
                              u_int8_t t8_face_vertex_id)
  {
    /* While t8code uses the z-order for quads, OpenFOAM uses an anti-clockwise numeration.
     * So if the face is a quad, we have to switch indices in accordance with this LUT.
     * This is not necessary for triangles. */
    if (face_class == T8_ECLASS_QUAD) {
      t8_face_vertex_id = quad_conversion[t8_face_vertex_id];
    }

    /* The OpenFOAM face point belonging to the current t8_face_vertex is calculated via the orientation
     * (orientation.second tells us if we count up +1 or down -1 depending on the fact that the face normals point
     * in the same direction and orientation.first tells us which OF point is located at t8 vertex 0) and the current t8_face_vertex.
     * We add t8_eclass_num_vertices[face_class] once, so that we never take a mod of a negative value. */
    const int order = orientation.second ? 1 : -1;
    return (order * t8_face_vertex_id + orientation.first + t8_eclass_num_vertices[face_class])
           % t8_eclass_num_vertices[face_class];
  }

  /**
   * Finds a list inside a list which does not contain given values
   * \param [in] listlist     A list of lists (range of ranges)
   * \param [in] forbidden    A range of forbidden values
   * \return                  The location of the first sublist, which does not contain any of \a forbidden. Empty on failure.
   */
  template <typename ListList, typename Forbidden>
  static std::optional<size_t>
  find_list_not_containing (const ListList& listlist, const Forbidden& forbidden)
  {
    /* Find the first index where the corresponding list contains none of the forbidden values. */
    auto it = std::ranges::find_if (listlist, [&] (const auto& vec) {
      /* Return true if the range contains none of the forbidden values. */
      for (const auto& forbidden_val : forbidden) {
        if (std::ranges::find (vec, forbidden_val) != std::ranges::end (vec)) {
          return false;
        }
      }
      return true;
    });

    if (it != std::ranges::end (listlist)) {
      return std::distance (std::ranges::begin (listlist), it);
    }
    return std::nullopt;
  }

  /**
   * Finds a list inside a list which contains given values
   * \param [in] listlist     A list of lists (range of ranges)
   * \param [in] required     A range of required values
   * \return                  The location of the first sublist, which contains all values of \a required. Empty on failure.
   */
  template <typename ListList, typename Required>
  static std::optional<size_t>
  find_list_containing (const ListList& listlist, const Required& required)
  {
    /* Find the first index where the corresponding list contains all of the required values. */
    auto it = std::ranges::find_if (listlist, [&] (const auto& vec) {
      /* Return true if the range contains all of the required values. */
      return std::ranges::all_of (
        required, [&] (const auto& val) { return std::ranges::find (vec, val) != std::ranges::end (vec); });
    });

    if (it != std::ranges::end (listlist)) {
      return std::distance (std::ranges::begin (listlist), it);
    }
    return std::nullopt;
  }

  /**
   * Build the cmesh from the raw mesh data.
   * \return  True on success
   */
  bool
  build_cmesh ();

  /**
   * Use the gathered cell information to build a hex cell and add it to the cmesh.
   * \param [in] cell_id           The id of the cell.
   * \param [in] face_ids          The ids of the cell faces.
   * \param [in] face_point_ids    The ids of the points of the cell faces.
   * \param [in] face_normals      The normals of the cell faces.
   */
  void
  reconstruct_hex_cell (size_t cell_id, std::vector<size_t> face_ids, std::vector<std::span<size_t>> face_point_ids,
                        std::vector<char> face_normals);

  /** Path to the OpenFOAM case. */
  t8_path m_case_dir;
  /** The cmesh to build. */
  t8_cmesh_t m_cmesh;
  /** The assigned communicator. */
  sc_MPI_Comm m_comm;

  /** Conversion table between OpenFOAM and t8code quad numeration.
   * We only need this for quads, since triangles numerated equally.
   */
  static constexpr std::array<int, 4> quad_conversion = { { 0, 1, 3, 2 } };

  /** Holds all points of the mesh. point_id -> x, y, z */
  std::vector<t8_3D_vec> m_points;
  /** Holds all faces of the mesh. face_id -> point_id 0, ..., point_id n */
  std::vector<std::vector<size_t>> m_face_points;
  /** Holds all cells and their faces. Second value shows of normal points outwards (1) or inwards (0).
   * cell_id -> face_id 0, ..., face_id n */
  std::vector<std::vector<std::pair<size_t, bool>>> m_cell_faces;
};
