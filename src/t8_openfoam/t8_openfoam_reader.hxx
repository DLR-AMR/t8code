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
   * Reader for OpenFOAM cases.
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
  read ()
  {

    /* Build paths to the standard OpenFOAM mesh files. */
    const t8_path case_faces_dir = m_case_dir / "constant/polyMesh/faces";
    const t8_path case_boundary_dir = m_case_dir / "constant/polyMesh/boundary";
    const t8_path case_neighbor_dir = m_case_dir / "constant/polyMesh/neighbour";
    const t8_path case_owner_dir = m_case_dir / "constant/polyMesh/owner";
    const t8_path case_points_dir = m_case_dir / "constant/polyMesh/points";

    /* The file reading needs to happen in this order, since these functions depend on each other. */
    bool error = 0;
    error = !read_points (case_points_dir);
    error &= !read_faces (case_faces_dir);
    error &= !read_owner (case_owner_dir);
    error &= !read_neighbor (case_neighbor_dir);

    if (error) {
      /* Return the uninitialized cmesh (nullptr) */
      t8_errorf ("ERROR during OpenFOAM case reading.\n");
      return nullptr;
    }
    if (!build_cmesh ()) {
      t8_errorf ("ERROR during OpenFOAM cmesh build.\n");
      return nullptr;
    }
    return m_cmesh;
  }

 private:
  /**
   * Compute the path to the OpenFOAM case (parent folder of *.foam file)
   * \param [in] foamfile   Path to the *.foam file
   * \return                Path to the OpenFOAM case
   */
  t8_path
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
  skip_openfoam_header (std::istream& input_stream)
  {

    /* Typical OpenFOAM header example */

    /*--------------------------------*- C++ -*----------------------------------*\
    //    =========                 |
    //    \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
    //     \\    /   O peration     | Website:  https://openfoam.org
    //      \\  /    A and           | Version:  13
    //       \\/     M anipulation  |
    //  \*---------------------------------------------------------------------------*/
    //  FoamFile
    //  {
    //      format      ascii;                            <-- format string we are searching for
    //      class       vectorField;
    //      location    "constant/polyMesh";
    //      object      points;
    //  }
    //  // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
    //
    //
    //  18        <-- first numeric value we are searching for
    //  (
    //  (0 0 0)
    //  (0.5 0 0)
    //  (1 0 0)
    //  (0 0.5 0)
    //  [ ... ]

    std::string line;
    bool format_checked = false;

    /* Loop through file till end of file or thill we break manually at the right position. */
    while (std::getline (input_stream, line)) {
      std::istringstream line_stream (line);

      /* If the format was not checked yet try to check it now. */
      if (!format_checked) {
        /* Check for ascii format */
        std::string field;
        line_stream >> field;

        /* If the line contains the word "format" it is the right line */
        if (field == "format") {
          std::string format_value;
          line_stream >> format_value;

          /* Remove trailing semicolon if present */
          if (!format_value.empty () && format_value.back () == ';')
            format_value.pop_back ();

          if (format_value != "ascii") {
            t8_errorf ("ERROR: OpenFOAM file format is not ascii: %s\n", format_value.c_str ());
            return false;
          }
          format_checked = true;
          continue;
        }
      }

      /* We do not need to check for numeric values when the format string did not appear yet. */
      if (format_checked) {
        /* Skip to first numeric value and then revert to start of numeric value */
        size_t numeric_value;
        if (line_stream >> numeric_value) {
          /* Numeric value found, move stream to start of line again */
          input_stream.seekg (-static_cast<std::streamoff> (line.size ()) - 1, std::ios_base::cur);
          return true;
        }
      }
    }

    /* Reaching this code means something went wrong. */
    if (!format_checked) {
      t8_errorf ("ERROR: Could not find format specification for OpenFOAM file.\n");
      return false;
    }
    t8_errorf ("ERROR: Could not find numeric list size after header.\n");
    return false;
  }

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
  read_openfoam_label_list (std::istream& input_stream, std::vector<size_t>& values)
  {
    /* Get size of list */
    size_t list_size = 0;
    input_stream >> list_size;

    /* Get list opening delimiter */
    char opening_delimiter = '\0';
    input_stream >> opening_delimiter;

    values.clear ();
    values.reserve (list_size);

    /* Expanded or compact list: <num_values>(val_0 ... val_n) with or without line break */
    if (opening_delimiter == '(') {
      /* Read values */
      while (values.size () < list_size) {
        size_t list_value;
        input_stream >> list_value;
        values.emplace_back (list_value);
      }

      /* Check if closing delimiter is at the expected place. */
      char closing_delimiter = '\0';
      input_stream >> closing_delimiter;
      if (closing_delimiter != ')') {
        return false;
      }

      return true;
    }

    /* Uniform encoding: <num_values>{val} */
    if (opening_delimiter == '{') {
      /* Get the uniform value. */
      size_t uniform_value = 0;
      input_stream >> uniform_value;

      /* Check for expected closing delimiter. */
      char closing_delimiter = '\0';
      input_stream >> closing_delimiter;
      if (closing_delimiter != '}') {
        return false;
      }

      /* Assign uniform values to list. */
      values.assign (list_size, uniform_value);
      return true;
    }

    return false;
  }

  /**
   * Reads in an OpenFOAM points file and stores the points.
   * \param [in] points_dir The path to the points file.
   * \return                True on success.
   */
  bool
  read_points (const t8_path& points_dir)
  {
    /* TODO: Implement compact lists for points (mostly happen when mesh only has one element) */

    std::ifstream file { points_dir };
    if (!file) {
      t8_errorf ("ERROR: File not found: %s\n", points_dir.c_str ());
      m_points.clear ();
      return false;
    }
    /* Skip OpenFOAM header. */
    if (!skip_openfoam_header (file)) {
      t8_errorf ("ERROR: Unable to parse file: %s\n", points_dir.c_str ());
      return false;
    }

    /* Get the point count */
    std::string line;
    size_t n_points = 0;
    while (std::getline (file, line)) {
      std::istringstream iss (line);
      if (iss >> n_points)
        break;
    }
    T8_ASSERTF (n_points > 0, "ERROR: Mesh contains no coordinates.");

    std::vector<std::array<double, 3>> points;
    points.reserve (n_points);

    /* Skip opening '(' */
    std::getline (file, line);

    /* Read all points */
    for (size_t i_point = 0; i_point < n_points; ++i_point) {
      std::getline (file, line);
      std::istringstream iss (line);
      char c;
      t8_3D_vec point;
      if (!(iss >> c >> point[0] >> point[1] >> point[2] >> c)) {
        t8_errorf ("ERROR: Unable to parse file: %s\n", points_dir.c_str ());
        m_points.clear ();
        return false;
      }
      m_points.emplace_back (point);
    }
    return true;
  }

  /**
   * Reads in an OpenFOAM faces file stores the faces.
   * \param [in] faces_dir  The path to the faces file.
   * \return                True on success.
   */
  bool
  read_faces (const t8_path& faces_dir)
  {
    std::ifstream file (faces_dir);
    /* Check if file exists. */
    if (!file) {
      t8_errorf ("ERROR: File not found: %s\n", faces_dir.c_str ());
      m_face_points.clear ();
      return false;
    }
    /* Skip OpenFOAM header. */
    if (!skip_openfoam_header (file)) {
      t8_errorf ("ERROR: Unable to parse file: %s\n", faces_dir.c_str ());
      return false;
    }

    /* Get number of faces */
    size_t n_faces = 0;
    std::string line;
    std::getline (file, line);
    std::istringstream iss (line);
    if (!(iss >> n_faces)) {
      t8_errorf ("ERROR: Unable to parse file: %s\n", faces_dir.c_str ());
    }

    /* skip opening '(' */
    std::getline (file, line);

    m_face_points.reserve (n_faces);

    /* Read all faces into m_face_points */
    for (size_t i_face = 0; i_face < n_faces; ++i_face) {
      std::vector<size_t> points;
      if (!read_openfoam_label_list (file, points)) {
        t8_errorf ("ERROR: Unable to parse file: %s\n", faces_dir.c_str ());
        return false;
      }
      m_face_points.emplace_back (std::move (points));
    }
    return true;
  }

  /**
   * Reads in an OpenFOAM owner file, assigns the faces to their cells and
   * also fills one half of the neighborhoods.
   * \param [in] owner_dir  The path to the owner file.
   * \return                True on success.
   */
  bool
  read_owner (const t8_path& owner_dir)
  {
    T8_ASSERT (!m_face_points.empty ());

    std::ifstream file (owner_dir);
    /* Check if file can be opened */
    if (!file) {
      t8_errorf ("ERROR: File not found: %s\n", owner_dir.c_str ());
      return false;
    }
    if (!skip_openfoam_header (file)) {
      t8_errorf ("ERROR: Unable to parse file: %s\n", owner_dir.c_str ());
      return false;
    }

    /** List of read values from owners file. Position is equal to face id and value
     * is the owner of the face. face_id -> cell_id */
    std::vector<size_t> owner_list;
    if (!read_openfoam_label_list (file, owner_list)) {
      t8_errorf ("ERROR: Unable to parse file: %s\n", owner_dir.c_str ());
      return false;
    }

    /* The number of cells in the mesh is the highest owner_id + 1 (since it starts at 0). */
    const size_t num_cells = *std::max_element (owner_list.begin (), owner_list.end ()) + 1;
    t8_debugf ("Reading OpenFOAM mesh with %li cells.\n", num_cells);
    m_cell_faces.resize (num_cells);

    /* Assign faces to their cells and fill neighborhoods. */
    for (size_t face_id = 0; face_id < owner_list.size (); ++face_id) {
      const size_t cell_id = owner_list[face_id];
      /* Add face to owner cell. Normal points outwards (1). */
      m_cell_faces[cell_id].emplace_back (std::make_pair (face_id, 1));
    }

    return true;
  }

  /**
   * Reads in an OpenFOAM owner file, assigns the faces to their cells and
   * also fills the other half of the neighborhoods.
   * \param [in] neighbor_dir  The path to the points file.
   * \return                    True on success.
   */
  bool
  read_neighbor (const t8_path& neighbor_dir)
  {
    T8_ASSERT (!m_cell_faces.empty ());
    /* Check if file can be opened */
    std::ifstream file (neighbor_dir);
    if (!file) {
      t8_errorf ("ERROR: File not found: %s\n", neighbor_dir.c_str ());
      return false;
    }

    if (!skip_openfoam_header (file)) {
      t8_errorf ("ERROR: Unable to parse file: %s\n", neighbor_dir.c_str ());
      return false;
    }

    /** List of read values from neighbor file. Position is equal to face id and value
     * is the neighbor of the face. face_id -> cell_id */
    std::vector<size_t> neighbor_list;
    if (!read_openfoam_label_list (file, neighbor_list)) {
      t8_errorf ("ERROR: Unable to parse file: %s\n", neighbor_dir.c_str ());
      return false;
    }

    /* Assign faces to their cells and fill neighborhoods. */
    for (size_t face_id = 0; face_id < neighbor_list.size (); ++face_id) {
      const size_t cell_id = neighbor_list[face_id];
      /* Add face to owner cell. Normal points inwards (-1). */
      m_cell_faces[cell_id].emplace_back (std::make_pair (face_id, 0));
    }

    return true;
  }

  /**
   * Build the cmesh from the raw mesh data.
   * \return  True on success
   */
  bool
  build_cmesh ()
  {
    t8_global_errorf ("ERROR: Not implemented yet. \n");
    return 0;
  }

  /** Path to the OpenFOAM case. */
  t8_path m_case_dir;
  /** The cmesh to build. */
  t8_cmesh_t m_cmesh;
  /** The assigned communicator. */
  sc_MPI_Comm m_comm;

  /** Holds all points of the mesh. point_id -> x, y, z */
  std::vector<t8_3D_vec> m_points;
  /** Holds all faces of the mesh. face_id -> point_id 0, ..., point_id n */
  std::vector<std::vector<size_t>> m_face_points;
  /** Holds all cells and their faces. Second value shows of normal points outwards (1) or inwards (0).
   * cell_id -> face_id 0, ..., face_id n */
  std::vector<std::vector<std::pair<size_t, bool>>> m_cell_faces;
};
