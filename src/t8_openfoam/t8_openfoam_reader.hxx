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

  /** Datatype to describe the relation between a t8 and an OpenFOAM face.
   * orientation.first: t8 face vertex 0 -> OF point id
   * orientation.second: t8 face normal == OF face normal
   */
  using t8_orientation = std::pair<u_int8_t, bool>;

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
  t8_face_vertex_to_of_point (t8_eclass_t face_class, t8_orientation orientation, u_int8_t t8_face_vertex_id)
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
  build_cmesh ()
  {
    T8_ASSERT (m_cell_faces.size () != 0);
    T8_ASSERT (m_face_points.size () != 0);
    T8_ASSERT (m_points.size () != 0);

    /* Make the neighborhood vector the right size and default construct all elements. */
    m_neighborhoods.resize (m_face_points.size ());

    t8_cmesh_init (&m_cmesh);
    t8_cmesh_register_geometry<t8_geometry_linear> (m_cmesh);
    /* Reconstruct cells */
    size_t cell_id = 0;
    for (auto& cell : m_cell_faces) {
      const t8_eclass_t eclass = get_cell_eclass (cell_id);
      if (eclass == T8_ECLASS_INVALID) {
        t8_errorf ("ERROR: Encountered invalid polyhedral cell with id %li.\n", cell_id);
        return false;
      }

      /* To make this code more readable, we gather all faces with points for the current cell. */
      std::vector<size_t> current_cell_face_ids;
      current_cell_face_ids.reserve (cell.size ());
      std::vector<std::span<size_t>> current_cell_faces;
      current_cell_faces.reserve (cell.size ());
      std::vector<char> current_cell_face_normals;
      current_cell_face_normals.reserve (cell.size ());
      for (const auto& i_face : cell) {
        current_cell_face_ids.emplace_back (i_face.first);
        current_cell_faces.emplace_back (m_face_points[i_face.first]);
        current_cell_face_normals.emplace_back (i_face.second);
      }

      switch (eclass) {
      case T8_ECLASS_HEX:
        reconstruct_hex_cell (cell_id, current_cell_face_ids, current_cell_faces, current_cell_face_normals);
        break;
      default:
        SC_ABORT_NOT_REACHED ();
      }
      ++cell_id;
    }
    t8_cmesh_commit (m_cmesh, m_comm);
    return 1;
  }

  /**
   * Computes the orientation of two faces according to \ref t8_cmesh_set_join.
   * \tparam TFaceClass         The eclass of the connecting face
   * \param [in] eclass1        Cell 1 eclass
   * \param [in] face_id1       Cell 1 face id
   * \param [in] orient1        Cell 1 orientation
   * \param [in] eclass2        Cell 2 eclass
   * \param [in] face_id2       Cell 2 face id
   * \param [in] orient2        cell 2 orientation
   * \return                    The orientation in accordance with \ref t8_cmesh_set_join
   */
  int
  compute_face_orientation (t8_eclass eclass1, u_int8_t face_id1, t8_orientation orient1, t8_eclass eclass2,
                            u_int8_t face_id2, t8_orientation orient2)
  {
    const t8_eclass face_class = (t8_eclass) t8_eclass_face_types[eclass1][face_id1];
    const u_int8_t num_vertices = t8_eclass_num_vertices[face_class];

    /* First, we determine the order of cells according to the documentation of
     * \ref t8_cmesh_set_join */
    bool main_cell; /** true = cell 1, false = cell 2 */
    if (eclass1 == eclass2) {
      main_cell = face_id1 <= face_id2;
    }
    else {
      main_cell = eclass1 < eclass2;
    }

    /* Save main and other cell */
    const t8_orientation main_orient = main_cell ? orient1 : orient2;
    const t8_orientation other_orient = main_cell ? orient2 : orient1;

    /* Convert orientation */
    const int8_t steps = other_orient.first - main_orient.first;
    const int8_t direction = other_orient.second ? -1 : 1;
    const int8_t final_orientation = (num_vertices + steps * direction) % num_vertices;

    /* For quads we have to convert between OF and t8 numeration */
    if (face_class == T8_ECLASS_QUAD)
      return quad_conversion[final_orientation];
    return final_orientation;
  }

  /**
   * Compute the neighbors of a cell.
   * \tparam TNumFaces The number of faces of the cell.
   * \param [in] cell_id            The id of the cell.
   * \param [in] cell_eclass        The eclass of the cell.
   * \param [in] faces              Array containing the face ids of the cell.
   * \param [in] face_orientations  Array containing the orientation of each face (t8 vertex 0 -> OF point id)
   */
  template <size_t TNumFaces>
  void
  compute_cell_neighbors (const size_t cell_id, const t8_eclass cell_eclass,
                          const std::array<std::optional<size_t>, TNumFaces> faces,
                          const std::array<t8_orientation, TNumFaces> face_orientations)
  {
    /* We are using \ref m_neighborhoods to link neighbors together.
     * Each face has a unique id. For each face of the cell we will save the cell id, t8 face id and orientation
     * in \ref m_neighborhoods. If we encounter a face, which is already filled then we know, that the already
     * saved face and cell is the neighbor. We can now join these cells via their shared face. */
    for (u_int8_t t8_face_id = 0; t8_face_id < TNumFaces; ++t8_face_id) {
      const size_t OF_face_id = faces[t8_face_id].value ();
      if (m_neighborhoods[OF_face_id].has_value ()) {
        /* There already is a cell saved with this face. This has to be the neighbor of the current cell. */
        const size_t other_cell_id = std::get<0> (*(m_neighborhoods[OF_face_id]));
        const t8_eclass other_cell_eclass = std::get<1> (*(m_neighborhoods[OF_face_id]));
        const u_int8_t other_cell_t8_face_id = std::get<2> (*(m_neighborhoods[OF_face_id]));
        const t8_orientation other_cell_face_orientation = std::get<3> (*(m_neighborhoods[OF_face_id]));

        /* Calculate the orientation between both cells. */
        const int orientation
          = compute_face_orientation (cell_eclass, t8_face_id, face_orientations[t8_face_id], other_cell_eclass,
                                      other_cell_t8_face_id, other_cell_face_orientation);
        t8_cmesh_set_join (m_cmesh, cell_id, other_cell_id, t8_face_id, other_cell_t8_face_id, orientation);
      }
      else {
        /* There is no cell with this face saved yet. Save this cell with this face. */
        m_neighborhoods[OF_face_id].emplace (
          std::make_tuple (cell_id, cell_eclass, t8_face_id, face_orientations[t8_face_id]));
      }
    }
  }

  /**
   * Use the gathered cell information to build a hex cell and add it to the cmesh.
   * \param [in] cell_id           The id of the cell.
   * \param [in] face_ids          The ids of the cell faces.
   * \param [in] face_point_ids    The ids of the points of the cell faces.
   * \param [in] face_normals      The normals of the cell faces.
   */
  void
  reconstruct_hex_cell (size_t cell_id, std::vector<size_t> face_ids, std::vector<std::span<size_t>> face_point_ids,
                        std::vector<char> face_normals)
  {
    /* To reconstruct the Hexahedron we will follow these stages:
     * 0. Initialize some needed variables.
     * 1. Assign the first OpenFOAM face to face 0 of the t8 hexahedron.
     * 2. Search for the opposite OpenFOAM face (does not share any vertices with face 0) and make it face 1.
     * 3. Search for a third face to find out the orientation of face 1.
     * 4. We now have all 8 vertices and can use them to find and orient the last four faces.
     * 5. Use the vertices and eclass to set the cmesh cell.
     * 6. Use the face information to link the cell to its neighbors and apply the boundary conditions. */

    /* ------------------------- 0. Initialization ------------------------- */

    std::array<std::optional<size_t>, t8_eclass_num_faces[T8_ECLASS_HEX]> faces {};
    std::array<std::optional<size_t>, t8_eclass_num_vertices[T8_ECLASS_HEX]> points {};
    std::array<t8_orientation, t8_eclass_num_faces[T8_ECLASS_HEX]> face_orientations {};
    constexpr t8_eclass_t eclass = T8_ECLASS_HEX;

    /** Orientation of the OpenFOAM face with regard to the t8code face.
     * t8 corner 0 -> OF corner <orientation.first>
     * same normal direction <orientation.second> = 1, different normal direction<orientation.second> = -1
     */
    t8_orientation orientation;
    /** Normal of the t8 face. 1 points outward, 0 inward. */
    bool t8_face_normal;

    /* ------------------------- 1. Assign face 0 ------------------------- */
    t8_face_normal = t8_eclass_face_orientation[eclass][0];
    /* The orientation is 1, if both face normals point in the same direction */
    orientation.second = t8_face_normal == face_normals[0];
    orientation.first = 0;
    face_orientations[0] = orientation;

    /* Now, we assign the first face to be face 0. */
    faces[0].emplace (face_ids[0]);
    const std::span<size_t> face_0_points = face_point_ids[0];
    for (u_int8_t i_face_point = 0; i_face_point < t8_eclass_num_vertices[T8_ECLASS_QUAD]; ++i_face_point) {
      const u_int8_t of_point_id = t8_face_vertex_to_of_point (T8_ECLASS_QUAD, orientation, i_face_point);
      points[t8_face_vertex_to_tree_vertex[T8_ECLASS_HEX][0][i_face_point]].emplace (face_0_points[of_point_id]);
    }

    /* ------------------------- 2. Find face 1 ------------------------- */

    /* After assigning face 0, we have to find the opposite face, face 1.
     * It should not share any points with face 0. */
    const auto opposite_face = find_list_not_containing (face_point_ids, face_point_ids[0]);
    T8_ASSERTF (
      opposite_face.has_value (),
      "ERROR: Could not find the opposite face for hex reconstruction. Maybe this cell was falsely classified "
      "as a hex\n");
    t8_face_normal = t8_eclass_face_orientation[eclass][1];
    orientation.second = t8_face_normal == face_normals[*opposite_face];
    faces[1].emplace (face_ids[*opposite_face]);
    const auto face_1_points = face_point_ids[*opposite_face];

    /* Now we can drop face 0 and the opposite face, since we do not need it anymore.
     * We delete the opposite face first, since deleting the 0th face would shift the ids. */
    face_ids.erase (face_ids.begin () + *opposite_face);
    face_point_ids.erase (face_point_ids.begin () + *opposite_face);
    face_normals.erase (face_normals.begin () + *opposite_face);
    face_ids.erase (face_ids.begin ());
    face_point_ids.erase (face_point_ids.begin ());
    face_normals.erase (face_normals.begin ());

    /* ------------------------- 3. Compute orientation of face 1 ------------------------- */

    /* We now have to find the orientation of the face. For this we can take any OpenFOAM face which contains the 0th vertex.
     * The 1st vertex will be saved right next to it. */
    const size_t face_0_vertex_0 = points[0].value ();
    std::span<size_t>::iterator found_face_0_vertex_0;
    const auto check_face = std::find_if (
      face_point_ids.begin (), face_point_ids.end (),
      [&face_0_vertex_0, &found_face_0_vertex_0] (std::span<size_t> current_face_points) {
        found_face_0_vertex_0 = std::find (current_face_points.begin (), current_face_points.end (), face_0_vertex_0);
        return found_face_0_vertex_0 != current_face_points.end ();
      });
    T8_ASSERT (check_face != face_point_ids.end ());

    /* We found a face containing the 0th vertex. Now we can check if the next vertex is not part of face 0. If true,
     * the next vertex is vertex 1 and therefore vertex 0 of face 1. Otherwise, the previous point is vertex 0 of face 1.
     * We iterate cyclic over the face vertices. So if we arrive at the end, we begin at the start. */
    auto next_vertex
      = (found_face_0_vertex_0 + 1) == check_face->end () ? check_face->begin () : found_face_0_vertex_0 + 1;
    /* Check if face 0 contains this vertex */
    size_t face_1_vertex_0;
    if (std::find (face_1_points.begin (), face_1_points.end (), *next_vertex) != face_1_points.end ()) {
      face_1_vertex_0 = *next_vertex;
    }
    else {
      const auto previous_vertex
        = found_face_0_vertex_0 == check_face->begin () ? check_face->end () : found_face_0_vertex_0 - 1;
      face_1_vertex_0 = *previous_vertex;
    }

    /* Now we can compute the complete orientation of face 1. */
    const auto face_1_vertex_0_iterator = std::find (face_1_points.begin (), face_1_points.end (), face_1_vertex_0);
    T8_ASSERT (face_1_vertex_0_iterator != face_1_points.end ());
    orientation.first = face_1_vertex_0_iterator - face_1_points.begin ();

    /* After that we add the remaining four points, completing the vertices of the cell */
    face_orientations[1] = orientation;
    for (u_int8_t i_face_point = 0; i_face_point < t8_eclass_num_vertices[T8_ECLASS_QUAD]; ++i_face_point) {
      const u_int8_t of_point_id = t8_face_vertex_to_of_point (T8_ECLASS_QUAD, orientation, i_face_point);
      points[t8_face_vertex_to_tree_vertex[T8_ECLASS_HEX][1][i_face_point]].emplace (face_1_points[of_point_id]);
    }

    /* ------------------------- 4. Add remaining four faces ------------------------- */

    /* Lastly, we have to match the rest of the faces so we can add the boundaries to the cmesh and
     * to get the orientations and faces for neighbor linkage. */
    for (size_t i_t8_face = 2; i_t8_face < t8_eclass_num_faces[T8_ECLASS_HEX]; ++i_t8_face) {
      std::array<size_t, t8_eclass_num_vertices[T8_ECLASS_QUAD]> face_vertices;
      for (size_t i_face_vertex = 0; i_face_vertex < t8_eclass_num_vertices[T8_ECLASS_QUAD]; ++i_face_vertex) {
        face_vertices[i_face_vertex]
          = points[t8_face_vertex_to_tree_vertex[T8_ECLASS_HEX][i_t8_face][i_face_vertex]].value ();
      }
      const auto face_id = find_list_containing (face_point_ids, face_vertices);
      T8_ASSERTF (
        face_id.has_value (),
        "ERROR: Could not find a specific face for hex reconstruction. Maybe this cell was falsely classified "
        "as a hex\n");
      t8_face_normal = t8_eclass_face_orientation[eclass][1];
      orientation.second = (face_normals[*face_id] == t8_face_normal);
      const size_t vertex_0 = *points[t8_face_vertex_to_tree_vertex[T8_ECLASS_HEX][i_t8_face][0]];
      orientation.first
        = std::distance (face_point_ids[*face_id].begin (),
                         std::find (face_point_ids[*face_id].begin (), face_point_ids[*face_id].end (), vertex_0));
      face_orientations[i_t8_face] = orientation;
    }

    /* ------------------------- 5. Set cmesh cell ------------------------- */

    std::array<double, t8_eclass_num_vertices[T8_ECLASS_HEX] * 3> vertices {};
    for (size_t i_vertex = 0; i_vertex < t8_eclass_num_vertices[T8_ECLASS_HEX]; ++i_vertex) {
      for (size_t i_coord = 0; i_coord < 3; ++i_coord) {
        vertices[i_vertex * 3 + i_coord] = m_points[points[i_vertex].value ()][i_coord];
      }
    }
    t8_cmesh_set_tree_class (m_cmesh, cell_id, T8_ECLASS_HEX);
    t8_cmesh_set_tree_vertices (m_cmesh, cell_id, vertices.data (), t8_eclass_num_vertices[T8_ECLASS_HEX]);

    /* ------------------------- 6. Set neighbors and boundary conditions ------------------------- */
    compute_cell_neighbors (cell_id, T8_ECLASS_HEX, faces, face_orientations);
    // TODO: Boundary conditions
  }

  /**
   * Struct to map face indices to patch/group labels efficiently.
   *
   * Internally uses a std::map of startIndex -> hash(label).
   * Each label is valid from its startIndex up to (but not including) the start of the next label.
   */
  struct face_to_boundary_group_hash_map
  {
   public:
    /**
     * Tag for boundary group hash strong type.
     */
    struct boundary_group_hash_tag
    {
    };
    /**
     * Strong type for boundary group hashes.
     */
    using boundary_group_hash = T8Type<size_t, boundary_group_hash_tag, EqualityComparable>;

    /**
     * Add a new range starting at startIndex with a label.
     * \param[in] start_index   Index where this boundary group starts.
     * \param[in] hash          The label name (patch or group).
     */
    inline void
    add_range (const size_t start_index, const boundary_group_hash hash)
    {
      m_ranges[start_index] = hash;
    }

    static inline boundary_group_hash
    hash (const std::string_view label)
    {
      return boundary_group_hash (std::hash<std::string_view> {}(label));
    }

    /**
     * Get the boundary_group_hash for a given face index.
     * \param[in] face_index  The face index to query.
     * \return                The corresponding label, or empty string if not found.
     */
    inline std::optional<size_t>
    operator[] (size_t face_index) const
    {
      /* upper_bound gives first key > face_index */
      auto it = m_ranges.upper_bound (face_index);
      if (it == m_ranges.begin ())
        return std::nullopt; /* no label found (face before first range) */
      --it;                  /* use previous key */
      return it->second;
    }

#if T8_ENABLE_DEBUG
    /**
     * Print all ranges for debugging.
     */
    void
    print_ranges () const
    {
      t8_debugf ("OpenFOAM boundary ranges:\n");
      for (auto range = m_ranges.begin (); range != m_ranges.end (); ++range) {
        t8_debugf ("StartFace %lu -> boundary hash %lu\n", range->first, range->second.get ());
      }
    }
#endif /* T8_ENABLE_DEBUG */

   private:
    /** Map for boundary ranges: startIndex -> hash(label) */
    std::map<size_t, boundary_group_hash> m_ranges;
  };

  /** Path to the OpenFOAM case. */
  t8_path m_case_dir;
  /** The cmesh to build. */
  t8_cmesh_t m_cmesh;
  /** The assigned communicator. */
  sc_MPI_Comm m_comm;

  /** Conversion table between OpenFOAM and t8code quad numeration.
   * We only need this for quads, since triangles numerated equally.
   */
  static constexpr std::array<u_int8_t, 4> quad_conversion = { { 0, 1, 3, 2 } };

  /** Holds all points of the mesh. point_id -> x, y, z */
  std::vector<t8_3D_vec> m_points;
  /** Holds all faces of the mesh. face_id -> point_id 0, ..., point_id n */
  std::vector<std::vector<size_t>> m_face_points;
  /** Holds all cells and their faces. Second value shows of normal points outwards (1) or inwards (0).
   * cell_id -> face_id 0, ..., face_id n */
  std::vector<std::vector<std::pair<size_t, bool>>> m_cell_faces;
  /** Holds the face neighbors of the mesh and is used for joining the faces of the cmesh.
   * face id -> tree id, eclass, face of tree, orientation
   * During reconstruction, the values will be filled by the first cell owning that face.
   * If the second cell then also encounters the same face, the two cells will be linked together. */
  std::vector<std::optional<std::tuple<size_t, t8_eclass, u_int8_t, t8_orientation>>> m_neighborhoods;
  /** Assigns a hash to each boundary face. */
  face_to_boundary_group_hash_map boundary_map;
};
