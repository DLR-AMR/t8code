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

/** \file t8_openfoam_reader.cxx
 * Implementation details for \ref t8_openfoam_reader.hxx
 */

#include <t8_openfoam/t8_openfoam_reader.hxx>

t8_cmesh_t
t8_openfoam_reader::read ()
{

  /* Build paths to the standard OpenFOAM mesh files. */
  const t8_path case_faces_file = m_case_dir / "constant/polyMesh/faces";
  const t8_path case_boundary_file = m_case_dir / "constant/polyMesh/boundary";
  const t8_path case_neighbor_file = m_case_dir / "constant/polyMesh/neighbour";
  const t8_path case_owner_file = m_case_dir / "constant/polyMesh/owner";
  const t8_path case_points_file = m_case_dir / "constant/polyMesh/points";

  /* The file reading needs to happen in this order, since these functions depend on each other. */
  bool error = 0;
  error = !read_points (case_points_file);
  error = error && !read_faces (case_faces_file);
  error = error && !read_owner (case_owner_file);
  error = error && !read_neighbor (case_neighbor_file);

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

bool
t8_openfoam_reader::skip_openfoam_header (std::istream& input_stream)
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

  /* Loop through file till end of file or till we break manually at the right position. */
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

bool
t8_openfoam_reader::read_openfoam_label_list (std::istream& input_stream, std::vector<size_t>& values)
{
  /* Get size of list */
  size_t list_size = 0;
  input_stream >> list_size;

  /* Get list opening delimiter */
  char opening_delimiter = '\0';
  input_stream >> opening_delimiter;

  values.clear ();
  values.reserve (list_size);

  switch (opening_delimiter) {
  /* Expanded or compact list: <num_values>(val_0 ... val_n) with or without line break */
  case '(': {
    /* Read values */
    while (values.size () < list_size) {
      size_t list_value;
      input_stream >> list_value;
      values.emplace_back (list_value);
    }

    /* Check if closing delimiter is at the expected place. */
    char closing_delimiter = '\0';
    input_stream >> closing_delimiter;
    if (closing_delimiter == ')') {
      return true;
    }
    [[fallthrough]];
  }

  /* Uniform encoding: <num_values>{val} */
  case '{': {
    /* Get the uniform value. */
    size_t uniform_value = 0;
    input_stream >> uniform_value;

    /* Assign uniform values to list. */
    values.assign (list_size, uniform_value);
    return true;

    /* Check for expected closing delimiter. */
    char closing_delimiter = '\0';
    input_stream >> closing_delimiter;
    if (closing_delimiter == '}') {
      return true;
    }
    [[fallthrough]];
  }

  default:
    t8_errorf ("Unrecognized encoding of OpenFOAM list.\n");
    return false;
  }
}

bool
t8_openfoam_reader::read_points (const t8_path& points_file)
{
  /* TODO: Implement compact lists for points (mostly happen when mesh only has one element) */

  std::ifstream file { points_file };
  if (!file) {
    t8_errorf ("ERROR: File not found: %s\n", points_file.c_str ());
    m_points.clear ();
    return false;
  }
  /* Skip OpenFOAM header. */
  if (!skip_openfoam_header (file)) {
    t8_errorf ("ERROR: Unable to parse file: %s\n", points_file.c_str ());
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
      t8_errorf ("ERROR: Unable to parse file: %s\n", points_file.c_str ());
      m_points.clear ();
      return false;
    }
    m_points.emplace_back (point);
  }
  return true;
}

bool
t8_openfoam_reader::read_faces (const t8_path& faces_file)
{
  std::ifstream file (faces_file);
  /* Check if file exists. */
  if (!file) {
    t8_errorf ("ERROR: File not found: %s\n", faces_file.c_str ());
    m_face_points.clear ();
    return false;
  }
  /* Skip OpenFOAM header. */
  if (!skip_openfoam_header (file)) {
    t8_errorf ("ERROR: Unable to parse file: %s\n", faces_file.c_str ());
    return false;
  }

  /* Get number of faces */
  size_t n_faces = 0;
  std::string line;
  std::getline (file, line);
  std::istringstream iss (line);
  if (!(iss >> n_faces)) {
    t8_errorf ("ERROR: Unable to parse file: %s\n", faces_file.c_str ());
  }

  /* skip opening '(' */
  std::getline (file, line);

  m_face_points.reserve (n_faces);

  /* Read all faces into m_face_points */
  for (size_t i_face = 0; i_face < n_faces; ++i_face) {
    std::vector<size_t> points;
    if (!read_openfoam_label_list (file, points)) {
      t8_errorf ("ERROR: Unable to parse file: %s\n", faces_file.c_str ());
      return false;
    }
    m_face_points.emplace_back (std::move (points));
  }
  return true;
}

bool
t8_openfoam_reader::read_owner (const t8_path& owner_file)
{
  T8_ASSERT (!m_face_points.empty ());

  std::ifstream file (owner_file);
  /* Check if file can be opened */
  if (!file) {
    t8_errorf ("ERROR: File not found: %s\n", owner_file.c_str ());
    return false;
  }
  if (!skip_openfoam_header (file)) {
    t8_errorf ("ERROR: Unable to parse file: %s\n", owner_file.c_str ());
    return false;
  }

  /** List of read values from owners file. Position is equal to face id and value
     * is the owner of the face. face_id -> cell_id */
  std::vector<size_t> owner_list;
  if (!read_openfoam_label_list (file, owner_list)) {
    t8_errorf ("ERROR: Unable to parse file: %s\n", owner_file.c_str ());
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

bool
t8_openfoam_reader::read_neighbor (const t8_path& neighbor_file)
{
  T8_ASSERT (!m_cell_faces.empty ());
  /* Check if file can be opened */
  std::ifstream file (neighbor_file);
  if (!file) {
    t8_errorf ("ERROR: File not found: %s\n", neighbor_file.c_str ());
    return false;
  }

  if (!skip_openfoam_header (file)) {
    t8_errorf ("ERROR: Unable to parse file: %s\n", neighbor_file.c_str ());
    return false;
  }

  /** List of read values from neighbor file. Position is equal to face id and value
     * is the neighbor of the face. face_id -> cell_id */
  std::vector<size_t> neighbor_list;
  if (!read_openfoam_label_list (file, neighbor_list)) {
    t8_errorf ("ERROR: Unable to parse file: %s\n", neighbor_file.c_str ());
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

bool
t8_openfoam_reader::build_cmesh ()
{
  T8_ASSERT (!m_cell_faces.empty ());
  T8_ASSERT (!m_face_points.empty ());
  T8_ASSERT (!m_points.empty ());

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

void
t8_openfoam_reader::reconstruct_hex_cell (size_t cell_id, std::vector<size_t> face_ids,
                                          std::vector<std::span<size_t>> face_point_ids, std::vector<char> face_normals)
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
  std::array<int8_t, t8_eclass_num_faces[T8_ECLASS_HEX]> face_orientations {};
  constexpr t8_eclass_t eclass = T8_ECLASS_HEX;

  /** Orientation of the OpenFOAM face with regard to the t8code face.
     * t8 corner 0 -> OF corner <orientation.first>
     * same normal direction <orientation.second> = 1, different normal direction<orientation.second> = -1
     */
  std::pair<u_int8_t, bool> orientation;
  /** Normal of the t8 face. 1 points outward, 0 inward. */
  bool t8_face_normal;

  /* ------------------------- 1. Assign face 0 ------------------------- */
  t8_face_normal = t8_eclass_face_orientation[eclass][0];
  /* The orientation is 1, if both face normals point in the same direction */
  orientation.second = t8_face_normal == face_normals[0];
  orientation.first = 0;

  /* Now, we assign the first face to be face 0. */
  faces[0].emplace (face_ids[0]);
  face_orientations[0] = 0;
  const std::span<size_t> face_0_points = face_point_ids[0];
  for (u_int8_t i_face_point = 0; i_face_point < t8_eclass_num_vertices[T8_ECLASS_QUAD]; ++i_face_point) {
    const u_int8_t of_point_id = t8_face_vertex_to_of_point (T8_ECLASS_QUAD, orientation, i_face_point);
    points[t8_face_vertex_to_tree_vertex[T8_ECLASS_HEX][0][i_face_point]].emplace (face_0_points[of_point_id]);
  }

  /* ------------------------- 2. Find face 1 ------------------------- */

  /* After assigning face 0, we have to find the opposite face, face 1.
     * It should not share any points with face 0. */
  const auto opposite_face = find_list_not_containing (face_point_ids, face_point_ids[0]);
  T8_ASSERTF (opposite_face.has_value (),
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
  face_orientations[1] = orientation.first;
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
    T8_ASSERTF (face_id.has_value (),
                "ERROR: Could not find a specific face for hex reconstruction. Maybe this cell was falsely classified "
                "as a hex\n");
    //t8_face_normal = t8_eclass_face_orientation[eclass][1];
    //orientation.second = (face_normals[*face_id] == t8_face_normal);
    const size_t vertex_0 = *points[t8_face_vertex_to_tree_vertex[T8_ECLASS_HEX][i_t8_face][0]];
    face_orientations[i_t8_face]
      = std::distance (face_point_ids[*face_id].begin (),
                       std::find (face_point_ids[*face_id].begin (), face_point_ids[*face_id].end (), vertex_0));
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
  //TODO
}
