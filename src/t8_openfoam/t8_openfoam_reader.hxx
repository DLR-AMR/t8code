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
#include <t8_vector/t8_vec.hxx>
#include <t8_geometry/t8_geometry_implementations/t8_geometry_linear.hxx>
#include <t8_eclass/t8_eclass.h>
#include <t8_types/t8_type.hxx>
#include <t8_types/t8_operators.hxx>

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
#include <map>
#include <tuple>
#include <array>

/**
 * Orientation of an OpenFOAM face with regard to the corresponding t8code face.
 */
struct t8_foam_face_orientation
{
  u_int8_t first_node;  /**< The foam face-vertex, which corresponds to the t8code face-vertex 0. */
  bool same_normal_dir; /**< True, if the face normals of the foam and the t8code face point in the same direction. */
};

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
  get_cell_eclass (const size_t cell_id);

  /**
   * Given a face eclass and orientation, this function converts the face id of a t8code face-vertex id to the corresponding
   * OpenFOAM (foam) face-vertex id.
   * \param [in] face_class         The eclass of the face (triangle or quad)
   * \param [in] orientation        The orientation of the foam face with regard to the t8code face
   * \param [in] t8_face_vertex_id  The t8 face-vertex id to convert
   * \return                        The converted foam face-vertex id
   */
  constexpr static u_int8_t
  t8_face_vertex_to_foam_point (t8_eclass_t face_class, const t8_foam_face_orientation orientation,
                                u_int8_t t8_face_vertex_id)
  {
    /* While t8code uses the z-order for quads, OpenFOAM uses an anti-clockwise numeration.
     * So if the face is a quad, we have to switch indices in accordance with this LUT.
     * This is not necessary for triangles. */
    if (face_class == T8_ECLASS_QUAD) {
      t8_face_vertex_id = quad_conversion[t8_face_vertex_id];
    }

    /* The OpenFOAM face point belonging to the current t8_face_vertex is calculated via the orientation
     * (orientation.same_normal_dir tells us if we count up +1 or down -1 depending on the fact that the face normals
     * point in the same direction and orientation.first_node tells us which foam point is located at t8 vertex 0)
     * and the current t8_face_vertex.
     * We add t8_eclass_num_vertices[face_class] once, so that we never take a mod of a negative value. */
    const int order = orientation.same_normal_dir ? 1 : -1;
    return (order * t8_face_vertex_id + orientation.first_node + t8_eclass_num_vertices[face_class])
           % t8_eclass_num_vertices[face_class];
  }

  /**
   * Search through a list of lists and find the first one that does not contain any "forbidden" values.
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
    /* If the search result is not the end iterator, return the position as distance from the start. */
    if (it != std::ranges::end (listlist)) {
      return std::distance (std::ranges::begin (listlist), it);
    }
    /* If no suitable list was found, return empty value. */
    return std::nullopt;
  }

  /**
   * Search through a list of lists and find the first one that contains all given required values.
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
    /* If the search result is not the end iterator, return the position as distance from the start. */
    if (it != std::ranges::end (listlist)) {
      return std::distance (std::ranges::begin (listlist), it);
    }
    /* If no suitable list was found, return empty value. */
    return std::nullopt;
  }

  /**
   * Build the cmesh from the raw mesh data.
   * \return  True on success
   */
  bool
  build_cmesh ();

  /**
   * Use the gathered cell information to build a hex cell and add it to the underlying cmesh.
   * \param [in] cell_id           The id of the cell.
   * \param [in] face_ids          The ids of the cell faces.
   * \param [in] face_point_ids    The ids of the points of the cell faces.
   * \param [in] face_normals      The normals of the cell faces. Using std::vector<char> instead of bools
   *                               because a vector of bools does not fulfill the std library requirements.
   */
  void
  reconstruct_hex_cell (size_t cell_id, std::vector<size_t> face_ids, std::vector<std::span<size_t>> face_point_ids,
                        std::vector<char> face_normals);

  /**
   * Computes the orientation of two faces according to \ref t8_cmesh_set_join.
   * \param [in] eclass1        Cell 1 eclass
   * \param [in] face_id1       Cell 1 face id
   * \param [in] orient1        Cell 1 orientation
   * \param [in] eclass2        Cell 2 eclass
   * \param [in] face_id2       Cell 2 face id
   * \param [in] orient2        cell 2 orientation
   * \return                    The orientation in accordance with \ref t8_cmesh_set_join
   */
  int
  compute_face_orientation (t8_eclass eclass1, u_int8_t face_id1, t8_foam_face_orientation orient1, t8_eclass eclass2,
                            u_int8_t face_id2, t8_foam_face_orientation orient2);

  /**
   * Compute the neighbors of a cell.
   * \tparam TNumFaces The number of faces of the cell.
   * \param [in] cell_id            The id of the cell.
   * \param [in] cell_eclass        The eclass of the cell.
   * \param [in] faces              Array containing the face ids of the cell.
   * \param [in] face_orientations  Array containing the orientation of each face (t8 vertex 0 -> foam point id)
   */
  template <size_t TNumFaces>
  void
  compute_cell_neighbors (const size_t cell_id, const t8_eclass cell_eclass,
                          const std::array<std::optional<size_t>, TNumFaces> faces,
                          const std::array<t8_foam_face_orientation, TNumFaces> face_orientations)
  {
    /* We are using \ref m_neighborhoods to link neighbors together.
     * Each face has a unique id. For each face of the cell we will save the cell id, t8 face id and orientation
     * in \ref m_neighborhoods. If we encounter a face, which is already filled then we know, that the already
     * saved face and cell is the neighbor. We can now join these cells via their shared face. */
    for (u_int8_t t8_face_id = 0; t8_face_id < TNumFaces; ++t8_face_id) {
      const size_t foam_face_id = faces[t8_face_id].value ();
      if (m_neighborhoods[foam_face_id].has_value ()) {
        /* There already is a cell saved with this face. This has to be the neighbor of the current cell. */
        const auto& [other_cell_id, other_cell_eclass, other_cell_t8_face_id, other_cell_face_orientation]
          = *m_neighborhoods[foam_face_id];

        /* Calculate the orientation between both cells. */
        const int orientation
          = compute_face_orientation (cell_eclass, t8_face_id, face_orientations[t8_face_id], other_cell_eclass,
                                      other_cell_t8_face_id, other_cell_face_orientation);
        t8_cmesh_set_join (m_cmesh, cell_id, other_cell_id, t8_face_id, other_cell_t8_face_id, orientation);
      }
      else {
        /* There is no cell with this face saved yet. Save this cell with this face. */
        m_neighborhoods[foam_face_id].emplace (cell_id, cell_eclass, t8_face_id, face_orientations[t8_face_id]);
      }
    }
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
   * We only need this for quads, since a triangle's corners are numbered equally.
   */
  static constexpr std::array quad_conversion = { 0, 1, 3, 2 };

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
  std::vector<std::optional<std::tuple<size_t, t8_eclass, u_int8_t, t8_foam_face_orientation>>> m_neighborhoods;
  /** Assigns a hash to each boundary face. */
  face_to_boundary_group_hash_map boundary_map;
};
