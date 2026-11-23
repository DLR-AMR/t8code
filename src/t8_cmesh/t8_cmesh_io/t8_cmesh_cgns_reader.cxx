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

/** \file t8_cmesh_cgns_reader.cxx
 * Implementation of the CGNS reader, see t8_cmesh_cgns_reader.hxx.
 *
 * The reader works in these steps:
 *  1. Read the coordinates and elements of all zones. Elements of the cell dimension become trees.
 *     In structured zones the cells of the index grid become trees.
 *     Every vertex gets a global id (zone offset + zone local id).
 *  2. Read the grid connectivities of all zones. Each of them yields pairs of corresponding vertices.
 *     Vertices of zone interfaces are the same vertex and are merged in a union-find structure.
 *     Vertices of periodic interfaces are kept apart and their correspondence is stored.
 *  3. Hash all tree faces by the (merged) ids of their vertices. Faces with the same key are neighbors.
 *     Boundary faces whose vertices correspond to the vertices of another boundary face across a periodic
 *     interface are neighbors as well.
 *  4. Resolve the boundary conditions of all zones to tree faces using the same face hash.
 *  5. Build and commit the cmesh.
 */

#include <t8_cmesh/t8_cmesh_io/t8_cmesh_cgns_reader.hxx>
#include <t8_cmesh/t8_cmesh.hxx>
#include <t8_cmesh/t8_cmesh_boundary_conditions/t8_cmesh_boundary_conditions.hxx>
#include <t8_cmesh/t8_cmesh_internal/t8_cmesh_types.h>
#include <t8_cmesh/t8_cmesh_vertex_connectivity/t8_cmesh_vertex_connectivity.h>
#include <t8_geometry/t8_geometry_with_vertices.h>
#include <t8_geometry/t8_geometry_implementations/t8_geometry_linear.hxx>
#include <t8_eclass/t8_eclass.h>

#include <cgnslib.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <cstring>
#include <limits>
#include <numeric>
#include <optional>
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <vector>

namespace
{

/** Length of name buffers for the CGNS mid level library (32 characters plus terminating zero). */
constexpr int T8_CGNS_NAME_LENGTH = 33;

/** Translate the CGNS corner vertex number of an element to the t8code vertex number.
 * The CGNS SIDS element ordering of the corner nodes coincides with the Gmsh ordering for all
 * supported element types, see https://cgns.github.io/standard/SIDS/convention.html#unstructured-grid-element-numbering-conventions
 */
constexpr int t8_cgns_vertex_to_t8_vertex[T8_ECLASS_COUNT][T8_ECLASS_MAX_CORNERS] = {
  { 0 },                      /* VERTEX */
  { 0, 1 },                   /* LINE */
  { 0, 1, 3, 2 },             /* QUAD */
  { 0, 1, 2 },                /* TRIANGLE */
  { 0, 1, 3, 2, 4, 5, 7, 6 }, /* HEX */
  { 0, 1, 3, 2 },             /* TET */
  { 0, 1, 2, 3, 4, 5 },       /* PRISM */
  { 0, 1, 3, 2, 4 }           /* PYRAMID */
};

/** Map a CGNS element type to a t8code element class.
 * For higher order elements only the corner nodes are used, which are always listed first in CGNS.
 * \return The eclass or T8_ECLASS_COUNT if the type is not supported.
 */
t8_eclass_t
t8_cgns_element_type_to_eclass (const CGNS_ENUMT (ElementType_t) type)
{
  switch (type) {
  case CGNS_ENUMV (NODE):
    return T8_ECLASS_VERTEX;
  case CGNS_ENUMV (BAR_2):
  case CGNS_ENUMV (BAR_3):
  case CGNS_ENUMV (BAR_4):
    return T8_ECLASS_LINE;
  case CGNS_ENUMV (TRI_3):
  case CGNS_ENUMV (TRI_6):
  case CGNS_ENUMV (TRI_9):
  case CGNS_ENUMV (TRI_10):
    return T8_ECLASS_TRIANGLE;
  case CGNS_ENUMV (QUAD_4):
  case CGNS_ENUMV (QUAD_8):
  case CGNS_ENUMV (QUAD_9):
  case CGNS_ENUMV (QUAD_12):
  case CGNS_ENUMV (QUAD_16):
    return T8_ECLASS_QUAD;
  case CGNS_ENUMV (TETRA_4):
  case CGNS_ENUMV (TETRA_10):
  case CGNS_ENUMV (TETRA_16):
  case CGNS_ENUMV (TETRA_20):
    return T8_ECLASS_TET;
  case CGNS_ENUMV (PYRA_5):
  case CGNS_ENUMV (PYRA_13):
  case CGNS_ENUMV (PYRA_14):
  case CGNS_ENUMV (PYRA_21):
  case CGNS_ENUMV (PYRA_29):
    return T8_ECLASS_PYRAMID;
  case CGNS_ENUMV (PENTA_6):
  case CGNS_ENUMV (PENTA_15):
  case CGNS_ENUMV (PENTA_18):
  case CGNS_ENUMV (PENTA_24):
  case CGNS_ENUMV (PENTA_38):
  case CGNS_ENUMV (PENTA_40):
    return T8_ECLASS_PRISM;
  case CGNS_ENUMV (HEXA_8):
  case CGNS_ENUMV (HEXA_20):
  case CGNS_ENUMV (HEXA_27):
  case CGNS_ENUMV (HEXA_32):
  case CGNS_ENUMV (HEXA_56):
  case CGNS_ENUMV (HEXA_64):
    return T8_ECLASS_HEX;
  default:
    return T8_ECLASS_COUNT;
  }
}

/** Check the return value of a CGNS call. On error print a message and return false from the calling function. */
#define T8_CGNS_CHECK(call) \
  do { \
    if ((call) != CG_OK) { \
      t8_global_errorf ("CGNS reader: %s failed: %s\n", #call, cg_get_error ()); \
      return false; \
    } \
  } while (0)

/** One element section (Elements_t node) of an unstructured zone. */
struct t8_cgns_section
{
  std::string name;                   /**< Name of the section. */
  CGNS_ENUMT (ElementType_t) type;    /**< Element type of the section, may be MIXED. */
  cgsize_t start;                     /**< First element id of the section (1-based, zone-wide). */
  cgsize_t end;                       /**< Last element id of the section. */
  std::vector<cgsize_t> connectivity; /**< ElementConnectivity array. Empty if the section was skipped. */
  std::vector<cgsize_t> offsets;      /**< ElementStartOffset array, only used for MIXED sections. */

  /** \return True if the section carries elements the reader can interpret. */
  bool
  is_readable () const
  {
    return !connectivity.empty ();
  }
};

/** The element class and the (1-based, zone local) vertex ids of one element of an unstructured zone. */
struct t8_cgns_element
{
  t8_eclass_t eclass;             /**< Element class, T8_ECLASS_COUNT if unsupported. */
  std::vector<cgsize_t> vertices; /**< Corner vertex ids in CGNS order. */
};

/** All data of one zone that the reader needs. */
struct t8_cgns_zone_data
{
  std::string name;                                  /**< Zone name. */
  bool structured = false;                           /**< True for structured zones. */
  int index_dim = 1;                                 /**< Number of indices per point (1 for unstructured zones). */
  std::array<cgsize_t, 3> vertex_dims = { 1, 1, 1 }; /**< Number of vertices per index direction (structured zones). */
  cgsize_t num_vertices = 0;                         /**< Number of vertices. */
  t8_gloidx_t vertex_offset = 0;                     /**< Global id of the zone local vertex 1. */
  t8_gloidx_t first_tree = 0;                        /**< Global id of the first tree of this zone. */
  t8_gloidx_t num_trees = 0;                         /**< Number of trees of this zone. */
  std::vector<t8_cgns_section> sections;             /**< Element sections (unstructured zones). */

  /** Convert a structured index tuple (1-based) to the zone local vertex id (1-based). */
  cgsize_t
  structured_vertex_id (const cgsize_t *ijk) const
  {
    T8_ASSERT (structured);
    cgsize_t id = ijk[0];
    cgsize_t stride = vertex_dims[0];
    for (int dim = 1; dim < index_dim; ++dim) {
      id += (ijk[dim] - 1) * stride;
      stride *= vertex_dims[dim];
    }
    return id;
  }

  /** Get the class and vertices of an element of an unstructured zone.
   * \param [in] element_id  The zone wide element id (1-based).
   * \return                 The element, or nothing if the id is not part of a readable section.
   */
  std::optional<t8_cgns_element>
  get_element (const cgsize_t element_id) const
  {
    for (const t8_cgns_section &section : sections) {
      if (element_id < section.start || element_id > section.end || !section.is_readable ()) {
        continue;
      }
      const cgsize_t local = element_id - section.start;
      CGNS_ENUMT (ElementType_t) type = section.type;
      const cgsize_t *nodes = nullptr;
      if (type == CGNS_ENUMV (MIXED)) {
        const cgsize_t offset = section.offsets[local];
        type = static_cast<CGNS_ENUMT (ElementType_t)> (section.connectivity[offset]);
        nodes = section.connectivity.data () + offset + 1;
      }
      else {
        int npe = 0;
        cg_npe (type, &npe);
        nodes = section.connectivity.data () + local * npe;
      }
      const t8_eclass_t eclass = t8_cgns_element_type_to_eclass (type);
      t8_cgns_element element;
      element.eclass = eclass;
      if (eclass != T8_ECLASS_COUNT) {
        const int num_corners = t8_eclass_num_vertices[eclass];
        element.vertices.assign (nodes, nodes + num_corners);
      }
      return element;
    }
    return std::nullopt;
  }

  /** Convert a point set of this zone to a list of zone local ids (1-based).
   * For unstructured zones the ids are the vertex or element ids as stored.
   * For structured zones the index tuples are converted to vertex ids, which is only meaningful
   * for point sets with GridLocation Vertex.
   * \param [in] type        PointList or PointRange (also the donor variants).
   * \param [in] raw_points  The point set as read from the file, \a index_dim entries per point.
   * \param [in] num_points  Number of points in \a raw_points.
   */
  std::vector<cgsize_t>
  point_set_to_ids (const CGNS_ENUMT (PointSetType_t) type, const std::vector<cgsize_t> &raw_points,
                    const cgsize_t num_points) const
  {
    const bool is_range = type == CGNS_ENUMV (PointRange) || type == CGNS_ENUMV (PointRangeDonor);
    std::vector<cgsize_t> ids;
    if (!structured) {
      if (is_range) {
        T8_ASSERT (num_points >= 2);
        for (cgsize_t id = raw_points[0]; id <= raw_points[1]; ++id) {
          ids.push_back (id);
        }
      }
      else {
        ids.assign (raw_points.begin (), raw_points.begin () + num_points);
      }
      return ids;
    }
    if (is_range) {
      T8_ASSERT (num_points >= 2);
      /* The range is a box between the two index tuples. */
      std::array<cgsize_t, 3> low = { 1, 1, 1 }, high = { 1, 1, 1 };
      for (int dim = 0; dim < index_dim; ++dim) {
        low[dim] = std::min (raw_points[dim], raw_points[index_dim + dim]);
        high[dim] = std::max (raw_points[dim], raw_points[index_dim + dim]);
      }
      for (cgsize_t k = low[2]; k <= high[2]; ++k) {
        for (cgsize_t j = low[1]; j <= high[1]; ++j) {
          for (cgsize_t i = low[0]; i <= high[0]; ++i) {
            const cgsize_t ijk[3] = { i, j, k };
            ids.push_back (structured_vertex_id (ijk));
          }
        }
      }
      return ids;
    }
    for (cgsize_t ipoint = 0; ipoint < num_points; ++ipoint) {
      ids.push_back (structured_vertex_id (raw_points.data () + ipoint * index_dim));
    }
    return ids;
  }
};

/** A tree of the cmesh under construction. */
struct t8_cgns_tree
{
  t8_eclass_t eclass;                                      /**< Element class. */
  int zone;                                                /**< Index of the zone the tree belongs to. */
  std::array<t8_gloidx_t, T8_ECLASS_MAX_CORNERS> vertices; /**< Global vertex ids in t8code vertex order. */
};

/** Union-find structure over the global vertex ids. */
class t8_cgns_union_find {
 public:
  explicit t8_cgns_union_find (const t8_gloidx_t num_elements): m_parent (num_elements)
  {
    std::iota (m_parent.begin (), m_parent.end (), 0);
  }

  t8_gloidx_t
  find (t8_gloidx_t id)
  {
    while (m_parent[id] != id) {
      m_parent[id] = m_parent[m_parent[id]];
      id = m_parent[id];
    }
    return id;
  }

  void
  unite (const t8_gloidx_t a, const t8_gloidx_t b)
  {
    const t8_gloidx_t root_a = find (a);
    const t8_gloidx_t root_b = find (b);
    if (root_a != root_b) {
      /* Attach the larger root to the smaller one so that representatives stay small. */
      m_parent[std::max (root_a, root_b)] = std::min (root_a, root_b);
    }
  }

 private:
  std::vector<t8_gloidx_t> m_parent;
};

/** Key of a tree face: the sorted (merged) vertex ids, padded with -1 at the front. */
using t8_cgns_face_key = std::array<t8_gloidx_t, T8_ECLASS_MAX_CORNERS_2D>;

struct t8_cgns_face_key_hash
{
  size_t
  operator() (const t8_cgns_face_key &key) const
  {
    size_t hash = 0;
    for (const t8_gloidx_t id : key) {
      hash ^= std::hash<t8_gloidx_t> {}(id) + 0x9e3779b97f4a7c15ULL + (hash << 6) + (hash >> 2);
    }
    return hash;
  }
};

/** One tree face in the face hash. */
struct t8_cgns_face_entry
{
  t8_gloidx_t tree;          /**< Global tree id. */
  int face;                  /**< Face number of the tree. */
  int num_vertices;          /**< Number of vertices of the face. */
  t8_cgns_face_key vertices; /**< Merged vertex ids in t8code face vertex order, padded with -1. */
};

using t8_cgns_face_map = std::unordered_map<t8_cgns_face_key, std::vector<t8_cgns_face_entry>, t8_cgns_face_key_hash>;

/** Build the key of a face from its vertex ids (padded with -1). */
t8_cgns_face_key
t8_cgns_make_face_key (const t8_cgns_face_key &vertices)
{
  t8_cgns_face_key key = vertices;
  std::sort (key.begin (), key.end ());
  return key;
}

/** \return True if two vertices of the face coincide. */
bool
t8_cgns_face_key_is_degenerate (const t8_cgns_face_key &key, const int num_vertices)
{
  /* The key is sorted, the padding of -1 is at the front. */
  const auto first = key.end () - num_vertices;
  return std::adjacent_find (first, key.end ()) != key.end ();
}

/** The vertex correspondence of one periodic interface.
 * Periodic vertices are not merged, since the two sides of a periodic interface may belong to the same tree
 * (for example in a domain that is only one cell wide). Instead the tree faces on both sides are joined
 * explicitly using this correspondence. */
struct t8_cgns_periodic_interface
{
  std::string name;                                              /**< Name of the connectivity node. */
  int zone;                                                      /**< Zone of the first vertex of each pair. */
  std::vector<std::pair<t8_gloidx_t, t8_gloidx_t>> vertex_pairs; /**< Global vertex ids (vertex, donor vertex). */
};

/** Everything the reader collects from the file. */
struct t8_cgns_file_data
{
  int cell_dim = 0;                             /**< Cell dimension of the base. */
  int phys_dim = 0;                             /**< Physical dimension of the base. */
  std::vector<t8_cgns_zone_data> zones;         /**< All zones of the base. */
  std::vector<double> coordinates;              /**< 3 coordinates per global vertex. */
  std::vector<t8_cgns_tree> trees;              /**< All trees. */
  std::optional<t8_cgns_union_find> merged_ids; /**< Merges vertices that coincide across zone interfaces. */
  std::vector<t8_cgns_periodic_interface> periodic_interfaces; /**< Vertex correspondences of periodic interfaces. */
  t8_cgns_face_map faces; /**< Hash of all tree faces by their merged vertex ids. */
  std::vector<std::array<std::string, T8_ECLASS_MAX_FACES>>
    face_labels;                        /**< Boundary condition label of every tree face. */
  bool has_boundary_conditions = false; /**< True if any BC_t node was found. */

  /** \return The global number of vertices of all zones. */
  t8_gloidx_t
  num_global_vertices () const
  {
    return static_cast<t8_gloidx_t> (coordinates.size () / 3);
  }

  /** Convert a zone local vertex id (1-based) to a global vertex id. */
  t8_gloidx_t
  global_vertex (const int zone, const cgsize_t local_id) const
  {
    T8_ASSERT (local_id >= 1 && local_id <= zones[zone].num_vertices);
    return zones[zone].vertex_offset + static_cast<t8_gloidx_t> (local_id) - 1;
  }

  /** Get the coordinates of a global vertex. */
  std::array<double, 3>
  vertex_coordinates (const t8_gloidx_t global_id) const
  {
    return { coordinates[3 * global_id], coordinates[3 * global_id + 1], coordinates[3 * global_id + 2] };
  }

  /** Find a zone by its name.
   * \return The zone index or -1 if no zone has that name.
   */
  int
  find_zone (const char *name) const
  {
    for (size_t izone = 0; izone < zones.size (); ++izone) {
      if (zones[izone].name == name) {
        return static_cast<int> (izone);
      }
    }
    return -1;
  }

  /** Record that the global vertex \a a of an interface corresponds to the global vertex \a b of the donor
   * interface. Vertices of non-periodic interfaces are the same vertex and are merged. Vertices of periodic
   * interfaces are kept apart and the correspondence is stored in \a periodic (if not nullptr). */
  void
  merge_vertices (const t8_gloidx_t a, const t8_gloidx_t b, t8_cgns_periodic_interface *periodic)
  {
    if (periodic == nullptr) {
      merged_ids->unite (a, b);
    }
    else {
      periodic->vertex_pairs.emplace_back (a, b);
    }
  }

  /** Store a periodic interface, if it has any vertex pairs. */
  void
  add_periodic_interface (std::optional<t8_cgns_periodic_interface> &periodic)
  {
    if (periodic && !periodic->vertex_pairs.empty ()) {
      periodic_interfaces.push_back (std::move (*periodic));
    }
  }
};

/** Read the coordinates of a zone into the global coordinate array. */
bool
t8_cgns_read_coordinates (const int fn, const int base, const int zone_index, const t8_cgns_zone_data &zone,
                          t8_cgns_file_data &data)
{
  int num_coords = 0;
  T8_CGNS_CHECK (cg_ncoords (fn, base, zone_index, &num_coords));
  const cgsize_t range_min[3] = { 1, 1, 1 };
  const cgsize_t range_max[3]
    = { zone.structured ? zone.vertex_dims[0] : zone.num_vertices, zone.vertex_dims[1], zone.vertex_dims[2] };
  const t8_gloidx_t offset = zone.vertex_offset;
  T8_ASSERT (data.coordinates.size () >= 3 * static_cast<size_t> (offset + zone.num_vertices));
  std::vector<double> buffer (zone.num_vertices);
  const std::array<const char *, 3> coordinate_names = { "CoordinateX", "CoordinateY", "CoordinateZ" };
  std::array<bool, 3> found = { false, false, false };
  for (int icoord = 1; icoord <= num_coords; ++icoord) {
    char name[T8_CGNS_NAME_LENGTH];
    CGNS_ENUMT (DataType_t) datatype;
    T8_CGNS_CHECK (cg_coord_info (fn, base, zone_index, icoord, &datatype, name));
    for (int axis = 0; axis < 3; ++axis) {
      if (strcmp (name, coordinate_names[axis]) != 0) {
        continue;
      }
      /* The library converts the stored precision to double. */
      T8_CGNS_CHECK (
        cg_coord_read (fn, base, zone_index, name, CGNS_ENUMV (RealDouble), range_min, range_max, buffer.data ()));
      for (cgsize_t ivertex = 0; ivertex < zone.num_vertices; ++ivertex) {
        data.coordinates[3 * (offset + ivertex) + axis] = buffer[ivertex];
      }
      found[axis] = true;
    }
  }
  for (int axis = 0; axis < data.phys_dim; ++axis) {
    if (!found[axis]) {
      t8_global_errorf ("CGNS reader: Zone %s does not provide %s. Only Cartesian coordinates are supported.\n",
                        zone.name.c_str (), coordinate_names[axis]);
      return false;
    }
  }
  return true;
}

/** Read all element sections of an unstructured zone. */
bool
t8_cgns_read_sections (const int fn, const int base, const int zone_index, t8_cgns_zone_data &zone)
{
  int num_sections = 0;
  T8_CGNS_CHECK (cg_nsections (fn, base, zone_index, &num_sections));
  zone.sections.resize (num_sections);
  for (int isection = 1; isection <= num_sections; ++isection) {
    t8_cgns_section &section = zone.sections[isection - 1];
    char name[T8_CGNS_NAME_LENGTH];
    int num_boundary = 0;
    int parent_flag = 0;
    T8_CGNS_CHECK (cg_section_read (fn, base, zone_index, isection, name, &section.type, &section.start, &section.end,
                                    &num_boundary, &parent_flag));
    section.name = name;
    if (section.type == CGNS_ENUMV (NGON_n) || section.type == CGNS_ENUMV (NFACE_n)) {
      t8_global_errorf ("CGNS reader: Warning: Skipping section %s of zone %s. Polyhedral sections (%s) are not "
                        "supported.\n",
                        section.name.c_str (), zone.name.c_str (), cg_ElementTypeName (section.type));
      continue;
    }
    if (section.type != CGNS_ENUMV (MIXED) && t8_cgns_element_type_to_eclass (section.type) == T8_ECLASS_COUNT) {
      t8_global_errorf ("CGNS reader: Warning: Skipping section %s of zone %s. Element type %s is not supported.\n",
                        section.name.c_str (), zone.name.c_str (), cg_ElementTypeName (section.type));
      continue;
    }
    cgsize_t data_size = 0;
    T8_CGNS_CHECK (cg_ElementDataSize (fn, base, zone_index, isection, &data_size));
    section.connectivity.resize (data_size);
    if (section.type == CGNS_ENUMV (MIXED)) {
      section.offsets.resize (section.end - section.start + 2);
      T8_CGNS_CHECK (cg_poly_elements_read (fn, base, zone_index, isection, section.connectivity.data (),
                                            section.offsets.data (), nullptr));
    }
    else {
      T8_CGNS_CHECK (cg_elements_read (fn, base, zone_index, isection, section.connectivity.data (), nullptr));
    }
    t8_debugf ("CGNS reader: Read section %s of zone %s, elements %lld to %lld of type %s.\n", section.name.c_str (),
               zone.name.c_str (), static_cast<long long> (section.start), static_cast<long long> (section.end),
               cg_ElementTypeName (section.type));
  }
  return true;
}

/** Create the trees of an unstructured zone from its sections. */
void
t8_cgns_create_unstructured_trees (const int zone_number, const t8_cgns_zone_data &zone, t8_cgns_file_data &data)
{
  for (const t8_cgns_section &section : zone.sections) {
    if (!section.is_readable ()) {
      continue;
    }
    for (cgsize_t element_id = section.start; element_id <= section.end; ++element_id) {
      const std::optional<t8_cgns_element> element = zone.get_element (element_id);
      T8_ASSERT (element.has_value ());
      if (element->eclass == T8_ECLASS_COUNT || t8_eclass_to_dimension[element->eclass] != data.cell_dim) {
        continue;
      }
      t8_cgns_tree tree;
      tree.eclass = element->eclass;
      tree.zone = zone_number;
      tree.vertices.fill (-1);
      const int num_corners = t8_eclass_num_vertices[tree.eclass];
      for (int icorner = 0; icorner < num_corners; ++icorner) {
        const int t8_corner = t8_cgns_vertex_to_t8_vertex[tree.eclass][icorner];
        tree.vertices[t8_corner] = data.global_vertex (zone_number, element->vertices[icorner]);
      }
      data.trees.push_back (tree);
    }
  }
}

/** Create the trees of a structured zone: one line, quad or hex per cell of the index grid.
 * The cells are numbered with the first index running fastest. */
void
t8_cgns_create_structured_trees (const int zone_number, const t8_cgns_zone_data &zone, t8_cgns_file_data &data)
{
  const int dim = zone.index_dim;
  const t8_eclass_t eclass = dim == 3 ? T8_ECLASS_HEX : (dim == 2 ? T8_ECLASS_QUAD : T8_ECLASS_LINE);
  std::array<cgsize_t, 3> num_cells = { 1, 1, 1 };
  for (int idim = 0; idim < dim; ++idim) {
    num_cells[idim] = zone.vertex_dims[idim] - 1;
  }
  const int num_corners = t8_eclass_num_vertices[eclass];
  for (cgsize_t k = 0; k < num_cells[2]; ++k) {
    for (cgsize_t j = 0; j < num_cells[1]; ++j) {
      for (cgsize_t i = 0; i < num_cells[0]; ++i) {
        t8_cgns_tree tree;
        tree.eclass = eclass;
        tree.zone = zone_number;
        tree.vertices.fill (-1);
        /* The t8code vertex numbering of lines, quads and hexes is the z-order of the corners. */
        for (int icorner = 0; icorner < num_corners; ++icorner) {
          const cgsize_t ijk[3] = { i + 1 + (icorner & 1), j + 1 + ((icorner >> 1) & 1), k + 1 + ((icorner >> 2) & 1) };
          tree.vertices[icorner] = data.global_vertex (zone_number, zone.structured_vertex_id (ijk));
        }
        data.trees.push_back (tree);
      }
    }
  }
}

/** Read the base and all zones with their coordinates and elements. Fills the trees of \a data. */
bool
t8_cgns_read_zones (const int fn, const int base, t8_cgns_file_data &data)
{
  char base_name[T8_CGNS_NAME_LENGTH];
  T8_CGNS_CHECK (cg_base_read (fn, base, base_name, &data.cell_dim, &data.phys_dim));
  if (data.cell_dim < 1 || data.cell_dim > 3 || data.phys_dim < data.cell_dim || data.phys_dim > 3) {
    t8_global_errorf ("CGNS reader: Base %s has cell dimension %d and physical dimension %d, which is not supported.\n",
                      base_name, data.cell_dim, data.phys_dim);
    return false;
  }
  int num_zones = 0;
  T8_CGNS_CHECK (cg_nzones (fn, base, &num_zones));
  if (num_zones < 1) {
    t8_global_errorf ("CGNS reader: Base %s does not contain any zone.\n", base_name);
    return false;
  }
  data.zones.resize (num_zones);
  /* First pass: zone names and sizes, so that the global vertex offsets are known. */
  t8_gloidx_t vertex_offset = 0;
  for (int izone = 1; izone <= num_zones; ++izone) {
    t8_cgns_zone_data &zone = data.zones[izone - 1];
    CGNS_ENUMT (ZoneType_t) zone_type;
    T8_CGNS_CHECK (cg_zone_type (fn, base, izone, &zone_type));
    T8_CGNS_CHECK (cg_index_dim (fn, base, izone, &zone.index_dim));
    char zone_name[T8_CGNS_NAME_LENGTH];
    cgsize_t size[9] = { 0 };
    T8_CGNS_CHECK (cg_zone_read (fn, base, izone, zone_name, size));
    zone.name = zone_name;
    if (zone_type == CGNS_ENUMV (Structured)) {
      zone.structured = true;
      if (zone.index_dim != data.cell_dim) {
        t8_global_errorf (
          "CGNS reader: Structured zone %s has index dimension %d but the base has cell dimension %d.\n", zone_name,
          zone.index_dim, data.cell_dim);
        return false;
      }
      zone.num_vertices = 1;
      for (int idim = 0; idim < zone.index_dim; ++idim) {
        zone.vertex_dims[idim] = size[idim];
        zone.num_vertices *= size[idim];
        if (size[idim] < 2) {
          t8_global_errorf ("CGNS reader: Structured zone %s has only %lld vertices in index direction %d.\n",
                            zone_name, static_cast<long long> (size[idim]), idim);
          return false;
        }
      }
    }
    else if (zone_type == CGNS_ENUMV (Unstructured)) {
      zone.structured = false;
      zone.index_dim = 1;
      zone.num_vertices = size[0];
    }
    else {
      t8_global_errorf ("CGNS reader: Zone %s is of type %s, which is not supported.\n", zone_name,
                        cg_ZoneTypeName (zone_type));
      return false;
    }
    zone.vertex_offset = vertex_offset;
    vertex_offset += zone.num_vertices;
  }
  data.coordinates.assign (3 * vertex_offset, 0.0);
  /* Second pass: coordinates, sections and trees. */
  for (int izone = 1; izone <= num_zones; ++izone) {
    t8_cgns_zone_data &zone = data.zones[izone - 1];
    if (!t8_cgns_read_coordinates (fn, base, izone, zone, data)) {
      return false;
    }
    zone.first_tree = static_cast<t8_gloidx_t> (data.trees.size ());
    if (zone.structured) {
      t8_cgns_create_structured_trees (izone - 1, zone, data);
    }
    else {
      if (!t8_cgns_read_sections (fn, base, izone, zone)) {
        return false;
      }
      t8_cgns_create_unstructured_trees (izone - 1, zone, data);
    }
    zone.num_trees = static_cast<t8_gloidx_t> (data.trees.size ()) - zone.first_tree;
    t8_debugf ("CGNS reader: %s zone %s has %lld vertices and %lld trees.\n",
               zone.structured ? "Structured" : "Unstructured", zone.name.c_str (),
               static_cast<long long> (zone.num_vertices), static_cast<long long> (zone.num_trees));
  }
  if (data.trees.empty ()) {
    t8_global_errorf ("CGNS reader: Base %s does not contain any element of dimension %d.\n", base_name, data.cell_dim);
    return false;
  }
  return true;
}

/** A rigid transformation p -> R (p - c) + c + t as stored in a periodic grid connectivity. */
struct t8_cgns_periodic_transform
{
  std::array<double, 3> center;      /**< RotationCenter */
  std::array<double, 3> angle;       /**< RotationAngle in radians about the x, y and z axis */
  std::array<double, 3> translation; /**< Translation */

  /** Fill from the single precision arrays of the CGNS library. */
  void
  set (const float *rotation_center, const float *rotation_angle, const float *translation_vector)
  {
    for (int i = 0; i < 3; ++i) {
      center[i] = rotation_center[i];
      angle[i] = rotation_angle[i];
      translation[i] = translation_vector[i];
    }
  }

  /** Build the rotation matrix R = Rz * Ry * Rx. */
  std::array<std::array<double, 3>, 3>
  rotation () const
  {
    const double cx = std::cos (angle[0]), sx = std::sin (angle[0]);
    const double cy = std::cos (angle[1]), sy = std::sin (angle[1]);
    const double cz = std::cos (angle[2]), sz = std::sin (angle[2]);
    return { { { cy * cz, sx * sy * cz - cx * sz, cx * sy * cz + sx * sz },
               { cy * sz, sx * sy * sz + cx * cz, cx * sy * sz - sx * cz },
               { -sy, sx * cy, cx * cy } } };
  }

  /** Apply the transformation. If \a inverse is true, apply the inverse transformation. */
  std::array<double, 3>
  apply (const std::array<double, 3> &point, const bool inverse) const
  {
    const auto R = rotation ();
    std::array<double, 3> result = { 0, 0, 0 };
    if (!inverse) {
      for (int i = 0; i < 3; ++i) {
        for (int j = 0; j < 3; ++j) {
          result[i] += R[i][j] * (point[j] - center[j]);
        }
        result[i] += center[i] + translation[i];
      }
    }
    else {
      for (int i = 0; i < 3; ++i) {
        for (int j = 0; j < 3; ++j) {
          result[i] += R[j][i] * (point[j] - center[j] - translation[j]);
        }
        result[i] += center[i];
      }
    }
    return result;
  }
};

double
t8_cgns_distance (const std::array<double, 3> &a, const std::array<double, 3> &b)
{
  return std::sqrt ((a[0] - b[0]) * (a[0] - b[0]) + (a[1] - b[1]) * (a[1] - b[1]) + (a[2] - b[2]) * (a[2] - b[2]));
}

/** Match the vertices of two face elements geometrically.
 * \param [in] vertices_a  Global vertex ids of the first face.
 * \param [in] vertices_b  Global vertex ids of the second face. Must have the same size as \a vertices_a.
 * \param [in] transform   If given, the vertices of the first face are transformed before matching.
 * \param [in] inverse     Apply the inverse transformation.
 * \param [in] data        The file data (for the coordinates).
 * \param [out] matching   On success matching[i] is the index into \a vertices_b that corresponds to vertices_a[i].
 * \return                 True if a one-to-one matching within tolerance was found.
 */
bool
t8_cgns_match_face_vertices (const std::vector<t8_gloidx_t> &vertices_a, const std::vector<t8_gloidx_t> &vertices_b,
                             const t8_cgns_periodic_transform *transform, const bool inverse,
                             const t8_cgns_file_data &data, std::vector<int> &matching)
{
  const size_t num_vertices = vertices_a.size ();
  T8_ASSERT (vertices_b.size () == num_vertices);
  /* Tolerance relative to the diameter of the second face. */
  double diameter = 0.0;
  for (size_t i = 0; i < num_vertices; ++i) {
    for (size_t j = i + 1; j < num_vertices; ++j) {
      diameter = std::max (
        diameter, t8_cgns_distance (data.vertex_coordinates (vertices_b[i]), data.vertex_coordinates (vertices_b[j])));
    }
  }
  /* Periodic transformations are stored in single precision in the file. */
  const double tolerance = 1e-4 * (diameter > 0.0 ? diameter : 1.0);
  matching.assign (num_vertices, -1);
  std::vector<bool> used (num_vertices, false);
  for (size_t i = 0; i < num_vertices; ++i) {
    std::array<double, 3> point = data.vertex_coordinates (vertices_a[i]);
    if (transform != nullptr) {
      point = transform->apply (point, inverse);
    }
    double best_distance = std::numeric_limits<double>::max ();
    int best = -1;
    for (size_t j = 0; j < num_vertices; ++j) {
      if (used[j]) {
        continue;
      }
      const double distance = t8_cgns_distance (point, data.vertex_coordinates (vertices_b[j]));
      if (distance < best_distance) {
        best_distance = distance;
        best = static_cast<int> (j);
      }
    }
    if (best < 0 || best_distance > tolerance) {
      return false;
    }
    used[best] = true;
    matching[i] = best;
  }
  return true;
}

/** Read the general grid connectivities (GridConnectivity_t) of a zone and merge the
 * corresponding vertices in the union-find structures. */
bool
t8_cgns_read_connectivities (const int fn, const int base, const int zone_index, t8_cgns_file_data &data)
{
  const int zone_number = zone_index - 1;
  const t8_cgns_zone_data &zone = data.zones[zone_number];
  int num_connectivities = 0;
  T8_CGNS_CHECK (cg_nconns (fn, base, zone_index, &num_connectivities));
  for (int iconn = 1; iconn <= num_connectivities; ++iconn) {
    char name[T8_CGNS_NAME_LENGTH];
    char donor_name[T8_CGNS_NAME_LENGTH];
    CGNS_ENUMT (GridLocation_t) location;
    CGNS_ENUMT (GridConnectivityType_t) conn_type;
    CGNS_ENUMT (PointSetType_t) point_set_type;
    CGNS_ENUMT (PointSetType_t) donor_point_set_type;
    CGNS_ENUMT (ZoneType_t) donor_zone_type;
    CGNS_ENUMT (DataType_t) donor_data_type;
    cgsize_t num_points = 0;
    cgsize_t num_donor_points = 0;
    T8_CGNS_CHECK (cg_conn_info (fn, base, zone_index, iconn, name, &location, &conn_type, &point_set_type, &num_points,
                                 donor_name, &donor_zone_type, &donor_point_set_type, &donor_data_type,
                                 &num_donor_points));
    if (conn_type != CGNS_ENUMV (Abutting1to1)) {
      t8_global_errorf ("CGNS reader: Warning: Skipping connectivity %s of zone %s of type %s. Only Abutting1to1 "
                        "interfaces are supported.\n",
                        name, zone.name.c_str (), cg_GridConnectivityTypeName (conn_type));
      continue;
    }
    const int donor_zone_number = data.find_zone (donor_name);
    if (donor_zone_number < 0) {
      t8_global_errorf ("CGNS reader: Warning: Skipping connectivity %s of zone %s. Donor zone %s not found.\n", name,
                        zone.name.c_str (), donor_name);
      continue;
    }
    const t8_cgns_zone_data &donor_zone = data.zones[donor_zone_number];
    const bool vertex_location = location == CGNS_ENUMV (Vertex);
    const bool face_location = location == CGNS_ENUMV (FaceCenter) || location == CGNS_ENUMV (EdgeCenter)
                               || location == CGNS_ENUMV (CellCenter);
    if (!vertex_location && !(face_location && !zone.structured && !donor_zone.structured)) {
      t8_global_errorf ("CGNS reader: Warning: Skipping connectivity %s of zone %s. GridLocation %s is not supported "
                        "for this zone type.\n",
                        name, zone.name.c_str (), cg_GridLocationName (location));
      continue;
    }
    std::vector<cgsize_t> raw_points (std::max<cgsize_t> (num_points * zone.index_dim, 1));
    std::vector<cgsize_t> raw_donor_points (std::max<cgsize_t> (num_donor_points * donor_zone.index_dim, 1));
    T8_CGNS_CHECK (
      cg_conn_read (fn, base, zone_index, iconn, raw_points.data (), donor_data_type, raw_donor_points.data ()));
    const std::vector<cgsize_t> points = zone.point_set_to_ids (point_set_type, raw_points, num_points);
    const std::vector<cgsize_t> donor_points
      = donor_zone.point_set_to_ids (donor_point_set_type, raw_donor_points, num_donor_points);
    if (points.size () != donor_points.size ()) {
      t8_global_errorf ("CGNS reader: Warning: Skipping connectivity %s of zone %s. It has %zu points but %zu donor "
                        "points.\n",
                        name, zone.name.c_str (), points.size (), donor_points.size ());
      continue;
    }
    /* Periodic interfaces carry a transformation from this interface to the donor interface. */
    std::optional<t8_cgns_periodic_transform> transform;
    {
      float center[3], angle[3], translation[3];
      if (cg_conn_periodic_read (fn, base, zone_index, iconn, center, angle, translation) == CG_OK) {
        transform.emplace ();
        transform->set (center, angle, translation);
      }
    }
    std::optional<t8_cgns_periodic_interface> periodic;
    if (transform) {
      periodic.emplace ();
      periodic->name = name;
      periodic->zone = zone_number;
    }
    t8_debugf ("CGNS reader: Reading %sconnectivity %s from zone %s to zone %s with %zu points at %s.\n",
               periodic ? "periodic " : "", name, zone.name.c_str (), donor_name, points.size (),
               cg_GridLocationName (location));
    if (vertex_location) {
      /* The point lists are vertex ids with a one-to-one correspondence. */
      for (size_t ipoint = 0; ipoint < points.size (); ++ipoint) {
        data.merge_vertices (data.global_vertex (zone_number, points[ipoint]),
                             data.global_vertex (donor_zone_number, donor_points[ipoint]),
                             periodic ? &*periodic : nullptr);
      }
      data.add_periodic_interface (periodic);
      continue;
    }
    /* The point lists are face element ids. Their vertices are matched geometrically. */
    size_t num_skipped = 0;
    for (size_t ipoint = 0; ipoint < points.size (); ++ipoint) {
      const std::optional<t8_cgns_element> element = zone.get_element (points[ipoint]);
      const std::optional<t8_cgns_element> donor_element = donor_zone.get_element (donor_points[ipoint]);
      if (!element || !donor_element || element->eclass == T8_ECLASS_COUNT || donor_element->eclass == T8_ECLASS_COUNT
          || t8_eclass_to_dimension[element->eclass] != data.cell_dim - 1
          || element->vertices.size () != donor_element->vertices.size ()) {
        ++num_skipped;
        continue;
      }
      std::vector<t8_gloidx_t> vertices, donor_vertices;
      for (const cgsize_t id : element->vertices) {
        vertices.push_back (data.global_vertex (zone_number, id));
      }
      for (const cgsize_t id : donor_element->vertices) {
        donor_vertices.push_back (data.global_vertex (donor_zone_number, id));
      }
      std::vector<int> matching;
      bool matched = t8_cgns_match_face_vertices (vertices, donor_vertices, transform ? &*transform : nullptr, false,
                                                  data, matching);
      if (!matched && transform) {
        /* The file might store the transformation in the opposite direction. */
        matched = t8_cgns_match_face_vertices (vertices, donor_vertices, &*transform, true, data, matching);
      }
      if (!matched) {
        ++num_skipped;
        continue;
      }
      for (size_t ivertex = 0; ivertex < vertices.size (); ++ivertex) {
        data.merge_vertices (vertices[ivertex], donor_vertices[matching[ivertex]], periodic ? &*periodic : nullptr);
      }
    }
    data.add_periodic_interface (periodic);
    if (num_skipped > 0) {
      t8_global_errorf ("CGNS reader: Warning: %zu of %zu faces of connectivity %s of zone %s could not be matched "
                        "with their donor faces and were skipped.\n",
                        num_skipped, points.size (), name, zone.name.c_str ());
    }
  }
  return true;
}

/** Read the one-to-one connectivities (GridConnectivity1to1_t) of a structured zone and merge the
 * corresponding vertices. The donor index of a point is
 *   donor_index = T * (index - begin) + donor_begin,
 * where the matrix T is given by the Transform vector, see the CGNS SIDS. */
bool
t8_cgns_read_1to1_connectivities (const int fn, const int base, const int zone_index, t8_cgns_file_data &data)
{
  const int zone_number = zone_index - 1;
  const t8_cgns_zone_data &zone = data.zones[zone_number];
  if (!zone.structured) {
    return true;
  }
  int num_connectivities = 0;
  T8_CGNS_CHECK (cg_n1to1 (fn, base, zone_index, &num_connectivities));
  const int dim = zone.index_dim;
  for (int iconn = 1; iconn <= num_connectivities; ++iconn) {
    char name[T8_CGNS_NAME_LENGTH];
    char donor_name[T8_CGNS_NAME_LENGTH];
    cgsize_t range[6] = { 0 };
    cgsize_t donor_range[6] = { 0 };
    int transform[3] = { 0 };
    T8_CGNS_CHECK (cg_1to1_read (fn, base, zone_index, iconn, name, donor_name, range, donor_range, transform));
    const int donor_zone_number = data.find_zone (donor_name);
    if (donor_zone_number < 0 || !data.zones[donor_zone_number].structured) {
      t8_global_errorf ("CGNS reader: Warning: Skipping connectivity %s of zone %s. Donor zone %s not found or not "
                        "structured.\n",
                        name, zone.name.c_str (), donor_name);
      continue;
    }
    const t8_cgns_zone_data &donor_zone = data.zones[donor_zone_number];
    bool valid_transform = true;
    for (int idim = 0; idim < dim; ++idim) {
      if (std::abs (transform[idim]) < 1 || std::abs (transform[idim]) > dim) {
        valid_transform = false;
      }
    }
    if (!valid_transform) {
      t8_global_errorf ("CGNS reader: Warning: Skipping connectivity %s of zone %s. Invalid Transform [%d, %d, %d].\n",
                        name, zone.name.c_str (), transform[0], transform[1], transform[2]);
      continue;
    }
    float center[3], angle[3], translation[3];
    std::optional<t8_cgns_periodic_interface> periodic;
    if (cg_1to1_periodic_read (fn, base, zone_index, iconn, center, angle, translation) == CG_OK) {
      periodic.emplace ();
      periodic->name = name;
      periodic->zone = zone_number;
    }
    t8_debugf ("CGNS reader: Reading %s1to1 connectivity %s from zone %s to zone %s.\n", periodic ? "periodic " : "",
               name, zone.name.c_str (), donor_name);
    /* Iterate over the box of the range. */
    std::array<cgsize_t, 3> low = { 1, 1, 1 }, high = { 1, 1, 1 };
    for (int idim = 0; idim < dim; ++idim) {
      low[idim] = std::min (range[idim], range[dim + idim]);
      high[idim] = std::max (range[idim], range[dim + idim]);
    }
    /* Collect the vertex pairs first, so that nothing is merged if the connectivity turns out to be invalid. */
    std::vector<std::pair<t8_gloidx_t, t8_gloidx_t>> pairs;
    bool inside = true;
    for (cgsize_t k = low[2]; k <= high[2] && inside; ++k) {
      for (cgsize_t j = low[1]; j <= high[1] && inside; ++j) {
        for (cgsize_t i = low[0]; i <= high[0] && inside; ++i) {
          const cgsize_t index[3] = { i, j, k };
          cgsize_t donor_index[3] = { 1, 1, 1 };
          for (int idim = 0; idim < dim; ++idim) {
            const int donor_dim = std::abs (transform[idim]) - 1;
            const int sign = transform[idim] > 0 ? 1 : -1;
            donor_index[donor_dim] = sign * (index[idim] - range[idim]) + donor_range[donor_dim];
          }
          for (int idim = 0; idim < dim; ++idim) {
            inside = inside && donor_index[idim] >= 1 && donor_index[idim] <= donor_zone.vertex_dims[idim];
          }
          if (inside) {
            pairs.emplace_back (data.global_vertex (zone_number, zone.structured_vertex_id (index)),
                                data.global_vertex (donor_zone_number, donor_zone.structured_vertex_id (donor_index)));
          }
        }
      }
    }
    if (!inside) {
      t8_global_errorf (
        "CGNS reader: Warning: Skipping connectivity %s of zone %s. It maps outside of donor zone %s.\n", name,
        zone.name.c_str (), donor_name);
      continue;
    }
    for (const auto &[vertex, donor_vertex] : pairs) {
      data.merge_vertices (vertex, donor_vertex, periodic ? &*periodic : nullptr);
    }
    data.add_periodic_interface (periodic);
  }
  return true;
}

/** Compute the merged vertex ids of a tree face in t8code face vertex order. */
t8_cgns_face_entry
t8_cgns_tree_face (const t8_gloidx_t tree_id, const int face, t8_cgns_file_data &data)
{
  const t8_cgns_tree &tree = data.trees[tree_id];
  const t8_eclass_t face_class = static_cast<t8_eclass_t> (t8_eclass_face_types[tree.eclass][face]);
  t8_cgns_face_entry entry;
  entry.tree = tree_id;
  entry.face = face;
  entry.num_vertices = t8_eclass_num_vertices[face_class];
  entry.vertices.fill (-1);
  for (int ivertex = 0; ivertex < entry.num_vertices; ++ivertex) {
    const int tree_vertex = t8_face_vertex_to_tree_vertex[tree.eclass][face][ivertex];
    entry.vertices[ivertex] = data.merged_ids->find (tree.vertices[tree_vertex]);
  }
  return entry;
}

/** Insert all tree faces into the face hash. */
void
t8_cgns_build_face_map (t8_cgns_file_data &data)
{
  size_t num_degenerate = 0;
  for (size_t itree = 0; itree < data.trees.size (); ++itree) {
    const t8_cgns_tree &tree = data.trees[itree];
    for (int face = 0; face < t8_eclass_num_faces[tree.eclass]; ++face) {
      const t8_cgns_face_entry entry = t8_cgns_tree_face (static_cast<t8_gloidx_t> (itree), face, data);
      const t8_cgns_face_key key = t8_cgns_make_face_key (entry.vertices);
      if (t8_cgns_face_key_is_degenerate (key, entry.num_vertices)) {
        /* Two vertices of the face were merged into one, e.g. by a periodic interface of a single cell. */
        ++num_degenerate;
        continue;
      }
      data.faces[key].push_back (entry);
    }
  }
  if (num_degenerate > 0) {
    t8_global_errorf ("CGNS reader: Warning: %zu tree faces have coinciding vertices after merging the interfaces. "
                      "They are treated as domain boundaries.\n",
                      num_degenerate);
  }
}

/** Compute the orientation of a face connection, see t8_cmesh_set_join.
 * The orientation is the index of the first vertex of the smaller face within the vertices of the bigger face.
 * Follows the convention of the gmsh reader. */
int
t8_cgns_face_orientation (const t8_cgns_face_entry &face_a, const t8_cgns_face_entry &face_b,
                          const t8_cgns_file_data &data)
{
  const t8_eclass_t class_a = data.trees[face_a.tree].eclass;
  const t8_eclass_t class_b = data.trees[face_b.tree].eclass;
  const int compare = t8_eclass_compare (class_a, class_b);
  const t8_cgns_face_entry *smaller = nullptr;
  const t8_cgns_face_entry *bigger = nullptr;
  if (compare > 0 || (compare == 0 && face_a.face >= face_b.face)) {
    smaller = &face_b;
    bigger = &face_a;
  }
  else {
    smaller = &face_a;
    bigger = &face_b;
  }
  const t8_gloidx_t vertex_zero = smaller->vertices[0];
  for (int ivertex = 0; ivertex < bigger->num_vertices; ++ivertex) {
    if (bigger->vertices[ivertex] == vertex_zero) {
      return ivertex;
    }
  }
  SC_ABORT_NOT_REACHED ();
  return -1;
}

/** Compute the orientation of a periodic face connection, see \ref t8_cgns_face_orientation.
 * The vertices of \a face_a are mapped to the vertices of \a face_b by \a vertex_map before comparison.
 * \return The orientation or -1 if the vertices do not correspond.
 */
int
t8_cgns_periodic_face_orientation (const t8_cgns_face_entry &face_a, const t8_cgns_face_entry &face_b,
                                   const std::unordered_map<t8_gloidx_t, t8_gloidx_t> &vertex_map,
                                   const t8_cgns_file_data &data)
{
  const t8_eclass_t class_a = data.trees[face_a.tree].eclass;
  const t8_eclass_t class_b = data.trees[face_b.tree].eclass;
  const int compare = t8_eclass_compare (class_a, class_b);
  const bool a_is_smaller = !(compare > 0 || (compare == 0 && face_a.face >= face_b.face));
  if (a_is_smaller) {
    /* Find the image of the first vertex of face a among the vertices of face b. */
    const t8_gloidx_t mapped = vertex_map.at (face_a.vertices[0]);
    for (int ivertex = 0; ivertex < face_b.num_vertices; ++ivertex) {
      if (face_b.vertices[ivertex] == mapped) {
        return ivertex;
      }
    }
  }
  else {
    /* Find the vertex of face a whose image is the first vertex of face b. */
    for (int ivertex = 0; ivertex < face_a.num_vertices; ++ivertex) {
      if (vertex_map.at (face_a.vertices[ivertex]) == face_b.vertices[0]) {
        return ivertex;
      }
    }
  }
  return -1;
}

/** Join the tree faces on both sides of the periodic interfaces.
 * A boundary face whose vertices are all part of the interface is joined with the boundary face
 * formed by the corresponding donor vertices.
 * \param [in,out] joined  Faces that were already joined. Faces listed here are skipped and joined faces are added.
 */
void
t8_cgns_set_periodic_joins (t8_cmesh_t cmesh, t8_cgns_file_data &data, std::unordered_set<t8_gloidx_t> &joined)
{
  /* Faces are identified by tree id * T8_ECLASS_MAX_FACES + face in \a joined. */
  const auto face_id = [] (const t8_cgns_face_entry &entry) { return entry.tree * T8_ECLASS_MAX_FACES + entry.face; };
  for (const t8_cgns_periodic_interface &interface : data.periodic_interfaces) {
    std::unordered_map<t8_gloidx_t, t8_gloidx_t> vertex_map;
    for (const auto &[vertex, donor_vertex] : interface.vertex_pairs) {
      vertex_map[data.merged_ids->find (vertex)] = data.merged_ids->find (donor_vertex);
    }
    size_t num_joins = 0;
    size_t num_failed = 0;
    const t8_cgns_zone_data &zone = data.zones[interface.zone];
    for (t8_gloidx_t itree = zone.first_tree; itree < zone.first_tree + zone.num_trees; ++itree) {
      const t8_cgns_tree &tree = data.trees[itree];
      for (int face = 0; face < t8_eclass_num_faces[tree.eclass]; ++face) {
        const t8_cgns_face_entry entry = t8_cgns_tree_face (itree, face, data);
        if (joined.count (face_id (entry)) > 0) {
          continue;
        }
        const auto found = data.faces.find (t8_cgns_make_face_key (entry.vertices));
        if (found == data.faces.end () || found->second.size () != 1) {
          /* Not a boundary face. */
          continue;
        }
        /* Map the vertices to the donor side. */
        t8_cgns_face_key donor_vertices;
        donor_vertices.fill (-1);
        bool all_mapped = true;
        for (int ivertex = 0; ivertex < entry.num_vertices && all_mapped; ++ivertex) {
          const auto mapped = vertex_map.find (entry.vertices[ivertex]);
          all_mapped = mapped != vertex_map.end ();
          if (all_mapped) {
            donor_vertices[ivertex] = mapped->second;
          }
        }
        if (!all_mapped) {
          continue;
        }
        const auto donor_found = data.faces.find (t8_cgns_make_face_key (donor_vertices));
        if (donor_found == data.faces.end () || donor_found->second.size () != 1) {
          ++num_failed;
          continue;
        }
        const t8_cgns_face_entry &donor_entry = donor_found->second[0];
        if (face_id (donor_entry) == face_id (entry) || joined.count (face_id (donor_entry)) > 0) {
          continue;
        }
        const int orientation = t8_cgns_periodic_face_orientation (entry, donor_entry, vertex_map, data);
        if (orientation < 0) {
          ++num_failed;
          continue;
        }
        t8_cmesh_set_join (cmesh, entry.tree, donor_entry.tree, entry.face, donor_entry.face, orientation);
        joined.insert (face_id (entry));
        joined.insert (face_id (donor_entry));
        ++num_joins;
      }
    }
    if (num_failed > 0) {
      t8_global_errorf ("CGNS reader: Warning: %zu faces of periodic interface %s of zone %s have no matching donor "
                        "face. They are treated as domain boundaries.\n",
                        num_failed, interface.name.c_str (), zone.name.c_str ());
    }
    t8_debugf ("CGNS reader: Joined %zu pairs of tree faces across periodic interface %s of zone %s.\n", num_joins,
               interface.name.c_str (), zone.name.c_str ());
  }
}

/** Join all pairs of tree faces with the same key and the faces of the periodic interfaces. */
void
t8_cgns_set_joins (t8_cmesh_t cmesh, t8_cgns_file_data &data)
{
  size_t num_joins = 0;
  size_t num_nonconforming = 0;
  std::unordered_set<t8_gloidx_t> joined;
  for (const auto &[key, entries] : data.faces) {
    if (entries.size () == 2) {
      const int orientation = t8_cgns_face_orientation (entries[0], entries[1], data);
      t8_cmesh_set_join (cmesh, entries[0].tree, entries[1].tree, entries[0].face, entries[1].face, orientation);
      ++num_joins;
    }
    else if (entries.size () > 2) {
      ++num_nonconforming;
    }
  }
  if (num_nonconforming > 0) {
    t8_global_errorf ("CGNS reader: Warning: %zu faces are shared by more than two trees. They are treated as domain "
                      "boundaries.\n",
                      num_nonconforming);
  }
  t8_debugf ("CGNS reader: Joined %zu pairs of tree faces.\n", num_joins);
  t8_cgns_set_periodic_joins (cmesh, data, joined);
}

/** Label a tree face with a boundary condition. Keeps the first label if the face already has one.
 * \return True if the label was set.
 */
bool
t8_cgns_label_face (const t8_cgns_face_entry &entry, const std::string &label, t8_cgns_file_data &data)
{
  std::string &current = data.face_labels[entry.tree][entry.face];
  if (current.empty ()) {
    current = label;
    return true;
  }
  if (current != label) {
    t8_global_errorf ("CGNS reader: Warning: Face %d of tree %lld is part of the boundary conditions %s and %s. "
                      "Keeping %s.\n",
                      entry.face, static_cast<long long> (entry.tree), current.c_str (), label.c_str (),
                      current.c_str ());
  }
  return false;
}

/** Read the boundary conditions of a zone and label the tree faces. */
bool
t8_cgns_read_boundary_conditions (const int fn, const int base, const int zone_index, t8_cgns_file_data &data,
                                  std::vector<t8_cmesh_cgns_boundary_condition> &boundary_conditions)
{
  const int zone_number = zone_index - 1;
  const t8_cgns_zone_data &zone = data.zones[zone_number];
  int num_bcs = 0;
  T8_CGNS_CHECK (cg_nbocos (fn, base, zone_index, &num_bcs));
  for (int ibc = 1; ibc <= num_bcs; ++ibc) {
    char name[T8_CGNS_NAME_LENGTH];
    CGNS_ENUMT (BCType_t) bc_type;
    CGNS_ENUMT (PointSetType_t) point_set_type;
    CGNS_ENUMT (DataType_t) normal_data_type;
    cgsize_t num_points = 0;
    cgsize_t normal_list_size = 0;
    int normal_index[3] = { 0, 0, 0 };
    int num_datasets = 0;
    T8_CGNS_CHECK (cg_boco_info (fn, base, zone_index, ibc, name, &bc_type, &point_set_type, &num_points, normal_index,
                                 &normal_list_size, &normal_data_type, &num_datasets));
    CGNS_ENUMT (GridLocation_t) location;
    T8_CGNS_CHECK (cg_boco_gridlocation_read (fn, base, zone_index, ibc, &location));
    std::vector<cgsize_t> raw_points (std::max<cgsize_t> (num_points * zone.index_dim, 1));
    /* The normal list is not needed, but the buffer must be large enough for the library to write into. */
    std::vector<double> normals (3 * static_cast<size_t> (std::max<cgsize_t> (normal_list_size, 1)));
    T8_CGNS_CHECK (cg_boco_read (fn, base, zone_index, ibc, raw_points.data (), normals.data ()));

    t8_cmesh_cgns_boundary_condition bc;
    bc.name = name;
    bc.type = static_cast<int> (bc_type);
    bc.type_name = cg_BCTypeName (bc_type);
    bc.zone = zone_number;
    bc.num_faces = 0;
    if (cg_goto (fn, base, "Zone_t", zone_index, "ZoneBC_t", 1, "BC_t", ibc, "end") == CG_OK) {
      char family_name[T8_CGNS_NAME_LENGTH];
      if (cg_famname_read (family_name) == CG_OK) {
        bc.family = family_name;
      }
    }
    data.has_boundary_conditions = true;

    if (location == CGNS_ENUMV (Vertex)) {
      /* The points are vertex ids. A boundary face belongs to the BC if all its vertices are listed. */
      const std::vector<cgsize_t> points = zone.point_set_to_ids (point_set_type, raw_points, num_points);
      std::unordered_set<t8_gloidx_t> bc_vertices;
      for (const cgsize_t id : points) {
        bc_vertices.insert (data.merged_ids->find (data.global_vertex (zone_number, id)));
      }
      for (t8_gloidx_t itree = zone.first_tree; itree < zone.first_tree + zone.num_trees; ++itree) {
        const t8_cgns_tree &tree = data.trees[itree];
        for (int face = 0; face < t8_eclass_num_faces[tree.eclass]; ++face) {
          const t8_cgns_face_entry entry = t8_cgns_tree_face (itree, face, data);
          const auto found = data.faces.find (t8_cgns_make_face_key (entry.vertices));
          if (found == data.faces.end () || found->second.size () != 1) {
            /* Not a boundary face. */
            continue;
          }
          const bool all_listed = std::all_of (entry.vertices.begin (), entry.vertices.begin () + entry.num_vertices,
                                               [&] (const t8_gloidx_t id) { return bc_vertices.count (id) > 0; });
          if (all_listed && t8_cgns_label_face (entry, bc.name, data)) {
            ++bc.num_faces;
          }
        }
      }
    }
    else if (zone.structured) {
      t8_global_errorf ("CGNS reader: Warning: Skipping boundary condition %s of structured zone %s. GridLocation %s "
                        "is not supported for structured zones, use Vertex.\n",
                        name, zone.name.c_str (), cg_GridLocationName (location));
    }
    else {
      /* The points are ids of face elements. */
      const std::vector<cgsize_t> points = zone.point_set_to_ids (point_set_type, raw_points, num_points);
      size_t num_skipped = 0;
      for (const cgsize_t element_id : points) {
        const std::optional<t8_cgns_element> element = zone.get_element (element_id);
        if (!element || element->eclass == T8_ECLASS_COUNT
            || t8_eclass_to_dimension[element->eclass] != data.cell_dim - 1) {
          ++num_skipped;
          continue;
        }
        t8_cgns_face_key vertices;
        vertices.fill (-1);
        const int num_vertices = static_cast<int> (element->vertices.size ());
        for (int ivertex = 0; ivertex < num_vertices; ++ivertex) {
          vertices[ivertex] = data.merged_ids->find (data.global_vertex (zone_number, element->vertices[ivertex]));
        }
        const auto found = data.faces.find (t8_cgns_make_face_key (vertices));
        if (found == data.faces.end ()) {
          ++num_skipped;
          continue;
        }
        for (const t8_cgns_face_entry &entry : found->second) {
          /* Only label faces of trees of this zone. The face hash may also contain the face of the
           * neighboring tree, if the boundary condition lies on an interface. */
          if (data.trees[entry.tree].zone == zone_number && t8_cgns_label_face (entry, bc.name, data)) {
            ++bc.num_faces;
          }
        }
      }
      if (num_skipped > 0) {
        t8_global_errorf ("CGNS reader: Warning: %zu of %zu face elements of boundary condition %s of zone %s do not "
                          "match a tree face and were skipped.\n",
                          num_skipped, points.size (), name, zone.name.c_str ());
      }
    }
    t8_debugf ("CGNS reader: Boundary condition %s (%s, family \"%s\") of zone %s covers %lld tree faces.\n", name,
               bc.type_name.c_str (), bc.family.c_str (), zone.name.c_str (), static_cast<long long> (bc.num_faces));
    boundary_conditions.push_back (std::move (bc));
  }
  return true;
}

/** Build the cmesh from the collected data. */
t8_cmesh_t
t8_cgns_build_cmesh (t8_cgns_file_data &data, sc_MPI_Comm comm)
{
  t8_cmesh_t cmesh;
  t8_cmesh_init (&cmesh);
  t8_cmesh_set_dimension (cmesh, data.cell_dim);
  const t8_geometry_c *geometry = t8_cmesh_register_geometry<t8_geometry_linear> (cmesh);

  std::array<double, 3 * T8_ECLASS_MAX_CORNERS> vertex_coordinates;
  std::array<t8_gloidx_t, T8_ECLASS_MAX_CORNERS> global_vertices;
  for (size_t itree = 0; itree < data.trees.size (); ++itree) {
    const t8_gloidx_t tree_id = static_cast<t8_gloidx_t> (itree);
    const t8_cgns_tree &tree = data.trees[itree];
    const int num_vertices = t8_eclass_num_vertices[tree.eclass];
    for (int ivertex = 0; ivertex < num_vertices; ++ivertex) {
      const std::array<double, 3> coordinates = data.vertex_coordinates (tree.vertices[ivertex]);
      std::copy (coordinates.begin (), coordinates.end (), vertex_coordinates.begin () + 3 * ivertex);
      global_vertices[ivertex] = data.merged_ids->find (tree.vertices[ivertex]);
    }
    t8_cmesh_set_tree_class (cmesh, tree_id, tree.eclass);
    t8_cmesh_set_tree_vertices (cmesh, tree_id, vertex_coordinates.data (), num_vertices);
    t8_cmesh_set_global_vertices_of_tree (cmesh, tree_id, global_vertices.data (), num_vertices);
    t8_cmesh_set_tree_geometry (cmesh, tree_id, geometry);
    int zone = tree.zone;
    t8_cmesh_set_attribute (cmesh, tree_id, t8_get_package_id (), T8_CMESH_CGNS_ZONE_ATTRIBUTE_KEY, &zone, sizeof (int),
                            0);
    if (data.has_boundary_conditions) {
      const int num_faces = t8_eclass_num_faces[tree.eclass];
      const std::vector<std::string> labels (data.face_labels[itree].begin (),
                                             data.face_labels[itree].begin () + num_faces);
      t8_cmesh_set_boundary_conditions (cmesh, tree_id, labels);
    }
  }
  t8_cgns_set_joins (cmesh, data);
  t8_cmesh_commit (cmesh, comm);
  return cmesh;
}

} /* namespace */

t8_cmesh_cgns_reader::t8_cmesh_cgns_reader (std::string filename): m_filename (std::move (filename))
{
}

t8_cmesh_t
t8_cmesh_cgns_reader::read (sc_MPI_Comm comm)
{
  m_zones.clear ();
  m_boundary_conditions.clear ();
  m_dimension = -1;

  int fn = -1;
  if (cg_open (m_filename.c_str (), CG_MODE_READ, &fn) != CG_OK) {
    t8_global_errorf ("CGNS reader: Could not open file %s: %s\n", m_filename.c_str (), cg_get_error ());
    return nullptr;
  }
  t8_cgns_file_data data;
  const int base = 1;
  const bool success = [&] () -> bool {
    int num_bases = 0;
    T8_CGNS_CHECK (cg_nbases (fn, &num_bases));
    if (num_bases < 1) {
      t8_global_errorf ("CGNS reader: File %s does not contain a base.\n", m_filename.c_str ());
      return false;
    }
    if (num_bases > 1) {
      t8_global_errorf ("CGNS reader: Warning: File %s contains %d bases. Only the first one is read.\n",
                        m_filename.c_str (), num_bases);
    }
    if (!t8_cgns_read_zones (fn, base, data)) {
      return false;
    }
    data.merged_ids.emplace (data.num_global_vertices ());
    for (size_t izone = 0; izone < data.zones.size (); ++izone) {
      const int zone_index = static_cast<int> (izone) + 1;
      if (!t8_cgns_read_connectivities (fn, base, zone_index, data)
          || !t8_cgns_read_1to1_connectivities (fn, base, zone_index, data)) {
        return false;
      }
    }
    t8_cgns_build_face_map (data);
    data.face_labels.resize (data.trees.size ());
    for (size_t izone = 0; izone < data.zones.size (); ++izone) {
      if (!t8_cgns_read_boundary_conditions (fn, base, static_cast<int> (izone) + 1, data, m_boundary_conditions)) {
        return false;
      }
    }
    return true;
  }();
  cg_close (fn);
  if (!success) {
    m_boundary_conditions.clear ();
    return nullptr;
  }

  for (const t8_cgns_zone_data &zone : data.zones) {
    m_zones.push_back ({ zone.name, zone.first_tree, zone.num_trees, static_cast<t8_gloidx_t> (zone.num_vertices) });
  }
  m_dimension = data.cell_dim;
  t8_global_productionf ("CGNS reader: Read %zu trees in %zu zones with %zu boundary conditions from %s.\n",
                         data.trees.size (), data.zones.size (), m_boundary_conditions.size (), m_filename.c_str ());
  return t8_cgns_build_cmesh (data, comm);
}

T8_EXTERN_C_BEGIN ();

t8_cmesh_t
t8_cmesh_from_cgns_file (const char *filename, sc_MPI_Comm comm)
{
  t8_cmesh_cgns_reader reader (filename);
  return reader.read (comm);
}

int
t8_cmesh_get_cgns_zone (t8_cmesh_t cmesh, t8_locidx_t ltreeid)
{
  T8_ASSERT (t8_cmesh_is_committed (cmesh, 0));
  const int *zone = static_cast<const int *> (
    t8_cmesh_get_attribute (cmesh, t8_get_package_id (), T8_CMESH_CGNS_ZONE_ATTRIBUTE_KEY, ltreeid));
  return zone == nullptr ? -1 : *zone;
}

T8_EXTERN_C_END ();
