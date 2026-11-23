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

/** \file t8_cmesh_cgns_reader.hxx
 * Reads a CGNS file (https://cgns.github.io) and constructs a cmesh from it.
 *
 * The reader is only available if t8code was configured with T8CODE_ENABLE_CGNS.
 *
 * Supported input:
 *  - The first base of the file with cell dimension 1, 2 or 3 and Cartesian coordinates.
 *  - Unstructured zones with element sections of the linear types NODE, BAR_2, TRI_3, QUAD_4, TETRA_4, PYRA_5,
 *    PENTA_6 and HEXA_8, their higher order variants (only the corner nodes are used) and MIXED sections of
 *    these types. Elements whose dimension equals the cell dimension of the base become the trees of the cmesh.
 *    Elements of one dimension less are used to resolve boundary conditions and zone interfaces.
 *    NGON_n and NFACE_n sections are skipped.
 *  - Structured zones. Every cell of the index grid becomes a line, quad or hex tree. The cells are numbered
 *    with the first index running fastest.
 *  - Zones: The trees of a zone get consecutive global tree ids in the order of the zones in the file.
 *    The index of the zone is stored at every tree as an attribute, see \ref t8_cmesh_get_cgns_zone.
 *  - Neighborships: Trees of one zone are connected via shared vertices. Interfaces between zones and
 *    periodic interfaces are read from the ZoneGridConnectivity nodes:
 *      - GridConnectivity1to1_t nodes of structured zones (index ranges with Transform).
 *      - GridConnectivity_t nodes of type Abutting1to1 with GridLocation Vertex (vertex to vertex
 *        correspondence, structured and unstructured zones) or FaceCenter/EdgeCenter (face element to face
 *        element correspondence, unstructured zones only, matched geometrically).
 *    Periodic transformations (RotationCenter, RotationAngle, Translation) are taken into account.
 *  - Boundary conditions: Every BC_t node of a zone (ZoneBC) is converted to a boundary condition label of the
 *    tree faces it covers. Point sets with GridLocation Vertex (all zone types) and FaceCenter/EdgeCenter
 *    (unstructured zones) are supported. The label is the name of the BC_t node. The labels are stored via
 *    \ref t8_cmesh_set_boundary_conditions and can be queried with \ref t8_cmesh_get_boundary_condition.
 *    Faces without a boundary condition in the file get the empty label "".
 *    The BC type and the family name are available via \ref t8_cmesh_cgns_reader::get_boundary_conditions.
 *
 * The cmesh is replicated, i.e. every process reads the complete file. Use \ref t8_cmesh_set_derive and
 * \ref t8_cmesh_set_partition_uniform afterwards to partition it.
 */

#ifndef T8_CMESH_CGNS_READER_HXX
#define T8_CMESH_CGNS_READER_HXX

#include <t8.h>
#include <t8_cmesh/t8_cmesh.h>

#include <string>
#include <vector>

/** Information about one zone of the CGNS file. */
struct t8_cmesh_cgns_zone
{
  std::string name;         /**< Name of the Zone_t node. */
  t8_gloidx_t first_tree;   /**< Global id of the first tree that was created from this zone. */
  t8_gloidx_t num_trees;    /**< Number of trees that were created from this zone. */
  t8_gloidx_t num_vertices; /**< Number of vertices of the zone in the file. */
};

/** Information about one boundary condition (BC_t node) of the CGNS file. */
struct t8_cmesh_cgns_boundary_condition
{
  std::string name;      /**< Name of the BC_t node. This is the label that is attached to the tree faces. */
  std::string family;    /**< FamilyName of the boundary condition, empty if none is given. */
  std::string type_name; /**< Name of the CGNS BCType_t, for example "BCWall". */
  int type;              /**< The CGNS BCType_t enum value. */
  int zone;              /**< Index of the zone (into \ref t8_cmesh_cgns_reader::get_zones) the BC belongs to. */
  t8_gloidx_t num_faces; /**< Number of tree faces that were labeled with this boundary condition. */
};

/** Reads a CGNS file and constructs a committed cmesh from it.
 * After \ref read was called successfully, the zones and boundary conditions found in the file
 * can be queried.
 *
 * \code
 * t8_cmesh_cgns_reader reader ("mesh.cgns");
 * t8_cmesh_t cmesh = reader.read (sc_MPI_COMM_WORLD);
 * for (const auto &zone : reader.get_zones ()) { ... }
 * \endcode
 */
class t8_cmesh_cgns_reader {
 public:
  /** Constructor.
   * \param [in] filename  Path of the CGNS file to read. The file is not opened before \ref read is called.
   */
  explicit t8_cmesh_cgns_reader (std::string filename);

  /** Read the file and construct a committed cmesh.
   * This function is collective over \a comm. Every process reads the complete file.
   * \param [in] comm  The MPI communicator to commit the cmesh with.
   * \return           The committed cmesh or nullptr if the file could not be read.
   */
  t8_cmesh_t
  read (sc_MPI_Comm comm);

  /** \return The zones of the file in the order of their tree ids. Only valid after a successful \ref read. */
  const std::vector<t8_cmesh_cgns_zone> &
  get_zones () const
  {
    return m_zones;
  }

  /** \return The boundary conditions of the file. Only valid after a successful \ref read. */
  const std::vector<t8_cmesh_cgns_boundary_condition> &
  get_boundary_conditions () const
  {
    return m_boundary_conditions;
  }

  /** \return The cell dimension of the base, i.e. the dimension of the cmesh. Only valid after a successful \ref read. */
  int
  get_dimension () const
  {
    return m_dimension;
  }

  /** \return The path of the file. */
  const std::string &
  get_filename () const
  {
    return m_filename;
  }

 private:
  std::string m_filename;                                              /**< Path of the CGNS file. */
  std::vector<t8_cmesh_cgns_zone> m_zones;                             /**< Zones found in the file. */
  std::vector<t8_cmesh_cgns_boundary_condition> m_boundary_conditions; /**< Boundary conditions found in the file. */
  int m_dimension = -1;                                                /**< Cell dimension of the base. */
};

T8_EXTERN_C_BEGIN ();

/** Read a CGNS file and construct a committed cmesh from it.
 * Convenience wrapper around \ref t8_cmesh_cgns_reader::read. The zone names and boundary condition
 * details of the file are not accessible via this function, use the class for that.
 * \param [in] filename  Path of the CGNS file.
 * \param [in] comm      The MPI communicator to commit the cmesh with.
 * \return               The committed cmesh or NULL if the file could not be read.
 */
t8_cmesh_t
t8_cmesh_from_cgns_file (const char *filename, sc_MPI_Comm comm);

/** Get the index of the CGNS zone a tree was read from.
 * \param [in] cmesh    A committed cmesh that was read with the CGNS reader.
 * \param [in] ltreeid  The local id of a tree of \a cmesh.
 * \return              The index of the zone (0-based, in the order of \ref t8_cmesh_cgns_reader::get_zones)
 *                      or -1 if the tree has no zone information.
 */
int
t8_cmesh_get_cgns_zone (t8_cmesh_t cmesh, t8_locidx_t ltreeid);

T8_EXTERN_C_END ();

#endif /* !T8_CMESH_CGNS_READER_HXX */
